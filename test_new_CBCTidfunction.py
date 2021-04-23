"""
date : tuesday 12 may 2020
author: Claudia Acquistapace
goal: build a function to identify all cloud base and cloud top of cloud entities occurring at the same time.
The idea is also to include human observations from the day to distinguish PBL from non-PBL clouds, but including all of
them in the clouds dataset, with and index identifying the cloud layer.
Concept of the code:
1) given the cloud mask, run the code to identify for each time, all cloud bases and all cloud tops using gradient method.
2) Identify cloud layers for each time by, given the cloud base value, find the min ( CT >= CBi). Assign an index for each
cloud layer identified in this way.
3) discriminate layers belonging to PBL from layers not belonging to PBL. Using the max cloud top height introduced by human
observations for the case study selected

"""
import xarray as xr
import numpy as np
import matplotlib
import scipy
import pylab
import netCDF4 as nc4
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import struct
import glob
import pandas as pd
import datetime as dt
import random
import datetime
import matplotlib.dates as mdates
import os
import atmos
import matplotlib as mpl
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
from myFunctions import f_closest
from myFunctions import f_resampleArrays2StandardData
from myFunctions import f_resample2StandardData
from myFunctions import f_plotTest

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']

def f_calculateCloudBaseTopThickness(cloudMask, time, height, humanInfo):
    """
    date : wednesday 13 may 2020
    author: Claudia Acquistapace
    goal: build a function to identify all cloud base and cloud top of clouds in the vertical profile at the same time.
    Human observations for the day distinguish manually PBL from non-PBL clouds. An additional dataset of
     PBL clouds is delivered based on this information.
    Concept of the code:
    step 1: given the cloud mask, find all cloud base and cloud tops.
    step 2: build cloud database with all clouds saved as xarray dataset
    step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
    input: cloudmask,
            time,
            height,
            humanInfo (dictionary including timeStart, timeEnd, PBLheight from human obs on the day)
    output: AllCloudDataset (xarray Dataset including cloud base, cloud top, cloud thickness, level number)
            PBLcloudDataset (xarray Dataset for PBL clouds with cloud base, cloud top, cloud thickness, level number)
    """
    dimTime   = len(time)
    dimHeight = len(height)
    heightPBL = humanInfo['heightPBL']
    timeStart = humanInfo['timeStart']
    timeEnd   = humanInfo['timeEnd']

    # STEP 1: identifying all cloud bases and tops
    # ---------------------------------------------------
    # converting cloud mask to 1 / 0 matrices
    BinaryMatrix = np.zeros((dimTime, dimHeight))
    for itime in range(dimTime):
        for iH in range(dimHeight):
            if cloudMask[itime, iH] != 0.:
                BinaryMatrix[itime, iH] = 1

    # calculating gradient of binary cloud mask
    gradBinary = np.diff(BinaryMatrix, axis=1)

    # counting max number of cloud base/cloud top found
    numberCB = []
    numberCT = []
    for itime in range(dimTime):
        column = gradBinary[itime, :]
        numberCB.append(len(np.where(column == -1.)[0][:]))
        numberCT.append(len(np.where(column == 1.)[0][:]))

    NCB = max(numberCB)
    NCT = max(numberCT)

    # generating cloud base and cloud top arrays
    CBarray = np.zeros((dimTime, NCB))
    CBarray.fill(np.nan)
    CTarray = np.zeros((dimTime, NCT))
    CTarray.fill(np.nan)
    NlayersArray = np.zeros((dimTime))
    NlayersArray.fill(np.nan)

    # if no cloud bases or no cloud tops are found, then CB and CT are assigned to nan
    if (NCB == 0) or (NCT == 0):
        CBarray[iTime, :] = np.nan
        CTarray[iTime, :] = np.nan
    else:
        # if some cloud base / cloud tops are found, all the found values are stored
        # storing cloud base and cloud top arrays
        for iTime in range(dimTime):
            column = gradBinary[iTime, :]
            indCB = np.where(column == -1.)[0][:]
            NfoundCB = len(indCB)
            indCT = np.where(column == 1.)[0][:]
            NfoundCT = len(indCT)
            CBarray[iTime, 0:NfoundCB] = height[indCB]
            CTarray[iTime, 0:NfoundCT] = height[indCT]
            NlayersArray[iTime] = numberCB[iTime]

    # calculating cloud thickness based on the cloud base and tops found ( 2d array (time, Nlevels))
    cloudThicknessDatabase = CTarray - CBarray
    # generating array of levels
    levels = np.arange(NCB)

    # step 2: build cloud database with all clouds saved as xarray dataset
    clouds = xr.Dataset(
        data_vars = {'cloudBase' : (('time', 'levels'), CBarray),
                     'cloudTop'  : (('time', 'levels'), CTarray),
                     'cloudThick': (('time', 'levels'), cloudThicknessDatabase)},
        coords    = {'levels': levels,
                     'time'  : time})

    # step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
    cloudsTimeWindow = clouds.sel(time=slice(timeStart, timeEnd), drop=False)
    #PBLclouds = cloudsTimeWindow.where(cloudsTimeWindow.cloudTop < heightPBL, drop=True)
    PBLclouds = cloudsTimeWindow.where(cloudsTimeWindow.cloudTop < heightPBL, drop=False)


    return(clouds, PBLclouds)

print('sto leggendo il codice')
# setting parameters for calculating averaging and domain size of the model:
NprofilesOut                    = 24  # hourly means
timeIncrement                   = 60  # hourly means
patch                           = 'patch003'
flagPlot                        = 0

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/old_before_13052020/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/old_before_13052020/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/debugging/newFunction/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)
date_arr                        = []
cloudMask_obs                   = np.zeros((Nfiles,9600, 150))
cloudMask_mod                   = np.zeros((Nfiles,9600, 150))

# days selected to test the algorithm
DaysArr                         = ['20130427']#,, '',]'20130428'
for iDayTest in range(len(DaysArr)):
    for indFile in range(Nfiles):
        #date = fileListMod[indFile][87:95]
        date = fileListMod[indFile][107:115]
        print(date)
        if date == DaysArr[iDayTest]:
            print('Date matched, processing')
            filenameMod = fileListMod[indFile]
            filenameObs = fileListObs[indFile]
            #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
            #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'
            #print(filenameMod)
            #print(filenameObs)
            yy = int(date[0:4])
            mm = int(date[4:6])
            dd = int(date[6:8])
            timeStart = datetime.datetime(yy, mm, dd, 6)
            timeEnd = datetime.datetime(yy, mm, dd, 23)
            heightPBL = 2500.#2500.

            humanInfo = {'timeStart':timeStart, 'timeEnd':timeEnd, 'heightPBL':heightPBL}
            date_arr.append(date)
            print('processing date ' + date)

            # reading time and height from ncdf file (grid of ICON LEM
            # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
            ncdata = Dataset(filenameMod, mode='r')
            time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                                ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
            height = ncdata.groups['Temp_data'].variables['height2'][:]

            # opening the file containing all the data
            infile = open(filenameObs, 'rb')
            new_dict = pickle.load(infile, encoding='latin1')

            cloudMask_obs[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[8]['cloudMask']), time[:], \
                                                                 new_dict[3]['height'], date)
            cloudMask_mod[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[7]['cloudMask']), time[:], \
                                                                    new_dict[3]['height'], date)
            cloudMaskModPlot = cloudMask_mod[indFile,:,:]
            cloudMaskObsPlot = cloudMask_obs[indFile,:,:]


            # call function to calculate CB, CT, thickness
            clouds, PBLclouds = f_calculateCloudBaseTopThickness(cloudMask_obs[indFile,:,:], time, height, humanInfo)

            # deriving lowest cloud base and corresponding cloud top for PBL clouds
            CBarr = np.zeros(len(time))
            CBarr.fill(np.nan)
            CTarr = np.zeros(len(time))
            CTarr.fill(np.nan)
            iPBL = 0
            for itime in range(len(time)):
                if iPBL < len(PBLclouds.time.values):
                    if clouds.time.values[itime] == PBLclouds.time.values[iPBL]:
                        print(iPBL)
                        CBarray = PBLclouds.cloudBase.values[iPBL, :]
                        if CBarray.size - np.count_nonzero(np.isnan(CBarray)) != 0:
                            minCB = np.nanmin(PBLclouds.cloudBase.values[iPBL, :])
                            CBarr[itime] = minCB
                            indexLevelMin = np.nanargmin(PBLclouds.cloudBase.values[iPBL, :])
                            CTarr[itime] = PBLclouds.cloudTop[iPBL, indexLevelMin]
                        else:
                            CBarr[itime] = np.nan
                            CTarr[itime] = np.nan
                        iPBL=iPBL+1
            print('cloud base and cloud top for ICON-LEM calculated ')


            # reading variabls to produce histograms
            geomThickness = clouds.cloudThick.values.flatten()
            geomThicknessPBL = PBLclouds.cloudThick.values.flatten()
            CBall = clouds.cloudBase.values.flatten()
            CBPBL = PBLclouds.cloudBase.values.flatten()
            CTall = clouds.cloudTop.values.flatten()
            CTPBL = PBLclouds.cloudTop.values.flatten()
            NlevelsPBL = len(PBLclouds.levels)
            NlevelsMax = len(clouds.levels)

            # plotting cloud base and cloud top selected
            fontSizeTitle = 12
            fontSizeX = 10
            fontSizeY = 10
            fs = 10
            fontSizeCbar = 10
            labelsizeaxes = 10
            cbarAspect = 10
            Ymax = 11000
            Ymin = 107.
            Nrows = 2
            Ncols = 1
            fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(10, 10))
            axes[0] = plt.subplot(Nrows, Ncols, 1)
            axes[0].spines["top"].set_visible(False)
            axes[0].spines["right"].set_visible(False)
            axes[0].get_xaxis().tick_bottom()
            axes[0].get_yaxis().tick_left()
            axes[0].grid(True, which="both")
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            axes[0].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
            colorArr = ['blue','green','orange','red','purple','pink']
            #for iplot in range(NlevelsMax):
            #    clouds.cloudBase[:,iplot].plot.line(marker='o', color=colorArr[iplot])
            #for iplot in range(NlevelsPBL):
            #    PBLclouds.cloudBase[:,iplot].plot.line(color='black', marker='o')
            plt.plot(time, CBarr, color='black')
            axes[1] = plt.subplot(Nrows, Ncols, 2)
            axes[1].spines["top"].set_visible(False)
            axes[1].spines["right"].set_visible(False)
            axes[1].get_xaxis().tick_bottom()
            axes[1].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            axes[1].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
            axes[1].grid(True, which="both")
            plt.plot(time, CTarr, color='black')

            #for iplot in range(NlevelsMax):
            #    clouds.cloudTop[:,iplot].plot.line(marker='o', color=colorArr[iplot])
            #for iplot in range(NlevelsPBL):
            #    PBLclouds.cloudTop[:,iplot].plot.line(color='black', marker='o')
            plt.tight_layout()
            plt.savefig(pathFig+date+'_timeSeries_CB_CT.png')

            # plotting distributions of CB, CT, cloud thickness for PBl and all clouds
            Nrows = 3
            Ncols = 2
            fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(15, 10))
            axes[0,0] = plt.subplot(Nrows, Ncols, 1)
            axes[0,0].spines["top"].set_visible(False)
            axes[0,0].spines["right"].set_visible(False)
            axes[0,0].get_xaxis().tick_bottom()
            axes[0,0].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[0,0].grid(True, which="both")
            axes[0,0].hist(CBall, range=[0., 8000.], bins=50, density=False, color='red', label='all clouds')
            plt.legend(frameon=False)
            axes[0,0].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[0,0].set_xlabel('occurrences [%]', fontsize=fs)

            axes[0,1] = plt.subplot(Nrows, Ncols, 2)
            axes[0,1].spines["top"].set_visible(False)
            axes[0,1].spines["right"].set_visible(False)
            axes[0,1].get_xaxis().tick_bottom()
            axes[0,1].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[0,1].grid(True, which="both")
            axes[0,1].hist(CBPBL, range=[0., 2500.], bins=20, density=False, color='black', label='PBL clouds')
            plt.legend(frameon=False)
            axes[0,1].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[0,1].set_xlabel('occurrences [%]', fontsize=fs)

            axes[1,0] = plt.subplot(Nrows, Ncols, 3)
            axes[1,0].spines["top"].set_visible(False)
            axes[1,0].spines["right"].set_visible(False)
            axes[1,0].get_xaxis().tick_bottom()
            axes[1,0].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[1,0].grid(True, which="both")
            axes[1,0].hist(CTall, range=[0., 10000.], bins=50, density=False, color='red', label='all clouds')
            plt.legend(frameon=False)
            axes[1,0].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[1,0].set_xlabel('occurrences [%]', fontsize=fs)

            axes[1,1] = plt.subplot(Nrows, Ncols, 4)
            axes[1,1].spines["top"].set_visible(False)
            axes[1,1].spines["right"].set_visible(False)
            axes[1,1].get_xaxis().tick_bottom()
            axes[1,1].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[1,1].grid(True, which="both")
            axes[1,1].hist(CTPBL, range=[0., 2500.], bins=20, density=False, color='black', label='PBL clouds')
            plt.legend(frameon=False)
            axes[1,1].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[1,1].set_xlabel('occurrences [%]', fontsize=fs)

            axes[2,0] = plt.subplot(Nrows, Ncols, 5)
            axes[2,0].spines["top"].set_visible(False)
            axes[2,0].spines["right"].set_visible(False)
            axes[2,0].get_xaxis().tick_bottom()
            axes[2,0].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[2,0].grid(True, which="both")
            axes[2,0].hist(geomThickness, range=[0., 10000.], bins=50, density=False, color='red', label='all clouds')
            plt.legend(frameon=False)
            axes[2,0].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[2,0].set_xlabel('occurrences [%]', fontsize=fs)

            axes[2,1] = plt.subplot(Nrows, Ncols, 6)
            axes[2,1].spines["top"].set_visible(False)
            axes[2,1].spines["right"].set_visible(False)
            axes[2,1].get_xaxis().tick_bottom()
            axes[2,1].get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            axes[2,1].grid(True, which="both")
            axes[2,1].hist(geomThicknessPBL, range=[0., 2000.], bins=20, density=False, color='black', label='PBL clouds')
            plt.legend(frameon=False)
            axes[2,1].set_ylabel('Cloud base [m]', fontsize=fs)
            axes[2,1].set_xlabel('occurrences [%]', fontsize=fs)
            plt.tight_layout()
            plt.savefig(pathFig+date+'_histograms_CB_CT_TK.png')




PlotFlag_cloudMask = 1
# plot of cloud fraction time height for all days of simulation/obs
if PlotFlag_cloudMask == 1:
    for indDay in range(len(date_arr)):
        date = date_arr[indDay]
        if date == DaysArr[iDayTest]:
            yy = int(date[0:4])
            mm = int(date[4:6])
            dd = int(date[6:8])
            print('plotting date :' + date)
            timePlotHalfHours = pd.date_range(date, periods=9600, freq='9s')
            timeStart = datetime.datetime(yy,mm,dd,6)
            timeEnd   = datetime.datetime(yy,mm,dd,23)
            Ncols = 1
            Nrows = 2
            Nplots = 2
            fontSizeTitle = 12
            fontSizeX = 10
            fontSizeY = 10
            fontSizeCbar = 10
            labelsizeaxes = 10
            cbarAspect = 10
            Ymax = 11000
            Ymin = 107.
            from matplotlib import rcParams
            rcParams['font.family'] = 'sans-serif'
            rcParams['font.sans-serif'] = ['Tahoma']
            fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(10,5))
            ax = plt.subplot(Nrows, Ncols, 1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
            cax = ax.pcolormesh(timePlotHalfHours, height, cloudMaskModPlot.transpose(), vmin=0, vmax=3,
                                cmap=plt.cm.get_cmap("GnBu", 4))
    
            ax.set_ylim(Ymin, Ymax)  # limits of the y-axes
            ax.set_xlim(timeStart,timeEnd)                                 # limits of the x-axes
            ax.set_title("cloud mask model", fontsize=fontSizeTitle)
            ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
            ax.set_ylabel("height [m]", fontsize=fontSizeY)
            cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical', aspect=cbarAspect)
            cbar.ticks = ([0, 1, 2, 3])
            cbar.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
            cbar.set_label(label="cloud type", size=fontSizeCbar)
            cbar.ax.tick_params(labelsize=fontSizeCbar)
    
            ax1 = plt.subplot(Nrows, Ncols, 2)
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_visible(False)
            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.set_xlim(timeStart,timeEnd)
            ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
            cax1 = ax1.pcolormesh(timePlotHalfHours, height, cloudMaskObsPlot.transpose(), vmin=0, vmax=3,
                                cmap=plt.cm.get_cmap("GnBu", 4))
            ax1.set_ylim(Ymin, Ymax)  # limits of the y-axes
            ax1.set_title("cloud mask obs", fontsize=fontSizeTitle)
            ax1.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
            ax1.set_ylabel("height [m]", fontsize=fontSizeY)
            cbar1 = fig.colorbar(cax1, ticks=[0, 1, 2, 3], orientation='vertical', aspect=cbarAspect)
            cbar1.ticks = ([0, 1, 2, 3])
            cbar1.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
            cbar1.set_label(label="cloud type", size=fontSizeCbar)
            cbar1.ax.tick_params(labelsize=fontSizeCbar)
            plt.tight_layout()
            plt.savefig(pathFig+date+'_cloudMask.png')


            ## step 1: find all cloud base and cloud tops.
            #cloudBAseTopAll = f_testCBCTId(cloudMask_obs[indFile,:,:], time, height)

            # step 2: build cloud database with all clouds saved as xarray dataset
            #cloudBaseDatabase = cloudBAseTopAll[0]
            #cloudTopDatabase  = cloudBAseTopAll[1]
            #NlevelsMax        = cloudBAseTopAll[2]
            #NlevelsArr        = cloudBAseTopAll[3]
            #gradBinary        = cloudBAseTopAll[4]
            #cloudThicknessDatabase    = cloudTopDatabase - cloudBaseDatabase
            #levels = np.arange(NlevelsMax)


            # step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
            #clouds = xr.Dataset(
            #    data_vars = {'cloudBase':(('time','levels'), cloudBaseDatabase),
            #                 'cloudTop':(('time','levels'), cloudTopDatabase),
            #                 'cloudThick': (('time', 'levels'), cloudThicknessDatabase)},
            #    coords    = {'levels':levels,
            #                 'time':time})

            #cloudsTimeWindow = clouds.sel(time=slice(timeStart, timeEnd))
            #PBLclouds = cloudsTimeWindow.where(cloudsTimeWindow.cloudTop < heightPBL, drop=True)
