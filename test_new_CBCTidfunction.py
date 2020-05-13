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

def f_testCBCTId(cloudMask, time, height):

    dimTime   = len(time)
    dimHeight = len(height)
    print(dimTime, dimHeight)
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
        #if len(np.where(column == 1.)[0][:]) == len(np.where(column == -1.)[0][:]):
        #    print('equal number of CB and CT found = ', len(np.where(column == -1.)[0][:]))
        #else:
        #    print('*********** something wrong *********')
    NCB = max(numberCB)
    NCT = max(numberCT)
    #print('max number of layers ', NCB)

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
    return(CBarray, CTarray, NCB, NlayersArray, gradBinary)


# setting parameters for calculating averaging and domain size of the model:
NprofilesOut                    = 24  # hourly means
timeIncrement                   = 60  # hourly means
patch                           = 'patch003'
flagPlot                        = 0

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/debugging/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)
date_arr                        = []
cloudMask_obs                   = np.zeros((Nfiles,9600, 150))
cloudMask_mod                   = np.zeros((Nfiles,9600, 150))

# days selected to test the algorithm
DaysArr                         = ['20130427']#,'20130428','20130518']
for iDayTest in range(len(DaysArr)):

    for indFile in range(Nfiles):
        date = fileListMod[indFile][87:95]

        if date == DaysArr[iDayTest]:
            print(indFile)

            filenameMod = fileListMod[indFile]
            filenameObs = fileListObs[indFile]
            #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
            #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'

            yy = int(date[0:4])
            mm = int(date[4:6])
            dd = int(date[6:8])
            timeStart = datetime.datetime(yy, mm, dd, 6)
            timeEnd = datetime.datetime(yy, mm, dd, 23)
            heightPBL = 2500.
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
            print(np.shape(cloudMask_obs[indFile,:,:]))

            # step 1: find all cloud base and cloud tops.
            cloudBAseTopAll = f_testCBCTId(cloudMask_obs[indFile,:,:], time, height)

            # step 2: build cloud database with all clouds saved as xarray dataset
            cloudBaseDatabase = cloudBAseTopAll[0]
            cloudTopDatabase  = cloudBAseTopAll[1]
            NlevelsMax        = cloudBAseTopAll[2]
            NlevelsArr        = cloudBAseTopAll[3]
            gradBinary        = cloudBAseTopAll[4]
            cloudThicknessDatabase    = cloudTopDatabase - cloudBaseDatabase

            # step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
            databaseCT = pd.DataFrame(cloudTopDatabase, index=time, columns=np.arange(5))
            mask_t = (databaseCT.index > timeStart) * (databaseCT.index < timeEnd)
            databaseCT_PBL = databaseCT.loc[mask_t,:]
            databaseCT_PBL.where(databaseCT_PBL.values < heightPBL, inplace = True)
            timePlotHalfHours = pd.date_range(date, periods=9600, freq='9s')

            fontSizeTitle = 12
            fontSizeX = 10
            fontSizeY = 10
            fontSizeCbar = 10
            labelsizeaxes = 10
            cbarAspect = 10
            Ymax = 11000
            Ymin = 107.
            fig, ax = plt.subplots(figsize=(10, 8))
            #ax = plt.subplot(1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
            cax = ax.pcolormesh(time, height, cloudMask_obs[0, :, :].transpose(), vmin=0, vmax=3,
                                  cmap=plt.cm.get_cmap("GnBu", 4))
            cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical', aspect=cbarAspect)
            cbar.ticks = ([0, 1, 2, 3])
            cbar.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
            cbar.set_label(label="cloud type", size=fontSizeCbar)
            cbar.ax.tick_params(labelsize=fontSizeCbar)
            plt.plot(databaseCT_PBL.index, databaseCT_PBL.values[:,0], color='red', label='CT1_PBL')
            plt.plot(databaseCT.index, databaseCT.values[:,0], color='black', label='CT1')
            plt.legend()

            plt.tight_layout()
            plt.savefig(pathFig+date+'_cloudMask_withCTFiltered.png')








PlotFlag_cloudMask = 0
# plot of cloud fraction time height for all days of simulation/obs
if PlotFlag_cloudMask == 1:
    for indDay in range(len(date_arr)):
        date = date_arr[indDay]
        yy = int(date[0:4])
        mm = int(date[4:6])
        dd = int(date[6:8])
        print('plotting date :' + date)
        timePlotHalfHours = pd.date_range(date, periods=9600, freq='9s')

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
        cax = ax.pcolormesh(timePlotHalfHours, height, cloudMask_mod[indDay,:,:].transpose(), vmin=0, vmax=3,
                            cmap=plt.cm.get_cmap("GnBu", 4))
        #ax.plot(timePlotHalfHours, CB_mod[indDay,:], color='black' )
        #ax.plot(timePlotHalfHours, CT_mod[indDay,:], color='red' )

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
        cax1 = ax1.pcolormesh(timePlotHalfHours, height, cloudMask_obs[indDay,:,:].transpose(), vmin=0, vmax=3,
                            cmap=plt.cm.get_cmap("GnBu", 4))
        #ax1.plot(timePlotHalfHours, CB_obs[indDay,:], color='black' )
        #ax1.plot(timePlotHalfHours, CT_obs[indDay,:], color='red' )
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


        fig, ax = plt.subplots(figsize=(10, 8))
        #ax = plt.subplot(1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        cax = ax.pcolormesh(time, height, gradBinary.transpose(), vmin=-1., vmax=1,
                            cmap=plt.cm.get_cmap("GnBu", 3))
        # ax.plot(timePlotHalfHours, CB_mod[indDay,:], color='black' )
        # ax.plot(timePlotHalfHours, CT_mod[indDay,:], color='red' )

        ax.set_ylim(Ymin, Ymax)  # limits of the y-axes
        ax.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        ax.set_title("binary model ", fontsize=fontSizeTitle)
        ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
        ax.set_ylabel("height [m]", fontsize=fontSizeY)
        cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='vertical', aspect=cbarAspect)
        cbar.ticks = ([-1, 0, 1])
        cbar.ax.set_yticklabels(['-1', '0', '1'])
        cbar.set_label(label='binary values', size=fontSizeCbar)
        cbar.ax.tick_params(labelsize=fontSizeCbar)
        plt.tight_layout()
        plt.savefig(pathFig + date + '_binaryMask.png')
