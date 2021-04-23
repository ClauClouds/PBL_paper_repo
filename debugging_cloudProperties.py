
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
import xarray as xr
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
from myFunctions import f_closest
from myFunctions2 import hourDecimal_to_datetime
from cloudnetFunctions import f_calculateCloudMaskCloudnet
from myFunctions import f_selectingPBLcloudWindow
from myFunctions import f_calculateCloudBaseTopThickness
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle


# setting parameters for calculating averaging and domain size of the model:
NprofilesOut                    = 24  # hourly means
timeIncrement                   = 60  # hourly means
patch                           = 'patch003'



# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
path_cloudnet_cat   = '/data/hatpro/jue/cloudnet/juelich/processed/categorize/2013/'
pathDebugFig = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/debugging/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListObs)
date_arr                        = []
#filenameObs = 'dataset_PBLcloudPaper_ModObs_20130502.p'
#filenameMod = 'icon_lem_derivedproperties20130502.nc'

fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)

# ----------------------------------------------------------------------------------------
# loop on the ensemble of days (statistics)
# ----------------------------------------------------------------------------------------
for indFile in range(Nfiles):
    print(indFile)

    date = fileListObs[indFile][81:89]
    print(date)

    filenameMod = pathMod+'icon_lem_derivedproperties'+date+'.nc'
    filenameObs = pathObs+'dictionaries_ModObs_'+date+'.p'
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])

    date_arr.append(date)
    print('processing date ' + date)

    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(filenameMod, mode='r')
    time = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]

    # -----------------------------------------------------------------------------------
    # ---- cloudnet classification
    # -----------------------------------------------------------------------------------
    pathIn = path_cloudnet_cat
    filename = date + '_juelich_categorize.nc'
    if (os.path.isfile(pathIn + filename)):
        print('cloudnet class found: reading and resampling data')
        CLOUDNET_data = Dataset(pathIn + filename, mode='r')

        # ----- reading CLOUDNET data variables
        time_CLOUDNET = CLOUDNET_data.variables['time'][:].copy()
        datetime_CLOUDNET     = hourDecimal_to_datetime(int(yy), int(mm), int(dd), time_CLOUDNET)

        height_CLOUDNET = CLOUDNET_data.variables['height'][:].copy()
        print(np.shape(time_CLOUDNET))
        print(np.shape(height_CLOUDNET))

        cloudnet = CLOUDNET_data.variables['category_bits'][:].copy()
        print(np.shape(cloudnet))
        P_CLOUDNET = CLOUDNET_data.variables['pressure'][:].copy()  # [Pa]
        T_CLOUDNET = CLOUDNET_data.variables['temperature'][:].copy()  # [K]
        Q_CLOUDNET = CLOUDNET_data.variables['specific_humidity'][:].copy()  # [Kg/Kg]
        model_Height_CLOUDNET = CLOUDNET_data.variables['model_height'][:].copy()
        Z_CLOUDNET = CLOUDNET_data.variables['Z'][:].copy()
        CB_CLOUDNET = CLOUDNET_data.variables['Z'][:].copy()

        # calculating cloud mask for observations
        cloudMask         = f_calculateCloudMaskCloudnet(time_CLOUDNET, height_CLOUDNET, \
                                        cloudnet.astype(int))


        # plotting cloud mask
        fig, ax1 = plt.subplots(figsize=(14, 10))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis_date()
        cax = ax1.pcolormesh(datetime_CLOUDNET, height_CLOUDNET, cloudMask.transpose(), vmin=0., vmax=np.nanmax(cloudMask), cmap='PiYG')
        ax1.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
        #ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        # ax1.set_title("mean skewness W obs", fontsize=16)#
        ax1.set_xlabel("time [hh:mm]", fontsize=12)
        ax1.set_ylabel("height [m]", fontsize=12)
        plt.tight_layout()
        plt.savefig(pathDebugFig+date+'_cloudMask_Cloudnet.png', format='png')

        # calculating cloud base , cloud top and cloud thickness for all clouds and for pbl clouds
        humanInfo = f_selectingPBLcloudWindow(date)
        print(humanInfo['timeStart'])
        print(humanInfo['timeEnd'])
        print(datetime_CLOUDNET[0])
        print(datetime_CLOUDNET[-1])
        clouds, PBLclouds = f_calculateCloudBaseTopThickness(cloudMask, datetime_CLOUDNET, height_CLOUDNET, humanInfo)

        # deriving lowest cloud base and corresponding cloud top for PBL clouds
        CBarr = np.zeros(len(datetime_CLOUDNET))
        CBarr.fill(np.nan)
        CTarr = np.zeros(len(datetime_CLOUDNET))
        CTarr.fill(np.nan)
        iPBL = 0
        for itime in range(len(datetime_CLOUDNET)):
            if iPBL < len(PBLclouds.time.values):
                if clouds.time.values[itime] == PBLclouds.time.values[iPBL]:
                    #print(iPBL)
                    CBarray = PBLclouds.cloudBase.values[iPBL, :]
                    if CBarray.size - np.count_nonzero(np.isnan(CBarray)) != 0:
                        minCB = np.nanmin(PBLclouds.cloudBase.values[iPBL, :])
                        CBarr[itime] = minCB
                        indexLevelMin = np.nanargmin(PBLclouds.cloudBase.values[iPBL, :])
                        CTarr[itime] = PBLclouds.cloudTop[iPBL, indexLevelMin]
                    else:
                        CBarr[itime] = np.nan
                        CTarr[itime] = np.nan
                    iPBL = iPBL + 1



    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_obs_nc = nc4.Dataset(pathObs+'PBLClouds_Obs_'+date+'.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_obs_nc))

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_mod_nc = nc4.Dataset(pathMod+'PBLClouds_iconlem_'+date+'.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_mod_nc))

    CB_mod_rescaled = PBLcloudDay_mod.cloudBase.values
    CB_obs_rescaled = PBLcloudDay_obs.cloudBase.values
    timeStart = humanInfo['timeStart']
    timeEnd = humanInfo['timeEnd']
    timePlot = pd.date_range(datetime.datetime(yy,mm,dd,7,0,0), datetime.datetime(yy,mm,dd,18,0,0), freq='9s')

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis_date()
    ax.plot(datetime_CLOUDNET, CBarr, 'o', label='cloudnet not resampled', color='blue')
    for indLev in range(8):
        ax.plot(timePlot, CB_mod_rescaled[:,indLev], 'o', label='model rescaled', color='red')
        ax.plot(timePlot, CB_obs_rescaled[:,indLev], 'o', label='obs rescaled', color='black')
        if indLev == 0:
            ax.legend(frameon=False)
    ax.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
    ax.set_xlim(timeStart, timeEnd)  # limits of the x-axes
    ax.set_title("cloud base [m]", fontsize=16)  #
    ax.set_xlabel("time [hh:mm]", fontsize=12)
    ax.set_ylabel("height [m]", fontsize=12)
    ax.grid(True, which="both")
    ax.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
    # ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
    ax.grid(True, which="both")
    ax.set_title("cloud base [m]", fontsize=16)  #
    ax.set_xlabel("time [hh:mm]", fontsize=12)
    ax.set_ylabel("height [m]", fontsize=12)
    plt.tight_layout()
    plt.savefig(pathDebugFig + date + '_cloudBase_mod_obs_debugging.png', format='png')




