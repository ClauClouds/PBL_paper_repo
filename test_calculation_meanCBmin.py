
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
import xarray as xr
import matplotlib as mpl
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
from myFunctions import f_closest
from myFunctions import f_resampleArrays2StandardData
from myFunctions import f_resample2StandardData


try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

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
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListMod)
dataset_mean_variance_obs       = []
dataset_mean_variance_mod       = []
dataset_mean_skewness_obs       = []
dataset_mean_skewness_mod       = []

date_arr                        = []
datetime_out                    = []
PBLcloud_obs                   = []
clouds_obs = []

from myFunctions import f_calculateMinCloudBaseTop
#def f_calculateMinCloudBaseTop(clouds, PBLclouds, dimTime, date_arr):
#    """author: claudia Acquistapace
#     date: 18/05/2020
#     goal: function to calculate the minimum cloud base for the PBL and the corresponding cloud top
#     input: clouds - list of xarray datasets of cloud properties
#            PBLclouds - list of xarray datasets of PBL cloud properties
#            date_arr - array of days to be processed
#            dimTime - time array dimension
#    output: minimum cloud base and corresponding cloud tops
#    """
#    # definition of output matrices
#    Nfiles = len(date_arr)
#    CBarr_obs = np.zeros((dimTime, Nfiles))
#    CTarr_obs = np.zeros((dimTime, Nfiles))


    # for each day, reading and saving minimum cloud base and corresponding cloud top
#   for indFile in range(Nfiles):
#        # calculating CB and CT for PBL
#        iPBL_obs = 0
#        print(date_arr[indFile])
#        for itime in range(dimTime):
#            Obs_dataset = PBLcloud_obs[indFile]
#            cloud_dataset = clouds_obs[indFile]
#            if iPBL_obs < len(Obs_dataset.time.values):
#                if cloud_dataset.time.values[itime] == Obs_dataset.time.values[iPBL_obs]:
#                    CBarr_obs1 = PBLcloud_obs[indFile].cloudBase.values[iPBL_obs, :]
##                    if CBarr_obs1.size - np.count_nonzero(np.isnan(CBarr_obs1)) != 0:
#                        minCB_obs = np.nanmin(Obs_dataset.cloudBase.values[iPBL_obs, :])
#                        CBarr_obs[itime, indFile] = minCB_obs
#                        indexLevelMin_obs = np.nanargmin(PBLcloud_obs[indFile].cloudBase.values[iPBL_obs, :])
#                        CTarr_obs[itime, indFile] = PBLcloud_obs[indFile].cloudTop[iPBL_obs, indexLevelMin_obs]
#                   iPBL_obs = iPBL_obs + 1

#    return (CBarr_obs, CTarr_obs)


for indFile in range(Nfiles):
    print(indFile)

    date = fileListObs[indFile][81:89]
    print('processing date ' + date)
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)
    timeReferenceFormat             = pd.date_range(date, periods=9600, freq='9s')

    print(date)
    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_obs_nc = nc4.Dataset(pathObs + 'Clouds_Obs_' + date + '.nc', mode='r')  # Or from siphon.ncss
    cloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_obs_nc))
    cloud_dataset = cloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    clouds_obs.append(cloud_dataset)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_obs_nc = nc4.Dataset(pathObs + 'PBLClouds_Obs_' + date + '.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_obs_nc))
    dataset = PBLcloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    PBLcloud_obs.append(dataset)


from myFunctions import f_calculateMinCloudBaseTop
CBarr_obs, CTarr_obs, TKarr_obs, CB_PBL_arr  =  f_calculateMinCloudBaseTop(clouds_obs, PBLcloud_obs, date_arr)
##
print(np.shape(CBarr_obs))
print(np.shape(CB_PBL_arr))

meanCB = np.nanmin(CBarr_obs, axis=1)
stdCB = np.nanstd(CBarr_obs, axis=1)
meanCB_PBL = np.nanmin(CB_PBL_arr, axis=1)
stdCB_PBL = np.nanstd(CB_PBL_arr, axis=1)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
axes = plt.subplot(1, 1, 1)
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)
axes.get_xaxis().tick_bottom()
axes.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=10)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=10)  # sets dimension of ticks in the plots
#plt.plot(timeReferenceFormat, meanCB+stdCB, linestyle=':', color='black')
plt.plot(timeReferenceFormat, meanCB, color='black', label="all")
plt.plot(timeReferenceFormat, meanCB_PBL, color='red', label='PBL')
#plt.plot(timeReferenceFormat, meanCB-stdCB, linestyle=':', color='black')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(pathFig + 'TestMinCB.png', format='png')