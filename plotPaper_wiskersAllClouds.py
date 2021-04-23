
"""
date : 04/06/2020
author: Claudia Acquistapace
goal: calculate distributions of cloud base/cloud tops/thickness values every 15 minutes to overplot on time series of other variables
"""

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
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListObs)
date_arr                        = []
cloudMask_obs                   = np.zeros((Nfiles,9600, 150))
cloudMask_mod                   = np.zeros((Nfiles,9600, 150))
cloudFractionTot_obs            = np.zeros((Nfiles,48, 150))
cloudFractionTot_mod            = np.zeros((Nfiles,48, 150))
clouds_obs                      = []
clouds_mod                      = []
PBLclouds_obs                   = []
PBLclouds_mod                   = []


cloudMatrix_base_obs = np.zeros((9600,8,Nfiles))
cloudMatrix_base_mod = np.zeros((9600,8,Nfiles))
Nlevels_example = np.arange(8)

for indFile in range(Nfiles):

    print(indFile)

    #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
    #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'
    date = fileListObs[indFile][81:89]
    filenameMod = pathMod+'icon_lem_derivedproperties'+date+'.nc'
    filenameObs = pathObs+'dictionaries_ModObs_'+date+'.p'
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])

    date_arr.append(date)
    # defining reference timearray for resizing cloud bases/tops/thicknesses
    print('processing date ' + date)
    # reading CB and CT files for model and obs

    timeReferenceFormat = pd.date_range(date, periods=9600, freq='9s')
    Dataset_example = xr.Dataset({'cloudBase': (['time', 'levels'], np.random.rand(len(timeReferenceFormat),len(Nlevels_example)))},
                                 coords={'time'  : timeReferenceFormat,
                                         'levels': Nlevels_example})
    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_obs_nc = nc4.Dataset(pathObs+'Clouds_Obs_'+date+'.nc', mode='r')  # Or from siphon.ncss
    cloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_obs_nc))
    #print(cloudDay_obs)
    cloudStandard_obs = cloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    cloudStandard_obs_reshaped = cloudStandard_obs.interp_like(Dataset_example)
    cloudMatrix_base_obs[:,:,indFile] = cloudStandard_obs_reshaped.cloudBase.values

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_mod_nc = nc4.Dataset(pathMod+'Clouds_iconlem_'+date+'.nc', mode='r')  # Or from siphon.ncss
    cloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_mod_nc))
    #clouds_mod.append(cloudDay_mod)
    cloudStandard_mod = cloudDay_mod.reindex({"time": timeReferenceFormat}, copy=True)
    cloudStandard_mod_reshaped = cloudStandard_mod.interp_like(Dataset_example)
    cloudMatrix_base_mod[:,:,indFile] = cloudStandard_mod_reshaped.cloudBase.values
    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    #PBLclouds_obs_nc = nc4.Dataset(pathObs+'PBLClouds_Obs_'+date+'.nc', mode='r')  # Or from siphon.ncss
    #PBLcloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_obs_nc))
    #PBLclouds_obs.append(PBLcloudDay_obs)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    #PBLclouds_mod_nc = nc4.Dataset(pathMod+'PBLClouds_iconlem_'+date+'.nc', mode='r')  # Or from siphon.ncss
    #PBLcloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_mod_nc))
    #PBLclouds_mod.append(PBLcloudDay_mod)

    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(filenameMod, mode='r')
    time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]

#%%
# define xr dataset for cloud base obs
cloudBaseAll_obs_xrDataset = xr.Dataset({'cloudBase': (['time', 'levels','Ndays'], cloudMatrix_base_obs)},
                              coords={'time'  : timeReferenceFormat,
                                      'levels': Nlevels_example,
                                      'Ndays': date_arr,
                                      'height': height})
cloudBaseAll_mod_xrDataset = xr.Dataset({'cloudBase': (['time', 'levels','Ndays'], cloudMatrix_base_mod)},
                              coords={'time'  : timeReferenceFormat,
                                      'levels': Nlevels_example,
                                      'Ndays': date_arr,
                                      'height': height})

# deriving cloud base PDF for every DeltaT interval.
dateSampler =  pd.date_range(timeReferenceFormat[0], timeReferenceFormat[-1], freq='15min')
Nbins = 30
matrix_2d_mod = np.zeros((Nbins,len(dateSampler)))
matrix_2d_obs = np.zeros((Nbins,len(dateSampler)))

for indMin in range(len(dateSampler)-1):
    selection_mod = cloudBaseAll_mod_xrDataset.sel(time=slice(dateSampler[indMin], dateSampler[indMin+1]))
    distrib_mod = selection_mod.cloudBase.values
    counts_mod, yBins_mod = np.histogram(distrib_mod.flatten(), range=[0., 3000.], bins=Nbins, density=True)

    selection_obs = cloudBaseAll_obs_xrDataset.sel(time=slice(dateSampler[indMin], dateSampler[indMin+1]))
    distrib_obs = selection_obs.cloudBase.values
    counts_obs, yBins_obs = np.histogram(distrib_obs.flatten(), range=[0., 3000.], bins=Nbins, density=True)
    yCenterBins_mod = np.diff(yBins_mod) + yBins_mod[:-1]
    yCenterBins_obs = np.diff(yBins_obs) + yBins_obs[:-1]

    matrix_2d_mod[:,indMin] = counts_mod
    matrix_2d_obs[:,indMin] = counts_obs
#%%

# plot minimum cloud base for every day to understand lower clouds at 16:00
minCloudBase_all_obs = np.zeros((len(timeReferenceFormat),len(date_arr)))
minCloudBase_all_mod = np.zeros((len(timeReferenceFormat),len(date_arr)))

for indDay in range(len(date_arr)):
    print('processing day ', indDay)
    for indTime in range(len(timeReferenceFormat)):
        minCloudBase_all_obs[indTime,indDay] = np.nanmin(cloudBaseAll_obs_xrDataset.cloudBase.values[indTime,:,indDay])
        minCloudBase_all_mod[indTime,indDay] = np.nanmin(cloudBaseAll_mod_xrDataset.cloudBase.values[indTime,:,indDay])

#%%
labelsizeaxes = 12
fig, ax = plt.subplots(figsize=(14, 5))
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax = 3000.
ymin = 107.

rcParams["legend.loc"] = 'upper right' 
timeStart = dateSampler[0]
timeEnd = dateSampler[-1]
ax = plt.subplot(1, 1, 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
for ind in range(len(date_arr)):
    ax.plot(timeReferenceFormat, minCloudBase_all_obs[:,ind], 's', label=date_arr[ind])
ax.set_ylim(ymin, ymax)
ax.grid(True, linestyle=':')

ax.set_xlabel('time [hh:mm]', fontsize=12)
ax.set_ylabel('Height [m]', fontsize=12)
ax.set_title('cloud base for 22 days- obs')
colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
colors = [colormap(i) for i in np.linspace(0, 1,len(ax.lines))]
for i,j in enumerate(ax.lines):
    j.set_color(colors[i])
plt.legend(frameon=True,facecolor='white', framealpha=1,ncol=2,handleheight=2.4, labelspacing=0.05)
plt.tight_layout()
plt.savefig(pathFig + 'allCB_obs.png', format='png')

fig, ax = plt.subplots(figsize=(14, 5))
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax = 3000.
ymin = 107.

rcParams["legend.loc"] = 'upper right' 
timeStart = dateSampler[0]
timeEnd = dateSampler[-1]
ax = plt.subplot(1, 1, 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
for ind in range(len(date_arr)):
    ax.plot(timeReferenceFormat, minCloudBase_all_mod[:,ind],'s', label=date_arr[ind])
ax.set_ylim(ymin, ymax)
ax.grid(True, linestyle=':')
ax.set_title('cloud base for 22 days- mod')
ax.set_xlabel('time [hh:mm]', fontsize=12)
ax.set_ylabel('Height [m]', fontsize=12)
colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
colors = [colormap(i) for i in np.linspace(0, 1,len(ax.lines))]
for i,j in enumerate(ax.lines):
    j.set_color(colors[i])
plt.legend(frameon=True,facecolor='white', framealpha=1, ncol=2,handleheight=2.4, labelspacing=0.05)
plt.tight_layout()
plt.savefig(pathFig + 'allCB_mod.png', format='png')
#%%
# plot 2d distribution of clodu bases for 15 min interval
Ncols = 1
Nrows = 2
Nplots = 2
fontSizeTitle = 12
fontSizeX = 10
fontSizeY = 10
fontSizeCbar = 10
labelsizeaxes = 10
cbarAspect = 10
cbarSeL_u = 'RdYlBu'#'bwr'
cbarSel_v = 'RdYlBu'#'bwr'
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(8, 10))
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax = 3000.
ymin = 107.
timeStart = dateSampler[0]
timeEnd = dateSampler[-1]
ax = plt.subplot(Nrows, Ncols, 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax = ax.pcolormesh(dateSampler, yCenterBins_obs, matrix_2d_obs, vmin=0., vmax=0.001, cmap='PuRd')
ax.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax.set_title('cloud base distribution for 15 min time interval (obs)', fontsize=fontSizeTitle, loc='left')
ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
cbar.set_label(label="Occurrences ",size=fontSizeCbar)
cbar.ax.tick_params(labelsize=labelsizeaxes)
#cbar.aspect=cbarAspect


ax = plt.subplot(Nrows, Ncols, 2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax = ax.pcolormesh(dateSampler, yCenterBins_mod, matrix_2d_mod, vmin=0., vmax=0.001, cmap='PuRd')
ax.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax.set_title('cloud base distribution for 15 min time interval (model)', fontsize=fontSizeTitle, loc='left')
ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
cbar.set_label(label="Occurrences ",size=fontSizeCbar)
cbar.ax.tick_params(labelsize=labelsizeaxes)

plt.tight_layout()
plt.savefig(pathFig + '2d_CB_distr.png', format='png')