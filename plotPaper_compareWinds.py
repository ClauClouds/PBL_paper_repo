

"""
Created on Wed 3 Jun 2020 10:57
@author: cacquist
goal : plot mean zonal and meridional wind maps for model and obs

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
from myFunctions import f_calculateMinCloudBaseTop

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
patch = 'patch003'

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListMod)
date_arr = []
merWind_mod_All = np.zeros((9600,150,Nfiles))
merWind_obs_All = np.zeros((9600,150,Nfiles))
zonWind_mod_All = np.zeros((9600,150,Nfiles))
zonWind_obs_All = np.zeros((9600,150,Nfiles))
vertWind_mod_All = np.zeros((9600,151,Nfiles))
vertWind_obs_All = np.zeros((9600,151,Nfiles))

for indFile in range(Nfiles):
    print(indFile)
    date = fileListObs[indFile][81:89]
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)

    print('processing date ' + date)
    ncdata = Dataset(pathMod + 'icon_lem_derivedproperties' + date + '.nc', mode='r')
    time = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]
    height_w = ncdata.groups['Temp_data'].variables['height'][:]
    # opening the file containing all the data
    infile = open(pathObs + 'dictionaries_ModObs_' + date + '.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')


    # reading meridional and zonal winds for model and obs
    merWind_mod = xr.Dataset({'Wind': (['time', 'height'], ncdata.groups['Temp_data'].variables['merWind'][:])},
                        coords = {'time': time,
                                  'height':height})
    zonWind_mod = xr.Dataset({'Wind': (['time', 'height'], ncdata.groups['Temp_data'].variables['zonalWind'][:])},
                        coords = {'time': time,
                                  'height':height})
    vertWind_mod = xr.Dataset({'Wind': (['time', 'height_w'], ncdata.groups['Temp_data'].variables['vertWind'][:])},
                        coords = {'time': time,
                                  'height':height_w})

    merWind_obs = xr.Dataset({'Wind': (['time', 'height'], new_dict[3]['meridionalWind'])},
                        coords = {'time': time,
                                  'height':height})
    zonWind_obs = xr.Dataset({'Wind': (['time', 'height'], new_dict[3]['zonalWind'])},
                        coords = {'time': time,
                                  'height':height})
    vertWind_obs = xr.Dataset({'Wind': (['time', 'height'], new_dict[3]['verticalWind'])},
                              coords={'time'  : time,
                                      'height': height})

    # reshaping on standard time format
    timeReferenceFormat = pd.date_range(date, periods=9600, freq='9s')

    merWind_mod_reshaped = merWind_mod.reindex({"time": timeReferenceFormat}, copy=True)
    merWind_mod_All[:,:,indFile] = merWind_mod_reshaped.Wind.values

    merWind_obs_reshaped = merWind_obs.reindex({"time": timeReferenceFormat}, copy=True)
    merWind_obs_All[:, :, indFile] = merWind_obs_reshaped.Wind.values

    zonWind_mod_reshaped = zonWind_mod.reindex({"time": timeReferenceFormat}, copy=True)
    zonWind_mod_All[:,:,indFile] = zonWind_mod_reshaped.Wind.values

    zonWind_obs_reshaped = zonWind_obs.reindex({"time": timeReferenceFormat}, copy=True)
    zonWind_obs_All[:, :, indFile] = zonWind_obs_reshaped.Wind.values

    vertWind_obs_reshaped = vertWind_obs.reindex({"time": timeReferenceFormat}, copy=True)
    vertWind_mod_reshaped = vertWind_mod.reindex({"time": timeReferenceFormat}, copy=True)



    # reshaping w wind obs on the model height format
    vertWind_obs_reshaped_h = vertWind_obs_reshaped.interp_like(vertWind_mod_reshaped)
    print('check data')
    print(np.nanmin(vertWind_obs_reshaped_h.Wind.values))
    print(np.nanmax(vertWind_obs_reshaped_h.Wind.values))

    vertWind_obs_All[:,:,indFile] = vertWind_obs_reshaped_h.Wind.values
    vertWind_mod_All[:,:,indFile] = vertWind_mod_reshaped.Wind.values



mean_merWind_mod = np.nanmean(merWind_mod_All, axis=2)
mean_merWind_obs = np.nanmean(merWind_obs_All, axis=2)
mean_zonWind_mod = np.nanmean(zonWind_mod_All, axis=2)
mean_zonWind_obs = np.nanmean(zonWind_obs_All, axis=2)
mean_vertWind_obs =  np.nanmean(vertWind_obs_All, axis=2)
mean_vertWind_mod =  np.nanmean(vertWind_mod_All, axis=2)
print(np.nanmin(mean_vertWind_obs))
print(np.nanmin(mean_vertWind_mod))
print(np.nanmax(mean_vertWind_obs))
print(np.nanmax(mean_vertWind_mod))

# plotting biases obs - mod
bias_merWind = mean_merWind_obs - mean_merWind_mod
bias_zonWind = mean_zonWind_obs - mean_zonWind_mod
bias_vertWind = mean_vertWind_obs - mean_vertWind_mod

print(bias_vertWind)
# calculating running means of cb arrays
timePlot = pd.date_range(date, periods=9600, freq='9s')
timePlotHours = pd.date_range(date, periods=24, freq='1h')
timeStart = datetime.datetime(yy,mm,dd,6)
timeEnd   = datetime.datetime(yy,mm,dd,23,59,59)
Ncols = 1
Nrows = 3
Nplots = 3
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
ymax = 2000.
ymin = 107.
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
cax = ax.pcolormesh(timePlot, height, bias_merWind.transpose(), vmin=-2., vmax=2., cmap=cbarSel_v)
ax.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax.set_title('a) bias meridional wind', fontsize=fontSizeTitle, loc='left')
ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
cbar.set_label(label="$v_{obs}-v_{mod} [ms^{-1}$]",size=fontSizeCbar)
cbar.ax.tick_params(labelsize=labelsizeaxes)
#cbar.aspect=cbarAspect

ax1 = plt.subplot(Nrows, Ncols, 2)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis_date()
cax1 = ax1.pcolormesh(timePlot, height, bias_zonWind.transpose(), vmin=-2., vmax=2., cmap=cbarSeL_u)
ax1.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax1.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax1.set_title('b) bias zonal wind', fontsize=fontSizeTitle, loc='left')
ax1.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax1.set_ylabel("height [m]", fontsize=fontSizeY)
cbar1 = fig.colorbar(cax1, orientation='vertical', aspect=cbarAspect)
cbar1.set_label(label="$u_{obs}-u_{mod} [ms^{-1}$]",size=fontSizeCbar)
cbar1.ax.tick_params(labelsize=labelsizeaxes)


ax2 = plt.subplot(Nrows, Ncols, 3)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis_date()
cax2 = ax2.pcolormesh(timePlot, height_w, bias_vertWind.transpose(), vmin=-2., vmax=2., cmap=cbarSeL_u)
ax2.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax2.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax2.set_title('b) bias vertical wind', fontsize=fontSizeTitle, loc='left')
ax2.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax2.set_ylabel("height [m]", fontsize=fontSizeY)
cbar2 = fig.colorbar(cax2, orientation='vertical', aspect=cbarAspect)
cbar2.set_label(label="$w_{obs}-w_{mod} [ms^{-1}$]",size=fontSizeCbar)
cbar2.ax.tick_params(labelsize=labelsizeaxes)
plt.tight_layout()
plt.savefig(pathFig + 'figure2_biasWinds.png', format='png')