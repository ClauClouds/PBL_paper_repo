"""
date : wednesday 8 april 2020
author: Claudia Acquistapace
goal: code extracted from PBLpaper_statisticAnalysis_obsmod.py to produce the figure 1 of the paper, hourly profiles of
variance of vertical velocity and skewness for model and observations from radiosondes.
modified on 18 May 2020 to include CB, CT from PBL clouds based on new derivations
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
import xarray as xr
import matplotlib as mpl
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
from myFunctions import f_closest
from myFunctions import f_resampleArrays2StandardData
from myFunctions import f_resample2StandardData
from myFunctions import f_calculateMinCloudBaseTop

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
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListMod)
dataset_mean_variance_obs       = []
dataset_mean_variance_mod       = []
dataset_mean_skewness_obs       = []
dataset_mean_skewness_mod       = []

# dynamical properties
varianceW_obs                   = np.zeros((Nfiles,9600, 150))
varianceW_mod                   = np.zeros((Nfiles,9600, 150))
#stdW_obs                        = []
#stdW_mod                        = []
skewnessW_obs                   = np.zeros((Nfiles,9600, 150))
skewnessW_mod                   = np.zeros((Nfiles,9600, 150))
date_arr                        = []
datetime_out                    = []
PBLcloud_obs                    = []
PBLcloud_mod                    = []
clouds_obs                      = []
clouds_mod                      = []

for indFile in range(Nfiles):
#for indFile in range(1):
    print(indFile)

    date = fileListObs[indFile][81:89]
    print('processing date ' + date)
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)
    timeReferenceFormat             = pd.date_range(date, periods=9600, freq='9s')

    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(pathMod + 'icon_lem_derivedproperties' + date + '.nc', mode='r')
    time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]

    # opening the file containing all the data
    infile = open(pathObs + 'dictionaries_ModObs_' + date + '.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')


    # reading variance and skewness of vertical velocity. the function f_resample2StandardData resamples data with
    # dimtime < 9600 to the standard grid size of 9600, 150 filling matrices with nans
    varianceW_obs[indFile,:,:] = f_resample2StandardData(np.asarray(new_dict[3]['varianceW']), time[:], \
                                                         new_dict[3]['height'], date)
    varianceW_mod[indFile,:,:] = f_resample2StandardData(np.asarray(new_dict[9]['varianceW']), time[:], \
                                 new_dict[3]['height'], date)

    skewnessW_obs[indFile,:,:] = f_resample2StandardData(np.asarray(new_dict[3]['skewnessW']), time[:], \
                                                         new_dict[3]['height'], date)
    skewnessW_mod[indFile,:,:] = f_resample2StandardData(np.asarray(ncdata.groups['Temp_data'].variables['skewnessW'][:]),
                                                         time[:], new_dict[3]['height'], date)


    # reading CB and CT files for model and obs
    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_obs_nc = nc4.Dataset(pathObs + 'Clouds_Obs_' + date + '.nc', mode='r')  # Or from siphon.ncss
    cloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_obs_nc))
    dataset_obs = cloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    clouds_obs.append(dataset_obs)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_mod_nc = nc4.Dataset(pathMod + 'Clouds_iconlem_' + date + '.nc', mode='r')  # Or from siphon.ncss
    cloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_mod_nc))
    dataset_mod = cloudDay_mod.reindex({"time": timeReferenceFormat}, copy=True)
    clouds_mod.append(dataset_mod)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_obs_nc = nc4.Dataset(pathObs + 'PBLClouds_Obs_' + date + '.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_obs_nc))
    PBLdataset_obs = PBLcloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    PBLcloud_obs.append(PBLdataset_obs)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_mod_nc = nc4.Dataset(pathMod + 'PBLClouds_iconlem_' + date + '.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_mod_nc))
    PBLdataset_mod = PBLcloudDay_mod.reindex({"time": timeReferenceFormat}, copy=True)
    PBLcloud_mod.append(PBLdataset_mod)

# calculating mean variance and skewness for all observations and modeling
meanVariance_obs = np.nanmean(varianceW_obs, axis=0)
meanVariance_mod = np.nanmean(varianceW_mod, axis=0)
meanSkewness_obs = np.nanmean(skewnessW_obs, axis=0)
meanSkewness_mod = np.nanmean(skewnessW_mod, axis=0)

print('calculating CB and CT for observations')
CBarr_obs, CTarr_obs, TKarr_obs, CB_PBL_obs =  f_calculateMinCloudBaseTop(clouds_obs, PBLcloud_obs, date_arr)
print('calculating CB and CT for model')
CBarr_mod, CTarr_mod, TKarr_mod, CB_PBL_mod  =  f_calculateMinCloudBaseTop(clouds_mod, PBLcloud_mod, date_arr)


#averaging cloud bases of all days together
meanCBobs = np.nanmin(CB_PBL_obs, axis=1)
stdCBobs  = np.nanstd(CB_PBL_obs, axis=1)
meanCBmod = np.nanmin(CB_PBL_mod, axis=1)
stdCBmod  = np.nanstd(CB_PBL_mod, axis=1)



# calculating running means of cb arrays
timePlot = pd.date_range(date, periods=9600, freq='9s')
cb_mod_DF = pd.DataFrame(meanCBmod, index=timePlot)
cb_mod_std_DF = pd.DataFrame(stdCBmod, index=timePlot)
cb_mod_roll = cb_mod_DF.rolling(window=100).mean()
cb_mod_roll_std = cb_mod_std_DF.rolling(window=100).mean()

cb_obs_DF = pd.DataFrame(meanCBobs, index=timePlot)
cb_obs_std_DF = pd.DataFrame(stdCBobs, index=timePlot)
cb_obs_roll = cb_obs_DF.rolling(window=100).mean()
cb_obs_roll_std = cb_obs_std_DF.rolling(window=100).mean()
#%%
timePlotHours = pd.date_range(date, periods=24, freq='1h')
#totMeanCB_mod    = np.nanmean(MeanCB_mod, axis=0)
#totMeanCB_obs    = np.nanmean(MeanCB_obs, axis=0)
#totMaxCB_mod    = np.nanmean(maxCB_mod, axis=0)
#totMaxCB_obs    = np.nanmean(maxCB_obs, axis=0)
#totMinCB_mod    = np.nanmean(minCB_mod, axis=0)
#totMinCB_obs    = np.nanmean(minCB_obs, axis=0)
timeStart = datetime.datetime(yy,mm,dd,6)
timeEnd   = datetime.datetime(yy,mm,dd,23,0,0)
Ncols = 1
Nrows = 4
Nplots = 4
fontSizeTitle = 12
fontSizeX = 10
fontSizeY = 10
fontSizeCbar = 10
labelsizeaxes = 10
cbarAspect = 10
cbarSeL_sigma = 'viridis'
cbarSel_skn = 'bwr'
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(8, 10))
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax = 2500.
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
cax = ax.pcolormesh(timePlot, height, meanVariance_obs.transpose(), vmin=0., vmax=1., cmap=cbarSeL_sigma)
#ax.plot(timePlot,meanCBobs, color='lightgrey', linestyle=':')
#ax.plot(timePlot,meanCTobs, color='lightgrey')
ax.plot(cb_obs_roll.index, cb_obs_roll.values, color='lightgrey')
#ax.plot(cb_obs_roll.index, cb_obs_roll.values+cb_obs_roll_std.values, linestyle=':', color='lightgrey')
#ax.plot(cb_obs_roll.index, cb_obs_roll.values-cb_obs_roll_std.values, linestyle=':', color='lightgrey')
ax.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax.set_title('a) observations', fontsize=fontSizeTitle, loc='left')
ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
cbar.set_label(label="${\sigma^{2}} [m^{2}s^{-2}$]",size=fontSizeCbar)
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
cax1 = ax1.pcolormesh(timePlot, height, meanVariance_mod.transpose(), vmin=0., vmax=1., cmap=cbarSeL_sigma)
ax1.plot(cb_mod_roll.index, cb_mod_roll, color='lightgrey')

#ax1.plot(cb_mod_roll_std.index, cb_mod_roll.values-cb_mod_roll_std.values, linestyle=':', color='lightgrey')
#ax1.plot(cb_mod_roll_std.index, cb_mod_roll.values+cb_mod_roll_std.values, linestyle=':', color='lightgrey')
##ax1.plot(timePlot,meanCBmod, color='lightgrey', linestyle=':')
#ax1.plot(timePlot,meanCTmod, color='lightgrey')
ax1.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax1.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax1.set_title('b) model', fontsize=fontSizeTitle, loc='left')
ax1.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax1.set_ylabel("height [m]", fontsize=fontSizeY)
cbar1 = fig.colorbar(cax1, orientation='vertical', aspect=cbarAspect)
cbar1.set_label(label="${\sigma^{2}} [m^{2}s^{-2}$]",size=fontSizeCbar)
cbar1.ax.tick_params(labelsize=labelsizeaxes)
#cbar1.aspect=cbarAspect

ax2 = plt.subplot(Nrows, Ncols, 3)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis_date()
cax2 = ax2.pcolormesh(timePlot, height, meanSkewness_obs.transpose(), vmin=-3., vmax=3., cmap=cbarSel_skn)
ax2.plot(cb_obs_roll.index, cb_obs_roll.values, color='black')
#ax2.plot(cb_obs_roll.index, cb_obs_roll.values+cb_obs_roll_std.values, linestyle=':', color='black')
#ax2.plot(cb_obs_roll.index, cb_obs_roll.values-cb_obs_roll_std.values, linestyle=':', color='black')
#ax2.plot(timePlot,meanCBobs, color='lightgrey', linestyle=':')
#ax2.plot(timePlot,meanCTobs, color='lightgrey')
ax2.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax2.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax2.set_title('c) observations', fontsize=fontSizeTitle, loc='left')#
ax2.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax2.set_ylabel("height [m]", fontsize=fontSizeY)
cbar2 = fig.colorbar(cax2, orientation='vertical', aspect=cbarAspect)
cbar2.set_label(label="vert. vel. skewness",size=fontSizeCbar)
cbar2.ax.tick_params(labelsize=labelsizeaxes)
#cbar2.aspect=cbarAspect

ax3 = plt.subplot(Nrows, Ncols, 4)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis_date()
cax3 = ax3.pcolormesh(timePlot, height, meanSkewness_mod.transpose(), vmin=-3, vmax=3, cmap=cbarSel_skn)
ax3.plot(cb_mod_roll.index, cb_mod_roll, color='black')
#ax3.plot(cb_mod_roll_std.index, cb_mod_roll.values-cb_mod_roll_std.values, linestyle=':', color='black')
#ax3.plot(cb_mod_roll_std.index, cb_mod_roll.values+cb_mod_roll_std.values, linestyle=':', color='black')
##ax3.plot(timePlot,meanCBobs, color='lightgrey', linestyle=':')
#ax3.plot(timePlot,meanCTobs, color='lightgrey')
ax3.set_ylim(ymin,ymax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax3.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax3.set_title('d) model', fontsize=fontSizeTitle, loc='left')#
ax3.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax3.set_ylabel("height [m]", fontsize=fontSizeY)
cbar3 = fig.colorbar(cax3, orientation='vertical', aspect=cbarAspect)
cbar3.set_label(label="vert. vel. skewness",size=fontSizeCbar)
cbar3.ax.tick_params(labelsize=labelsizeaxes)
#cbar3.aspect=cbarAspect
plt.tight_layout()
plt.savefig(pathFig + 'figure1_var_W_maps.png', format='png')



if flagPlot == 1:
    Ncols = 5
    Nrows = 2
    Nplots = 11
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = 107.
    xmax          = 1.5
    fontSizeTitle = 16
    fontSizeX     = 10
    fontSizeY     = 10

    #timeTitles = [']
    indHourPlotStart = 8 # corresponding to starting at 8 UTC 
    
    for indPlot in range(1, Nplots):
        ax        = plt.subplot(2,5,indPlot)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left() 
        #ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        plt.plot(meanVariance_obs[:,indHourPlotStart], height, label='obs',  color='black')
        plt.plot(meanVariance_mod[:,indHourPlotStart], height, label='icon-lem',  color='red')
        plt.hlines(totMeanCB_obs[indHourPlotStart], 0.,xmax, color='black', linestyle=':', label='CB-obs')
        plt.hlines(totMeanCB_mod[indHourPlotStart], 0.,xmax, color='red', linestyle=':', label='CB-mod')
        y1        = meanVariance_obs[:,indHourPlotStart]-stdVariance_obs[:,indHourPlotStart]
        y2        = meanVariance_obs[:,indHourPlotStart]+stdVariance_obs[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = meanVariance_mod[:,indHourPlotStart]-stdVariance_mod[:,indHourPlotStart]
        y2        = meanVariance_mod[:,indHourPlotStart]+stdVariance_mod[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        if indPlot == 1:
            plt.legend(loc="upper right", fontsize=12, frameon=False)
        plt.ylim(ymin,ymax)
        plt.xlim(0.,xmax)
        plt.title(str(indHourPlotStart-1)+' UTC', fontsize=fontSizeTitle)
        plt.xlabel('${\sigma^{2}} [m^{2}s^{-2}$]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart+1
    
    plt.savefig(pathFig+'figure1_var_W_profiles.png', format='png')


    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = 107.
    xmax          = 1.5
    xmin = -1.5
    indHourPlotStart = 8 # corresponding to starting at 8 UTC

    for indPlot in range(1, Nplots):
        ax = plt.subplot(2, 5, indPlot)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        # ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=10)  # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)  # sets dimension of ticks in the plots
        plt.plot(meanSkewness_obs[:, indHourPlotStart], height, label='obs', color='black')
        plt.plot(meanSkewness_mod[:, indHourPlotStart], height, label='icon-lem', color='red')
        plt.hlines(totMeanCB_obs[indHourPlotStart], xmin, xmax, color='black', linestyle=':', label='CB-obs')
        plt.hlines(totMeanCB_mod[indHourPlotStart], xmin, xmax, color='red', linestyle=':', label='CB-mod')
        plt.axvline(x=0., color='grey', linestyle='-.')
        y1 = meanSkewness_obs[:, indHourPlotStart] - stdVariance_obs[:, indHourPlotStart]
        y2 = meanSkewness_obs[:, indHourPlotStart] + stdVariance_obs[:, indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2 > y1, facecolor='black', alpha=0.2)
        y1 = meanSkewness_mod[:, indHourPlotStart] - stdVariance_mod[:, indHourPlotStart]
        y2 = meanSkewness_mod[:, indHourPlotStart] + stdVariance_mod[:, indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2 > y1, facecolor='red', alpha=0.2)
        if indPlot == 1:
            plt.legend(loc="upper right", fontsize=12, frameon=False)
        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)
        plt.title(str(indHourPlotStart - 1) + ' UTC', fontsize=fontSizeTitle)
        plt.xlabel('w skewness ', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart + 1

    plt.savefig(pathFig + 'figure2_skn_W_profiles.png', format='png')

