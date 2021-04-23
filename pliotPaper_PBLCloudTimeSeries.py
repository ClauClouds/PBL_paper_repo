"""
date : Friday 19 june 2020
author: Claudia Acquistapace
goal: calculate cloud fraction for PBL clouds and mean CB /CT time series with corresponding LWP time series for PBL clouds only
procedure:
method 1:
- identify max cloud top and min cloud base for every time of the PBL cloud time interval
- reshape cloud mask in the PBL cloud time interval
- set to zero all cloud fraction out of the height range established by cloud base/top
method 2:
- define a new cloud mask based on PBL cloud tops and bases.
- calculate cloud fraction providing the new cloud mask as imput


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
from myFunctions import f_calPBLcloudMask
from myFunctions import f_calculateCloudFractionPBLclouds
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
flag_plotCloudFraction = 0
flag_plotCloudMask = 0

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
pathDebugFig = pathFig+'/debugging/'
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
CF_mod_arr = []
CF_obs_arr = []
for indFile in range(Nfiles):

    print(indFile)
    #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
    #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'
    date = fileListObs[indFile][81:89]
    print(date)

    filenameMod = pathMod+'icon_lem_derivedproperties'+date+'.nc'
    filenameObs = pathObs+'dictionaries_ModObs_'+date+'.p'
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])

    date_arr.append(date)

    # reading CB and CT files for model and obs
    print('processing date ' + date)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_obs_nc = nc4.Dataset(pathObs+'PBLClouds_Obs_'+date+'.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_obs_nc))
    PBLclouds_obs.append(PBLcloudDay_obs)

    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    PBLclouds_mod_nc = nc4.Dataset(pathMod+'PBLClouds_iconlem_'+date+'.nc', mode='r')  # Or from siphon.ncss
    PBLcloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(PBLclouds_mod_nc))
    PBLclouds_mod.append(PBLcloudDay_mod)

    # reading height array for model and for observations
    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(filenameMod, mode='r')
    time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]
    PBLtime = PBLcloudDay_mod.time.values

    PBLcloudMask_mod = f_calPBLcloudMask(PBLcloudDay_mod,PBLtime,height)
    PBLcloudMask_obs = f_calPBLcloudMask(PBLcloudDay_obs,PBLtime,height)


    # defining a new cloud mask based on the PBl cloud bases and tops
    #CB_matrix = PBLcloudDay_mod.cloudBase.values
    #CT_matrix = PBLcloudDay_mod.cloudTop.values
    #PBLtime   = PBLcloudDay_mod.time.values

    #PBL_cloudMask = np.zeros((len(PBLtime),len(height)))
    #for indTime in range(len(PBLtime)):
    #    for indLev in range(8):
    #        if (~np.isnan(CB_matrix[indTime, indLev])):
    #            cb_height = CB_matrix[indTime, indLev]
    #            ct_height = CT_matrix[indTime, indLev]
    #            ind = (height >= cb_height) * (height <= ct_height)
    #            PBL_cloudMask[indTime,ind] = 1


    # plot of cloud mask for PBL clouds
    if flag_plotCloudMask == 1:
        fig, ax1 = plt.subplots(figsize=(14, 10))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis_date()
        cax = ax1.pcolormesh(PBLtime, height, PBLcloudMask_mod.transpose(), vmin=0., vmax=1., cmap='PiYG')
        ax1.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
        #ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        # ax1.set_title("mean skewness W obs", fontsize=16)#
        ax1.set_xlabel("time [hh:mm]", fontsize=12)
        ax1.set_ylabel("height [m]", fontsize=12)
        plt.tight_layout()
        plt.savefig(pathDebugFig+'cloudMaskPBL_mod_'+date+'.png', format='png')

        fig, ax1 = plt.subplots(figsize=(14, 10))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis_date()
        cax = ax1.pcolormesh(PBLtime, height, PBLcloudMask_obs.transpose(), vmin=0., vmax=1., cmap='PiYG')
        ax1.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
        #ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        # ax1.set_title("mean skewness W obs", fontsize=16)#
        ax1.set_xlabel("time [hh:mm]", fontsize=12)
        ax1.set_ylabel("height [m]", fontsize=12)
        plt.tight_layout()
        plt.savefig(pathDebugFig+'cloudMaskPBL_obs_'+date+'.png', format='png')

    # calculating cloud fraction every 15 minutes
    Nmin_string = '15'
    CF_mod_Dataset = f_calculateCloudFractionPBLclouds(PBLcloudMask_mod,PBLtime,height,Nmin_string)
    CF_obs_Dataset = f_calculateCloudFractionPBLclouds(PBLcloudMask_obs,PBLtime,height,Nmin_string)
    CF_mod_arr.append(CF_mod_Dataset)
    CF_obs_arr.append(CF_obs_Dataset)

    #cloudMask_DF = pd.DataFrame(PBLcloudMask_mod, index=PBLtime, columns=height)
    #datetime_CF = pd.date_range(start=PBLtime[0], end=PBLtime[-1], freq='15min')
    #datetime_CF = datetime_CF.to_pydatetime()

    #CF_PBL = np.zeros((len(datetime_CF), len(height)))
    #for indTime in range(len(datetime_CF)-1):
    #    mask_t = (cloudMask_DF.index > datetime_CF[indTime]) * (cloudMask_DF.index < datetime_CF[indTime+1])
    #    Selection_cloudMask = cloudMask_DF[mask_t]
    #    for indHeight in range(len(height)):
    #        CFArray = Selection_cloudMask.loc[:,Selection_cloudMask.columns[indHeight]]
    #        CF_PBL[indTime,indHeight] = len(CFArray[CFArray == 1])/len(CFArray)


    # plot of cloud fraction time height
    if flag_plotCloudFraction == 1:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(14, 8))
        ax1 = plt.subplot(2, 1, 1)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis_date()
        cax = ax1.pcolormesh(CF_obs_Dataset.time.values, CF_obs_Dataset.height.values, CF_obs_Dataset.CF.values.transpose(), vmin=0., vmax=0.5, cmap='BuPu')
        ax1.set_ylim(107., 4000.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
        #ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        ax1.set_title("cloud fraction (15min) obs", fontsize=16)#
        ax1.set_xlabel("time [hh:mm]", fontsize=12)
        ax1.set_ylabel("height [m]", fontsize=12)

        ax2 = plt.subplot(2, 1, 2)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.get_xaxis().tick_bottom()
        ax2.get_yaxis().tick_left()
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax2.xaxis_date()
        cax = ax2.pcolormesh(CF_mod_Dataset.time.values, CF_mod_Dataset.height.values, CF_mod_Dataset.CF.values.transpose(), vmin=0., vmax=0.5, cmap='BuPu')
        ax2.set_ylim(107., 4000.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
        #ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        ax2.set_title("cloud fraction (15min) model", fontsize=16)#
        ax2.set_xlabel("time [hh:mm]", fontsize=12)
        ax2.set_ylabel("height [m]", fontsize=12)
        plt.tight_layout()
        plt.savefig(pathDebugFig+'cloudFraction_PBL_mod_obs_'+date+'.png', format='png')

    # opening the file containing all the data
    ##infile = open(filenameObs, 'rb')
    #new_dict = pickle.load(infile, encoding='latin1')

    #cloudMask_obs[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[8]['cloudMask']), time[:], \
                                                         #new_dict[3]['height'], date)
    #cloudMask_mod[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[7]['cloudMask']), time[:], \
                                                          #  new_dict[3]['height'], date)
    #cloudFractionTot_obs[indFile,:,:] = np.asarray(new_dict[8]['totalCloudFraction'])
    #cloudFractionTot_mod[indFile,:,:] = np.asarray(new_dict[7]['totalCloudFraction'])
    #height_obs                            = np.asarray(new_dict[7]['heightCloudFraction'])
    #print(height_mod[0:5])
    #print(height_obs[0:5])
    #streasuka

# concatenating along a new dimension ( number of days)
CFAll_obs_dataset = xr.concat(CF_obs_arr, 'ndays')
CFAll_mod_dataset = xr.concat(CF_mod_arr, 'ndays')
PBLAll_clouds_obs_dataset = xr.concat(PBLclouds_obs, 'ndays')
PBLAll_clouds_mod_dataset = xr.concat(PBLclouds_mod, 'ndays')


# calculating average cloud fraction every 15 min
CFall_mod_mean = CFAll_mod_dataset.mean(dim='ndays', skipna=True)
CFall_obs_mean = CFAll_obs_dataset.mean(dim='ndays', skipna=True)
PBLAll_clouds_obs_mean = PBLAll_clouds_obs_dataset.mean(dim='ndays', skipna=True)
PBLAll_clouds_mod_mean = PBLAll_clouds_mod_dataset.mean(dim='ndays', skipna=True)

# calculating min CB/max CT over days
CB_obs_min = PBLAll_clouds_obs_dataset.cloudBase.min(dim='ndays', skipna=True)
CT_obs_max = PBLAll_clouds_obs_dataset.cloudTop.max(dim='ndays', skipna=True)
CB_mod_min = PBLAll_clouds_mod_dataset.cloudBase.min(dim='ndays', skipna=True)
CT_mod_max = PBLAll_clouds_mod_dataset.cloudTop.max(dim='ndays', skipna=True)

print(CB_obs_min)
# calculating max and min across the levels
CB_obs_minmin = CB_obs_min.min(dim='levels', skipna=True)
CT_obs_maxmax = CT_obs_max.max(dim='levels', skipna=True)
CB_mod_minmin = CB_mod_min.min(dim='levels', skipna=True)
CT_mod_maxmax = CT_mod_max.max(dim='levels', skipna=True)

CBCTTAll_obs_mean = PBLAll_clouds_obs_mean.mean(dim='levels', skipna=True)
CBCTTAll_mod_mean = PBLAll_clouds_mod_mean.mean(dim='levels', skipna=True)
CBCTTAll_obs_std = PBLAll_clouds_obs_mean.std(dim='levels', skipna=True)
CBCTTAll_mod_std = PBLAll_clouds_mod_mean.std(dim='levels', skipna=True)

# calculating running mean of the average time serie
CB_obs_pd = pd.Series(CB_obs_minmin, index=CBCTTAll_obs_mean.time.values)
CB_obs_runMean = CB_obs_pd.rolling('15min').mean()
CB_obs_runStd = CB_obs_pd.rolling('15min').std()
CB_mod_pd = pd.Series(CB_mod_minmin, index=CBCTTAll_mod_mean.time.values)
CB_mod_runMean = CB_mod_pd.rolling('15min').mean()
CB_mod_runStd = CB_mod_pd.rolling('15min').std()

CT_obs_pd = pd.Series(CT_obs_maxmax, index=CBCTTAll_obs_mean.time.values)
CT_obs_runMean = CT_obs_pd.rolling('15min').mean()
CT_obs_runStd = CT_obs_pd.rolling('15min').std()
CT_mod_pd = pd.Series(CT_mod_maxmax, index=CBCTTAll_mod_mean.time.values)
CT_mod_runMean = CT_mod_pd.rolling('15min').mean()
CT_mod_runStd = CT_mod_pd.rolling('15min').std()

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(14, 8))
ax1 = plt.subplot(2, 1, 1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis_date()
y1 = CB_obs_runMean - CB_obs_runStd
y2 = CB_obs_runMean + CB_obs_runStd
ax1.fill_between(CBCTTAll_obs_mean.time.values, y1, y2, where=y2 > y1, facecolor='black', alpha=0.2)
ax1.plot(CBCTTAll_obs_mean.time.values, CB_obs_runMean, color='black')

y11 = CB_mod_runMean - CB_mod_runStd
y22 = CB_mod_runMean + CB_mod_runStd
ax1.fill_between(CBCTTAll_mod_mean.time.values, y11, y22, where=y22 > y11, facecolor='red', alpha=0.2)
ax1.plot(CBCTTAll_mod_mean.time.values, CB_mod_runMean, color='red')
#ax1.plot(CBCTTAll_obs_mean.time.values, CB_obs_runMean+CB_obs_runStd, linestyle=':', color='black')
#ax1.plot(CBCTTAll_obs_mean.time.values, CB_obs_runMean-CB_obs_runStd, linestyle=':', color='black')

#ax1.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudBase.values+CBCTTAll_obs_std.cloudBase.values, linestyle=':', color='black')
#ax1.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudBase.values, color='black', label='obs')
#ax1.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudBase.values-CBCTTAll_obs_std.cloudBase.values, linestyle=':', color='black')
#ax1.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_mod_mean.cloudBase.values+CBCTTAll_mod_std.cloudBase.values, linestyle=':', color='red')
#ax1.plot(CBCTTAll_mod_mean.time.values, CBCTTAll_mod_mean.cloudBase.values, color='red', label='model')
#ax1.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_mod_mean.cloudBase.values-CBCTTAll_mod_std.cloudBase.values, linestyle=':', color='red')
ax1.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
#ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
ax1.set_title("cloud base [m]", fontsize=16)#
ax1.set_xlabel("time [hh:mm]", fontsize=12)
ax1.set_ylabel("height [m]", fontsize=12)
ax1.grid(True, which="both")

ax2 = plt.subplot(2, 1, 2)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis_date()
y1 = CT_obs_runMean - CT_obs_runStd
y2 = CT_obs_runMean + CT_obs_runStd
ax2.fill_between(CBCTTAll_obs_mean.time.values, y1, y2, where=y2 > y1, facecolor='black', alpha=0.2)
ax2.plot(CBCTTAll_obs_mean.time.values, CT_obs_runMean, color='black')

y11 = CT_mod_runMean - CT_mod_runStd
y22 = CT_mod_runMean + CT_mod_runStd
ax2.fill_between(CBCTTAll_mod_mean.time.values, y11, y22, where=y22 > y11, facecolor='red', alpha=0.2)
ax2.plot(CBCTTAll_mod_mean.time.values, CT_mod_runMean, color='red')
#ax2.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudTop.values+CBCTTAll_obs_std.cloudTop.values, linestyle=':', color='black')
#ax2.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudTop.values, color='black', label='obs')
#ax2.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_obs_mean.cloudTop.values-CBCTTAll_obs_std.cloudTop.values, linestyle=':', color='black')
#ax2.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_mod_mean.cloudTop.values+CBCTTAll_mod_std.cloudTop.values, linestyle=':', color='red')
#ax2.plot(CBCTTAll_mod_mean.time.values, CBCTTAll_mod_mean.cloudTop.values, color='red', label='model')
#ax2.plot(CBCTTAll_obs_mean.time.values, CBCTTAll_mod_mean.cloudTop.values-CBCTTAll_mod_std.cloudTop.values, linestyle=':', color='red')
ax2.set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
#ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
ax2.grid(True, which="both")

ax2.set_title("cloud top [m]", fontsize=16)#
ax2.set_xlabel("time [hh:mm]", fontsize=12)
ax2.set_ylabel("height [m]", fontsize=12)
plt.tight_layout()
plt.savefig(pathDebugFig+'cloudTop_Base_meanAll_mod_obs.png', format='png')









Ncols = 1
Nrows = 2
from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(10,6), constrained_layout=True)

# Design your figure properties
widths = [8,1]
gs = GridSpec(Nrows, Ncols + 1, figure=fig, width_ratios=widths)

axes = []
for i in range(Nrows):
    for j in range(Ncols):
        axes.append(fig.add_subplot(gs[i, j]))
        axes[-1].spines["top"].set_visible(False)
        axes[-1].spines["right"].set_visible(False)
        axes[-1].get_xaxis().tick_bottom()
        axes[-1].get_yaxis().tick_left()
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        axes[-1].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        axes[-1].xaxis_date()
        if i == 0:
            im = axes[-1].pcolormesh(CFall_obs_mean.time.values, CFall_obs_mean.height.values, CFall_obs_mean.CF.values.transpose(), vmin=0., vmax=0.2, cmap='BuPu')
            axes[-1].set_title("cloud fraction (15min) obs", fontsize=12)#
        if i == 1:
            im = axes[-1].pcolormesh(CFall_mod_mean.time.values, CFall_mod_mean.height.values, CFall_mod_mean.CF.values.transpose(), vmin=0., vmax=0.2, cmap='BuPu')
            axes[-1].set_title("cloud fraction (15min) mod", fontsize=12)#
#cax = ax1.pcolormesh(CFall_obs_mean.time.values, CFall_obs_mean.height.values, CFall_obs_mean.CF.values.transpose(), vmin=0., vmax=0.5, cmap='BuPu')
        axes[-1].set_ylim(107., 3500.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
#ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
        axes[-1].set_xlabel("time [hh:mm]", fontsize=12)
        axes[-1].set_ylabel("height [m]", fontsize=12)
        axes[-1].grid(True, which="both")

# Shared colorbar
axes.append(fig.add_subplot(gs[:, Ncols]))
cbar = fig.colorbar(im, cax=axes[-1], shrink=5.)
cbar.set_label(label="cloud fraction",size=12)
cbar.ax.tick_params(labelsize=12)
#cbar.aspect=40
plt.savefig(pathDebugFig+'cloudFraction_PBL_mod_obs_AllDays.png', format='png')

##cbar = fig.colorbar(cax, orientation='vertical')
#cbar.set_label(label="cloud fraction",size=16)
#cbar.ax.tick_params(labelsize=12)
#cbar.aspect=80



#ax2 = plt.subplot(2, 1, 2)
#ax2.spines["top"].set_visible(False)
#ax2.spines["right"].set_visible(False)
#ax2.get_xaxis().tick_bottom()
#ax2.get_yaxis().tick_left()
#ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#ax2.xaxis_date()
#cax = ax2.pcolormesh(CFall_mod_mean.time.values, CFall_mod_mean.height.values, CFall_mod_mean.CF.values.transpose(), vmin=0., vmax=0.5, cmap='BuPu')
#ax2.set_ylim(107., 4000.)  # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
#ax1.set_xlim(timeStart, timeEnd)  # limits of the x-axes
#ax2.set_title("cloud fraction (15min) model", fontsize=14)#
#ax2.set_xlabel("time [hh:mm]", fontsize=12)
#ax2.set_ylabel("height [m]", fontsize=12)
#cbar = fig.colorbar(cax, ax=ax.flat)
#cbar.set_label(label="cloud fraction",size=16)
#cbar.ax.tick_params(labelsize=12)
#cbar.aspect=80
#plt.tight_layout()

#print(ax.flat)
#fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
#for ax in axes.flat:
#    im = ax.pcolormesh(np.random.random((10,10)), vmin=0, vmax=1)#

#fig.colorbar(im, ax=axes.flat)
#plt.show()



#def f_calculateCloudFractionPBLclouds(PBLcloudmask,time,height,Nmin_string):
#    '''
#    author: Claudia Acquistapace
#    date: friday 19 June 2020
#    goal: calculate cloud fraction from cloud mask over a given amount of minutes Nmin
#    input: PBLcloudmask: cloud mask matrix
#            time: time array corresponding to the cloud mask
#            height: height array corresponding to the cloud mask
#            Nmin: number of minutes over which to calculate the cloud fraction
#    output: CF_dataset = xr.Dataset({'CF': (['time','height'], CF_PBL)},
#                        coords = {'time':datetime_CF,
#                                  'height':height})
#    '''
#    # calculating cloud fraction every 15 minutes
#    cloudMask_DF = pd.DataFrame(PBLcloudmask, index=time, columns=height)
##    datetime_CF = pd.date_range(start=PBLtime[0], end=PBLtime[-1], freq=Nmin_string+'min')
#    datetime_CF = datetime_CF.to_pydatetime()

#    CF_PBL = np.zeros((len(datetime_CF), len(height)))
#    for indTime in range(len(datetime_CF)-1):
#        mask_t = (cloudMask_DF.index > datetime_CF[indTime]) * (cloudMask_DF.index < datetime_CF[indTime+1])
#        Selection_cloudMask = cloudMask_DF[mask_t]
#        for indHeight in range(len(height)):
#            CFArray = Selection_cloudMask.loc[:,Selection_cloudMask.columns[indHeight]]
#            CF_PBL[indTime,indHeight] = len(CFArray[CFArray == 1])/len(CFArray)

#    CF_dataset = xr.Dataset({'CF': (['time','height'], CF_PBL)},
#                        coords = {'time':datetime_CF,
#                                  'height':height})
#    return(CF_dataset)