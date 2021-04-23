#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:57:03 2020

@author: cacquist
@goal : determine cloud generated by surface induced heating/motions via LCL.
Method: the idea is to calculate LCL height for every time stamp based on pressure, 
humidity and temperature at the surface and then identify as PBL clouds only those clouds 
whose cloud base is located within z_lcl+-x%z_lcl
The method is applied between 7:00 and 18:00 
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
from myFunctions import lcl
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']


# directories where data are stored
patch = 'patch003'
Nlevels_example = np.arange(8)
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathOut                         = pathObs
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/debugging/'
path_tower          = '/data/hatpro/jue/hdcp2/meteo_data/'#'/data/TR32/D2/data/juelich/met_mast/'
path_mwr_joyce      = '/data/hatpro/hps/data/level2/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListObs)
date_arr                        = []
percentage                      = 0.5 # percentage expressed as a number between 0 and 1, for searching for PBl clouds above and below Z_LCL
flagPlot                        = 1
# reading LCL ncdf files with LCL arrays for model and all types of observations
LCL_nc = nc4.Dataset(pathObs+'LCL_Mod_obs_allDataset.nc', mode='r')  # Or from siphon.ncss
LCL_dataset = xr.open_dataset(xr.backends.NetCDF4DataStore(LCL_nc))
LCL_dates = LCL_dataset.Nfiles.values


for indFile in range(Nfiles):
    print(indFile)
    # defining reference timearray for resizing cloud bases/tops/thicknesses
    date = fileListObs[indFile][81:89]
    print('processing date ' + date)

    filenameMod = pathMod+'icon_lem_derivedproperties'+date+'.nc'
    filenameObs = pathObs+'dictionaries_ModObs_'+date+'.p'
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)

    # slicing in time the xarray dataset between Tstart and Tend
    Tstart  = datetime.datetime(yy,5,28,7,0,0)
    Tend    = datetime.datetime(yy,5,28,18,0,0)

    # reading corresponding LCL height for model and obs corresponding to the selected day
    indDate = np.where(LCL_dates == date)
    LCL_mod = LCL_dataset.LCL_ICONLEM.sel(time=slice(Tstart, Tend))
    LCL_obs = LCL_dataset.LCL_EC.sel(time=slice(Tstart, Tend))#.values[:,indDate]
    PlotTime = LCL_dataset.time.sel(time=slice(Tstart, Tend))

    # calculating range of heights for accepting clouds as PBL clouds
    LCL_day_mod = LCL_mod.values[:,indDate][:,0,0]
    LCL_day_obs = LCL_obs.values[:,indDate][:,0,0]
    #print('max and min uncertainty for LCL height')
    print(np.nanmin(percentage * LCL_day_mod))
    print(np.nanmax(percentage * LCL_day_mod))

    top_mod    = LCL_day_mod + percentage * LCL_day_mod
    bottom_mod = LCL_day_mod - percentage * LCL_day_mod
    top_obs    = LCL_day_obs + percentage * LCL_day_obs
    bottom_obs = LCL_day_obs - percentage * LCL_day_obs


    # reading CB and CT files
    timeReferenceFormat = pd.date_range(date, periods=9600, freq='9s')
    Dataset_example = xr.Dataset({'cloudBase': (['time', 'levels'], np.random.rand(len(timeReferenceFormat),len(Nlevels_example)))},
                                 coords={'time'  : timeReferenceFormat,
                                         'levels': Nlevels_example})


    # read xarray Dataset from nc file and convert it from NDF4.Dataset to xarray.dataset type for allowing concatenation
    clouds_mod_nc = nc4.Dataset(pathMod+'Clouds_iconlem_'+date+'.nc', mode='r')  # Or from siphon.ncss
    cloudDay_mod = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_mod_nc))
    cloudStandard_mod = cloudDay_mod.reindex({"time": timeReferenceFormat}, copy=True)
    cloudStandard_mod_reshaped = cloudStandard_mod.interp_like(Dataset_example)

    clouds_obs_nc = nc4.Dataset(pathMod+'Clouds_Obs_'+date+'.nc', mode='r')  # Or from siphon.ncss
    cloudDay_obs = xr.open_dataset(xr.backends.NetCDF4DataStore(clouds_obs_nc))
    cloudStandard_obs = cloudDay_obs.reindex({"time": timeReferenceFormat}, copy=True)
    cloudStandard_obs_reshaped = cloudStandard_obs.interp_like(Dataset_example)

    # slicing between time min and time max for the cloud dataset
    Tstart_ref = datetime.datetime(yy,mm,dd,7,0,0)
    Tend_ref = datetime.datetime(yy,mm,dd,18,0,0)
    Convective_clouds_mod = cloudStandard_mod_reshaped.sel(time=slice(Tstart_ref, Tend_ref))
    Convective_clouds_obs = cloudStandard_obs_reshaped.sel(time=slice(Tstart_ref, Tend_ref))

    PBL_cloud_base_mod = np.zeros((len(Convective_clouds_obs.time.values),8))
    PBL_cloud_base_mod.fill(np.nan)
    PBL_cloud_top_mod = np.zeros((len(Convective_clouds_obs.time.values), 8))
    PBL_cloud_top_mod.fill(np.nan)
    PBL_cloud_thickness_mod = np.zeros((len(Convective_clouds_obs.time.values), 8))
    PBL_cloud_thickness_mod.fill(np.nan)
    PBL_cloud_base_obs = np.zeros((len(Convective_clouds_obs.time.values),8))
    PBL_cloud_base_obs.fill(np.nan)
    PBL_cloud_top_obs = np.zeros((len(Convective_clouds_obs.time.values), 8))
    PBL_cloud_top_obs.fill(np.nan)
    PBL_cloud_thickness_obs = np.zeros((len(Convective_clouds_obs.time.values), 8))
    PBL_cloud_thickness_obs.fill(np.nan)

    for indTime in range(len(Convective_clouds_mod.time.values)):
        for indLev in range(8):
            cb_value_mod = Convective_clouds_mod.cloudBase.values[indTime, indLev]
            ct_value_mod = Convective_clouds_mod.cloudTop.values[indTime, indLev]
            cloudThickness_value_mod = Convective_clouds_mod.cloudThick.values[indTime, indLev]
            cb_value_obs = Convective_clouds_obs.cloudBase.values[indTime, indLev]
            ct_value_obs = Convective_clouds_obs.cloudTop.values[indTime, indLev]
            cloudThickness_value_obs = Convective_clouds_obs.cloudThick.values[indTime, indLev]
            if (cb_value_mod > bottom_mod[indTime]) * (cb_value_mod < top_mod[indTime]):
                PBL_cloud_base_mod[indTime, indLev] = cb_value_mod
                PBL_cloud_top_mod[indTime, indLev] = ct_value_mod
                PBL_cloud_thickness_mod[indTime, indLev] = cloudThickness_value_mod
            if (cb_value_obs > bottom_obs[indTime]) * (cb_value_obs < top_obs[indTime]):
                PBL_cloud_base_obs[indTime, indLev] = cb_value_obs
                PBL_cloud_top_obs[indTime, indLev] = ct_value_obs
                PBL_cloud_thickness_obs[indTime, indLev] = cloudThickness_value_obs



    # define PBLcloud datasets and save them as ncdf files.
    PBL_cloud_Dataset_mod = xr.Dataset({'cloudBase': (['time', 'levels'], PBL_cloud_base_mod),
                                         'cloudTop': (['time', 'levels'], PBL_cloud_top_mod),
                                        'cloudThick': (['time', 'levels'], PBL_cloud_thickness_mod)},
                                 coords={'time'  : PlotTime,
                                         'levels': Nlevels_example})
    PBL_cloud_Dataset_obs = xr.Dataset({'cloudBase': (['time', 'levels'], PBL_cloud_base_obs),
                                         'cloudTop': (['time', 'levels'], PBL_cloud_top_obs),
                                        'cloudThick': (['time', 'levels'], PBL_cloud_thickness_obs)},
                                 coords={'time'  : PlotTime,
                                         'levels': Nlevels_example})

    PBL_cloud_Dataset_obs.to_netcdf(pathOut+'PBLClouds_Obs_'+date+'.nc')
    PBL_cloud_Dataset_mod.to_netcdf(pathOut+'PBLClouds_iconlem_'+date+'.nc')

    # plotting area for cloud base search
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8))
    fontSizeTitle = 12
    fontSizeX = 10
    fontSizeY = 10
    fontSizeCbar = 10
    labelsizeaxes = 10
    from matplotlib import rcParams
    plt.gcf().subplots_adjust(bottom=0.15)
    matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Tahoma']
    ax = plt.subplot(2,1,1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylim(107., 4000.)
    # ax[0].set_xlim(-100., 450.)
    # ax[0].set_ylim(-100., 450.)
    ax.set_xlabel('time [hh:mm]', fontsize=16)
    ax.set_ylabel('height [m]', fontsize=16)
    ax.set_title('model')
    ax.grid(b=True, which='major', color='#666666', linestyle=':')
    ax.plot(PlotTime, top_mod[:], linestyle=':', label='top LCL mod', color='blue')
    ax.plot(PlotTime, bottom_mod[:], linestyle=':', label='bottom LCL mod', color='blue')
    ax.plot(PlotTime, LCL_day_mod[:], linestyle=':', label='LCL', color='red')

    for indLev in range(8):
        ax.plot(PlotTime, Convective_clouds_mod.cloudBase.values[:,indLev], 'o', label='cb all mod', color=plt.cm.cool(indLev))
        #ax.plot(PlotTime, PBLclouds.cloudBase.values[:,indLev], 'o', label='cb', color=plt.cm.viridis(indLev))
        ax.plot(PlotTime, PBL_cloud_base_mod[:,indLev], marker=6, label='cb mod', color=plt.cm.seismic(indLev))
        ax.plot(PlotTime, PBL_cloud_top_mod[:,indLev], marker=11, label='ct mod', color=plt.cm.inferno(indLev))
    if indLev == 0:
        ax.legend(frameon=False)

    ax2 = plt.subplot(2,1,2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    ax2.set_ylim(107., 4000.)

    # ax[0].set_xlim(-100., 450.)
    # ax[0].set_ylim(-100., 450.)
    ax2.set_xlabel('time [hh:mm]', fontsize=16)
    ax2.set_ylabel('height [m]', fontsize=16)
    ax.set_title('obs')

    ax2.grid(b=True, which='major', color='#666666', linestyle=':')
    ax2.plot(PlotTime, top_obs[:], linestyle=':', label='top LCL obs', color='blue')
    ax2.plot(PlotTime, bottom_obs[:], linestyle=':', label='bottom LCL obs', color='blue')
    ax2.plot(PlotTime, LCL_day_obs[:], linestyle=':', label='LCL obs', color='red')

    for indLev in range(8):
        ax2.plot(PlotTime, Convective_clouds_obs.cloudBase.values[:,indLev], 'o', label='cb all obs', color=plt.cm.cool(indLev))
        #ax2.plot(PlotTime, PBLclouds.cloudBase.values[:,indLev], 'o', label='cb', color=plt.cm.viridis(indLev))
        ax2.plot(PlotTime, PBL_cloud_base_obs[:,indLev], marker=6, label='cb obs', color=plt.cm.seismic(indLev))
        ax2.plot(PlotTime, PBL_cloud_top_obs[:,indLev], marker=11, label='ct obs', color=plt.cm.inferno(indLev))
    if indLev == 0:
        ax2.legend(frameon=False)

    fig.tight_layout()
    plt.savefig(pathFig + date + '_LCL_edges.png', format='png')


