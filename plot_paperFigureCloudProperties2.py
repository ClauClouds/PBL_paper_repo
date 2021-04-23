#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:00:24 2020
@date; 19/05/2020
@author: Claudia Acquistapace
@goal: produce mean time series of cloud base, cloud top, cloud thickness, for the entire dataset
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



date_arr                        = []
datetime_out                    = []
PBLcloud_obs                    = []
PBLcloud_mod                    = []
clouds_obs                      = []
clouds_mod                      = []

# surface fluxes
NprofilesOut = 48 # half hourly data
SHF_mod_matrix                      = np.zeros((Nfiles,NprofilesOut))
SHF_obs_matrix                      = np.zeros((Nfiles,NprofilesOut))
LHF_mod_matrix                      = np.zeros((Nfiles,NprofilesOut))
LHF_obs_matrix                      = np.zeros((Nfiles,NprofilesOut))
datetime_out = pd.date_range(start=datetime.datetime(2000, 1, 1,0,0,0), end=datetime.datetime(2000, 1, 1,23,59,59), freq='30min')

#%%

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
    # w_mod         = ncdata.groups['Temp_data'].variables['vertWind']
    # varW_mod      = ncdata.groups['Temp_data'].variables['varianceW']

    # opening the file containing all the data
    infile = open(pathObs + 'dictionaries_ModObs_' + date + '.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')

    # fluxes
    SHF_mod = new_dict[10]['SHF_iconlem']
    SHF_obs = new_dict[10]['SHF_obs']
    LHF_mod = new_dict[10]['LHF_iconlem']
    LHF_obs = new_dict[10]['LHF_obs']
    timeFluxes = new_dict[10]['datetime_30m']
    SHF_mod_matrix[indFile,:] = SHF_mod[:-1]
    SHF_obs_matrix[indFile,:] = SHF_obs
    LHF_mod_matrix[indFile,:] = LHF_mod[:-1]
    LHF_obs_matrix[indFile,:] = LHF_obs
    



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

#%%
# calculating mean CB and CT for each time, averaging together all days
print('calculating CB and CT for observations')
CBarr_obs, CTarr_obs, TKarr_obs =  f_calculateMinCloudBaseTop(clouds_obs,PBLcloud_obs, date_arr)
print('calculating CB and CT for model')
CBarr_mod, CTarr_mod, TKarr_mod  =  f_calculateMinCloudBaseTop(clouds_mod,PBLcloud_mod, date_arr)

#averaging cloud bases of all days together
meanCBobs = np.nanmean(CBarr_obs, axis=1)
stdCBobs  = np.nanstd(CBarr_obs, axis=1)
meanCTobs = np.nanmean(CTarr_obs, axis=1)
stdCTobs  = np.nanstd(CTarr_obs, axis=1)
meanCBmod = np.nanmean(CBarr_mod, axis=1)
stdCBmod  = np.nanstd(CBarr_mod, axis=1)
meanCTmod = np.nanmean(CTarr_mod, axis=1)
stdCTmod  = np.nanstd(CTarr_mod, axis=1)
meanThickness_mod = np.nanmean(TKarr_mod, axis=1)
stdThickness_mod = np.nanstd(TKarr_mod, axis=1)
meanThickness_obs = np.nanmean(TKarr_obs, axis=1)
stdThickness_obs = np.nanstd(TKarr_obs, axis=1)
MeanSHF_mod = np.nanmean(SHF_mod_matrix, axis=0)
MeanSHF_obs = np.nanmean(SHF_obs_matrix, axis=0)
MeanLHF_mod = np.nanmean(LHF_mod_matrix, axis=0)
MeanLHF_obs = np.nanmean(LHF_obs_matrix, axis=0)

print(np.shape(meanCBobs))

#%%
# plotting
meanArr_mod = [meanCBmod, meanCTmod, meanThickness_mod, -MeanSHF_mod, -MeanLHF_mod]
meanArr_obs = [meanCBobs, meanCTobs, meanThickness_obs, MeanSHF_obs, MeanLHF_obs]
YlabelArr = ['cloud base [m]', 'cloud top [m]', 'geometrical thickness [m]', 'sensible heat flux', 'latent heat flux']
yminArr = [700., 0., 0., -50., 0.]
ymaxArr = [100., 2000., 2000., 300., 250.]
Nrows = 3
Ncols = 1
Nplots = 3
fontSizeTitle = 16
fontSizeX     = 14
fontSizeY     = 14
labelsizeaxes = 12
fs            = 12
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
timeStandard = pd.date_range(start=datetime.datetime(2000,1,1,0,0,0), \
                                    end=datetime.datetime(2000,1,1,23,59,59), freq='9s')
fig, axes       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(8,10))
#matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
axes[0] = plt.subplot(Nrows, Ncols, 1)
axes[0].set_ylabel('Cloud base [m] - PBL clouds -', fontsize=fs)
axes[0].set_xlabel('time [hh:mm]', fontsize=fs)
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
axes[0].get_xaxis().tick_bottom()
axes[0].get_yaxis().tick_left()
axes[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axes[0].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
axes[0].xaxis_date()
axes[0].plot(timeStandard, meanCBobs, label='mean obs', color='black')
#axes[0].plot(timeStandard, meanCBobs-stdCBobs, linestyle=':', linewidth=0.1, color='black')
#axes[0].plot(timeStandard, meanCBobs+stdCBobs, linestyle=':', linewidth=0.1, color='black')
#axes[0].plot(timeStandard, meanCBmod-stdCBmod, linestyle=':', linewidth=0.1, color='red')
#axes[0].plot(timeStandard, meanCBmod+stdCBmod, linestyle=':', linewidth=0.1, color='red')
axes[0].plot(timeStandard, meanCBmod, label='mean mod', color='red')
axes[0].set_xlim(datetime.datetime(2000, 1, 1,6,0,0), datetime.datetime(2000, 1, 1,23,59,59))
axes[0].grid(True, which="both")
axes[0].set_ylim(500., 2500.)

axes[1] = plt.subplot(Nrows, Ncols, 2)
axes[1].set_ylabel('Cloud top [m] - PBL clouds -', fontsize=fs)
axes[1].set_xlabel('time [hh:mm]', fontsize=fs)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
axes[1].get_xaxis().tick_bottom()
axes[1].get_yaxis().tick_left()
axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axes[1].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
axes[1].xaxis_date()
#axes[1].plot(timeStandard, meanCTobs-stdCTobs, linestyle=':', linewidth=0.1, color='black')
#axes[1].plot(timeStandard, meanCTobs+stdCTobs, linestyle=':', linewidth=0.1, color='black')
#axes[1].plot(timeStandard, meanCTmod-stdCTmod, linestyle=':', linewidth=0.1, color='red')
#axes[1].plot(timeStandard, meanCTmod+stdCTmod, linestyle=':', linewidth=0.1, color='red')
axes[1].plot(timeStandard, meanCTobs, label='mean obs', color='black')
axes[1].plot(timeStandard, meanCTmod, label='mean mod', color='red')
axes[1].set_xlim(datetime.datetime(2000, 1, 1,6,0,0), datetime.datetime(2000, 1, 1,23,59,59))
axes[1].grid(True, which="both")
axes[1].set_ylim(500., 3500.)

axes[2] = plt.subplot(Nrows, Ncols, 3)
axes[2].set_ylabel('Cloud thickness [m] - PBL clouds -', fontsize=fs)
axes[2].set_xlabel('time [hh:mm]', fontsize=fs)
axes[2].spines["top"].set_visible(False)
axes[2].spines["right"].set_visible(False)
axes[2].get_xaxis().tick_bottom()
axes[2].get_yaxis().tick_left()
axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axes[2].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
axes[2].xaxis_date()
#axes[2].plot(timeStandard, meanThickness_obs-stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_obs+stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_mod-stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
#axes[2].plot(timeStandard, meanThickness_mod+stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
axes[2].plot(timeStandard, meanThickness_obs, label='mean obs', color='black')
axes[2].plot(timeStandard, meanThickness_mod, label='mean mod', color='red')
axes[2].set_xlim(datetime.datetime(2000, 1, 1,6,0,0), datetime.datetime(2000, 1, 1,23,59,59))
axes[2].grid(True, which="both")
axes[2].set_ylim(0., 1500.)


#axes[3] = plt.subplot(Nrows, Ncols, 4)
#axes[3].set_ylabel('sensible heat flux [w/m^2]', fontsize=fs)
#axes[3].set_xlabel('time [hh:mm]', fontsize=fs)
#axes[3].spines["top"].set_visible(False)
#axes[3].spines["right"].set_visible(False)
#axes[3].get_xaxis().tick_bottom()
#axes[3].get_yaxis().tick_left()
#axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#axes[3].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#axes[3].xaxis_date()
#axes[2].plot(timeStandard, meanThickness_obs-stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_obs+stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_mod-stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
#axes[2].plot(timeStandard, meanThickness_mod+stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
#axes[3].plot(timeFluxes[:-1], MeanSHF_obs, label='mean obs', color='black')
#axes[3].plot(timeFluxes[:-1], -MeanSHF_mod, label='mean mod', color='red')
#axes[3].set_xlim(datetime.datetime(2013, 5, 28,6,0,0), datetime.datetime(2013, 5, 28,23,59,59))
#axes[3].grid(True, which="both")
#axes[3].set_ylim(-50., 250.)


#axes[4] = plt.subplot(Nrows, Ncols, 5)
#axes[4].set_ylabel('latent heat flux [w/m^2]', fontsize=fs)
#axes[4].set_xlabel('time [hh:mm]', fontsize=fs)
#axes[4].spines["top"].set_visible(False)
#axes[4].spines["right"].set_visible(False)
#axes[4].get_xaxis().tick_bottom()
#axes[4].get_yaxis().tick_left()
#axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#axes[4].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#axes[4].xaxis_date()
#axes[2].plot(timeStandard, meanThickness_obs-stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_obs+stdThickness_obs, linestyle=':', linewidth=0.1, color='black')
#axes[2].plot(timeStandard, meanThickness_mod-stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
#axes[2].plot(timeStandard, meanThickness_mod+stdThickness_mod, linestyle=':', linewidth=0.1, color='red')
#axes[4].plot(timeFluxes[:-1], MeanLHF_obs, label='mean obs', color='black')
#axes[4].plot(timeFluxes[:-1], -MeanLHF_mod, label='mean mod', color='red')
#axes[4].set_xlim(datetime.datetime(2013, 5, 28,6,0,0), datetime.datetime(2013, 5, 28,23,59,59))
#axes[4].grid(True, which="both")
#axes[4].set_ylim(-50., 250.)

plt.tight_layout()
plt.savefig(pathFig+'figure8_cloudProperties_vs_SHF.png', format='png')

