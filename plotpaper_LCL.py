#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:13:35 2020

@author: cacquist
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
patch                           = 'patch003'

flagThermodynVarScatterPlots  = 1
# thermodynamic variables
thetaV_radiosObs                = []
thetaV_iconlem                  = []
thetaV_mwrObs                   = []
rh_radiosObs                    = []
rh_mod                          = []
T_mod                           = []
P_radiosondes                   = []
T_radiosondes                   = []
theta_v_radiosondes             = []
height_radiosondes              = []
time_radiosondes                = []
theta_v_mod                     = []
date_arr                        = []
time_mod                        = []
height_mod                      = []
lcl_radiosondes                 = []
ccl_radiosondes                 = []
lts_radiosondes                 = []
pblHeight_radiosondes           = []
lcl_mod                         = []
ccl_mod                         = []
lts_mod                         = []
pblHeightRN_mod                 = []
pblHeightTW_mod                 = []
pblHeightRN_windLidarObs        = []
pblHeightTW_windLidarObs        = []


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


for indFile in range(Nfiles):

    print(indFile)
    date = fileListObs[indFile][81:89]
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)

    print('processing date ' + date)
    ncdata = Dataset(pathMod+'icon_lem_derivedproperties'+date+'.nc', mode='r')
    time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]
    
    # opening the file containing all the data
    infile = open(pathObs+'dictionaries_ModObs_'+date+'.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')
    
    # reading data radiosoundings
    from myFunctions import f_reshapeRadiosondes
    RadiosondeFormatted = f_reshapeRadiosondes(new_dict)
     
    P_radiosondes.append(RadiosondeFormatted['P'])
    T_radiosondes.append(RadiosondeFormatted['T'])
    theta_v_radiosondes.append(RadiosondeFormatted['theta_v'])
    height_radiosondes.append(RadiosondeFormatted['height'])
    time_radiosondes.append(RadiosondeFormatted['time'])
    lcl_radiosondes.append(RadiosondeFormatted['lcl'])
    ccl_radiosondes.append(RadiosondeFormatted['ccl'])
    lts_radiosondes.append(RadiosondeFormatted['lts'])
    pblHeight_radiosondes.append(RadiosondeFormatted['pblHeight'])
    rh_radiosObs.append(RadiosondeFormatted['RH'])
    
    # reading model data from lem for fluxes at surface and atm indeces
    theta_v_mod.append(new_dict[5]['virtualPotentialTemperature'])
    time_mod.append(new_dict[5]['time'])
    height_mod.append(new_dict[5]['height'])
    lcl_mod.append(new_dict[5]['lclHeight'])
    ccl_mod.append(new_dict[5]['cclHeight'])
    lts_mod.append(new_dict[5]['LTS'])
    pblHeightTW_mod.append(new_dict[9]['PBLHeightTW'])
    pblHeightRN_mod.append(new_dict[9]['PBLHeightRN'])
    rh_mod.append(new_dict[5]['relativeHumidity'])
    T_mod.append(new_dict[9]['T_iconlem'])



# =============================================================================
# calculating and plotting potential temperature , temperature and relative humidity profiles 
# =============================================================================
from myFunctions import f_calculateMeanThetaVModelProfiles
theta_v_dict_obs_mod_arr = f_calculateMeanThetaVModelProfiles(time_radiosondes, \
                                                              theta_v_radiosondes,\
                                                              T_radiosondes, \
                                                              rh_radiosObs, \
                                                              height_radiosondes, \
                                                              lcl_radiosondes, \
                                                              ccl_radiosondes, \
                                                              lts_radiosondes, \
                                                              pblHeight_radiosondes, \
                                                              theta_v_mod, \
                                                              T_mod, \
                                                              rh_mod, \
                                                              time_mod, \
                                                              height_mod, \
                                                              lcl_mod, \
                                                              ccl_mod, \
                                                              lts_mod, \
                                                              pblHeightRN_mod)






from myFunctions import f_calculateMeanProfilesPlotThetaVRadiosondes
result = f_calculateMeanProfilesPlotThetaVRadiosondes(theta_v_dict_obs_mod_arr, height_mod)

MatrixHourMeanProfileThetaRad = result[0]
MatrixHourStdProfileThetaRad  = result[1]
listHourDict                  = result[2]
MatrixHourMeanProfileTRad     = result[3]
MatrixHourStdProfileTRad      = result[4]
MatrixHourMeanProfileRHRad    = result[5]
MatrixHourStdProfileRHRad     = result[6]
gridHeight                    = height_mod[0]

#%%
# =============================================================================
# calculating and plotting Thermdynamic variables (PBL height, LCL, LTS)
# =============================================================================
print('calculating scatter plots of LCL, LTS, PBL height for the whole dataset')

hourList = []
LCLobs   = []
LCLmod   = []
PBLobs   = []
PBLmod   = []
LTSobs   = []
LTSmod   = []
CCLobs   = []
CCLmod   = []
for ind in range(len(listHourDict)):
    hour = listHourDict[ind]['hour']
    dim = listHourDict[ind]['n_lcl_hour']
    hourList.append(np.repeat(hour, dim))
    LCLobs.append(listHourDict[ind]['lcl_rad_hour'])
    LCLmod.append(listHourDict[ind]['lcl_mod_hour'])
    CCLobs.append(listHourDict[ind]['ccl_rad_hour'])
    CCLmod.append(listHourDict[ind]['ccl_mod_hour'])
    PBLobs.append(listHourDict[ind]['pblHeight_rad_hour'])
    PBLmod.append(listHourDict[ind]['pblHeight_mod_hour'])
    LTSobs.append(listHourDict[ind]['lts_rad_hour'])
    LTSmod.append(listHourDict[ind]['lts_mod_hour'])
    print(hour)
    
hourList = np.concatenate(hourList)
LCLobs = np.concatenate(LCLobs)
LCLmod = np.concatenate(LCLmod)
CCLobs = np.concatenate(CCLobs)
CCLmod = np.concatenate(CCLmod)
LTSobs = np.concatenate(LTSobs)
LTSmod = np.concatenate(LTSmod)
PBLobs = np.concatenate(PBLobs)
PBLmod = np.concatenate(PBLmod)
dti = pd.date_range('2018-01-01', periods=24, freq='H')
data = np.arange(-1000, 2000)
print(CCLobs)
print(CCLmod)

#%%
m, b = np.polyfit(LCLobs, LCLmod, 1)
data = np.arange(2500)
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8,12))
from matplotlib import rcParams
rcParams['font.sans-serif'] = ['Tahoma']
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(2,1,1)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = hourList
plt.plot(data, data, color='black', linestyle=':')
plt.plot(data, m*data + b, color='red')
plt.xlim(0., 2150.)
plt.ylim(0., 2150.)
plt.xlabel('Lifting condensation level [m] - obs ', fontsize=16)
plt.ylabel('Lifting condensation level [m] - mod ', fontsize=16)
cmap = plt.cm.get_cmap('jet', 19) 
cax = ax.scatter(LCLobs, LCLmod, c=colors, cmap=cmap, s=200, vmin=5, vmax=23., edgecolors='black')
#plt.grid(b=True, which='major', color='#666666', linestyle=':')
cbar= fig.colorbar(cax, cmap=cmap)
cbar.set_label(label='time [hh:mm]', size=15, family='Tahoma')

ax = plt.subplot(2,1,2)
m, b = np.polyfit(CCLobs, CCLmod, 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
colors = hourList
plt.plot(data, data, color='black', linestyle=':')
plt.plot(data, m*data + b, color='red')
plt.xlim(0., 4150.)
plt.ylim(0., 4150.)
plt.xlabel('Convective condensation level [m] - obs ', fontsize=16)
plt.ylabel('Convective condensation level [m] - mod ', fontsize=16)
cmap = plt.cm.get_cmap('jet', 19)
cax = ax.scatter(CCLobs, CCLmod, c=colors, cmap=cmap, s=200, vmin=5, vmax=23., edgecolors='black')
#plt.grid(b=True, which='major', color='#666666', linestyle=':')
cbar= fig.colorbar(cax, cmap=cmap)
cbar.set_label(label='time [hh:mm]', size=15, family='Tahoma')
plt.tight_layout()
plt.savefig(pathFig+'LCL_CCL_scatterplot_obs_mod_allDays.png', format='png')




#%%
m, b = np.polyfit(PBLobs, PBLmod, 1)

data = np.arange(2500)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
from matplotlib import rcParams
rcParams['font.sans-serif'] = ['Tahoma']
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = hourList
plt.plot(data, data, color='black', linestyle=':')
plt.plot(data, m*data + b, color='red')
plt.xlim(0., 2450.)
plt.ylim(0., 2450.)
plt.xlabel('PBL height obs [m]', fontsize=16)
plt.ylabel('PBL height mod [m]', fontsize=16)
cmap = plt.cm.get_cmap('jet', 19) 
cax = ax.scatter(PBLobs, PBLmod, c=colors, cmap=cmap, s=200, vmin=5, vmax=23., edgecolors='black')
#plt.grid(b=True, which='major', color='#666666', linestyle=':')
cbar= fig.colorbar(cax, cmap=cmap)
cbar.set_label(label='time [hh:mm]', size=15, family='Tahoma')
plt.tight_layout()
plt.savefig(pathFig+'PBL_scatterplot_obs_mod_allDays.png', format='png')


#%%
if flagThermodynVarScatterPlots == 1:
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    matplotlib.rcParams['savefig.dpi'] = 100
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ax = plt.subplot(1,1,1)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    colors = hourList
    plt.plot(data, data, color='black', linestyle=':')
    plt.xlim(0., 30.)
    plt.ylim(0., 30.)
    plt.xlabel('LTS obs [m]', fontsize=16)
    plt.ylabel('LTS mod [m]', fontsize=16)
    cmap = plt.cm.get_cmap('bwr', 19) 
    cax = ax.scatter(LTSobs, LTSmod, c=colors, cmap=cmap, s=100, vmin=5, vmax=23.)
    fig.colorbar(cax, label='time [hh]', cmap=cmap)
    plt.tight_layout()
    plt.savefig(pathFig+'lts_scatterplot_obs_mod.png', format='png')
    
    
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    matplotlib.rcParams['savefig.dpi'] = 100
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ax = plt.subplot(1,1,1)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    colors = hourList
    plt.plot(data, data, color='black', linestyle=':')
    plt.xlim(0., 2400.)
    plt.ylim(0., 2400.)
    plt.xlabel('PBL height obs [m]', fontsize=16)
    plt.ylabel('PBL height mod [m]', fontsize=16)
    cmap = plt.cm.get_cmap('bwr', 19) 
    cax = ax.scatter(PBLobs, PBLmod, c=colors, cmap=cmap, s=100, vmin=5, vmax=23.)
    fig.colorbar(cax, label='time [hh]', cmap=cmap)
    plt.tight_layout()
    plt.savefig(pathFig+'pblh_scatterplot_obs_mod.png', format='png')
        