#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 13:44:44 2019

@author: cacquist
"""



# ----- importing libraries needed
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
import math
    
from f_processModelOutput import f_processModelOutput
from myFunctions2 import f_readingTowerData
from myFunctions import f_resamplingfield
from myFunctions2 import hourDecimal_to_datetime
from myFunctions import f_resampling_twoD_Field
from myFunctions import f_downscaleScalarfield
from myFunctions import f_downscalevectorfield
#from myFunctions import f_processRadiosondesDay
from myFunctions import f_calcThermodynamics
from cloudnetFunctions import f_calculateCloudMaskCloudnet
from myFunctions import f_resamplingMatrixCloudnet
from myFunctions import f_calculateAllCloudQuantities
from myFunctions import f_selectingPBLcloudWindow


try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
pathDebugFig       = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/20130502/'
date               = '20130502'
filename           = 'wind_dbs-3_'+date+'.nc'
pathIn             =  '/data/TR32/D2/data/wind_lidar/data/nc/2013/05/02/'
windLidarScan_data = Dataset(pathIn+filename, mode='r')
time_windScan      = windLidarScan_data.variables['time'][:].copy()
T_unix             = (time_windScan-2440587.5)*86400
T_units            = 'seconds since 1970-01-01 00:00:00'
datetime_windScan  = nc4.num2date(T_unix[:],T_units)
height_windScan    = windLidarScan_data.variables['height'][:].copy()
speed_windScan     = windLidarScan_data.variables['speed'][:].copy()
dir_windScan       = windLidarScan_data.variables['dir'][:].copy()
windVec_windScan   = windLidarScan_data.variables['wind_vec'][:].copy()
u_windScan         = windVec_windScan[0,:,:].T
v_windScan         = windVec_windScan[1,:,:].T


iconLemData        = Dataset('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130502.nc', mode='r')


time_ICON          = iconLemData.groups['Temp_data'].variables['datetime_ICON'][:].copy()
datetime_ICON      = nc4.num2date(iconLemData.groups['Temp_data'].variables['datetime_ICON'][:], \
                                  iconLemData.groups['Temp_data'].variables['datetime_ICON'].units) 
height_ICON        = iconLemData.groups['Temp_data'].variables['height2'][:].copy()
cloudMask_ICON     = iconLemData.groups['Temp_data'].variables['cloudMask'][:].copy()

ICON_DF            = pd.DataFrame(cloudMask_ICON, index=datetime_ICON, columns=height_ICON)     
ICON_DF_T          = pd.DataFrame(cloudMask_ICON.transpose(), index=height_ICON, columns=datetime_ICON)
U_iconRes          = f_resampling_twoD_Field(u_windScan, datetime_windScan, \
                                              height_windScan, ICON_DF, ICON_DF_T)
V_iconRes          = f_resampling_twoD_Field(v_windScan, datetime_windScan, \
                                              height_windScan, ICON_DF, ICON_DF_T)

#%%
w_plot = np.ma.array(u_windScan, mask=np.isnan(u_windScan))
fig, ax = plt.subplots(figsize=(12,6))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax1 = ax.pcolormesh(datetime_windScan, height_windScan, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
ax.set_ylim(400.,3000.)                                               # limits of the y-axe
ax.set_xlim()                                                        # limits of the x-axes
ax.set_title("u wind - JOYCE", fontsize=16)
ax.set_xlabel("time [hh:mm]", fontsize=16)
ax.set_ylabel("height [m]", fontsize=16)
plt.legend(loc='upper left')
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.ticks=([0,1,2,3])
cbar.set_label(label="horizontal wind [m/s]",size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
fig.tight_layout()
plt.savefig(pathDebugFig+'uwindScan_originalRes_'+date+'.png', format='png')



w_plot = np.ma.array(U_iconRes.values.T, mask=np.isnan(U_iconRes.values.T))
fig, ax = plt.subplots(figsize=(12,6))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
ax.set_ylim(400.,2300.)                                               # limits of the y-axe
ax.set_xlim()                                                        # limits of the x-axes
ax.set_title("u wind - JOYCE", fontsize=16)
ax.set_xlabel("time [$hh:mm$]", fontsize=16)
ax.set_ylabel("height [m]", fontsize=16)
plt.legend(loc='upper left')
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.ticks=([0,1,2,3])
cbar.set_label(label="u wind [$ms^{-1}$]",size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
fig.tight_layout()
plt.savefig(pathDebugFig+'uwindScan_RescaledICONRes_'+date+'.png', format='png')


w_plot = np.ma.array(v_windScan, mask=np.isnan(v_windScan))
fig, ax = plt.subplots(figsize=(12,6))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax1 = ax.pcolormesh(datetime_windScan, height_windScan, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
ax.set_ylim(400.,3000.)                                               # limits of the y-axe
ax.set_xlim()                                                        # limits of the x-axes
ax.set_title("v wind - JOYCE", fontsize=16)
ax.set_xlabel("time [hh:mm]", fontsize=16)
ax.set_ylabel("height [m]", fontsize=16)
plt.legend(loc='upper left')
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.ticks=([0,1,2,3])
cbar.set_label(label="v wind [$ms^{-1}$]",size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
fig.tight_layout()
plt.savefig(pathDebugFig+'vwindScan_originalRes_'+date+'.png', format='png')



w_plot = np.ma.array(V_iconRes.values.T, mask=np.isnan(V_iconRes.values.T))
fig, ax = plt.subplots(figsize=(12,6))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
ax.set_ylim(400.,2300.)                                               # limits of the y-axe
ax.set_xlim()                                                        # limits of the x-axes
ax.set_title("v wind - JOYCE", fontsize=16)
ax.set_xlabel("time [$hh:mm$]", fontsize=16)
ax.set_ylabel("height [m]", fontsize=16)
plt.legend(loc='upper left')
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.ticks=([0,1,2,3])
cbar.set_label(label="v wind [$ms^{-1}$]",size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
fig.tight_layout()
plt.savefig(pathDebugFig+'vwindScan_RescaledICONres_'+date+'.png', format='png')