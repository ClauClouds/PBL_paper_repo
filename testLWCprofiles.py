#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 17:22:40 2019

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
import datetime as dt
import random
import datetime
import matplotlib.dates as mdates
import os  
import atmos
import matplotlib as mpl
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
path = '/data/hatpro/jue/cloudnet/juelich/products/lwc-scaled-adiabatic/2013/'
file = '20130502_juelich_lwc-scaled-adiabatic.nc'


data = Dataset(path+file, mode='r')
LWC = data.variables['lwc'][:]
LWC_err = data.variables['lwc_error'][:]
LWC_retrieval_status = data.variables['lwc_retrieval_status'][:]
time = nc4.num2date(data.variables['time'][:], \
                     data.variables['time'].units)
height = data.variables['height'][:]
LWC[LWC==0]= np.nan
LWC_plot = np.ma.array(LWC, mask=np.isnan(LWC))
#%%
fig, ax = plt.subplots(figsize=(12,6))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
LWC[LWC==0]= np.nan
LWC_plot = np.ma.array(LWC, mask=np.isnan(LWC))

cax1 = ax.pcolormesh(time, height, LWC_plot.transpose(), \
                     vmin=np.nanmin(LWC_plot), vmax=np.nanmax(LWC_plot), cmap='plasma')
ax.set_ylim(400.,3000.)                                               # limits of the y-axe
ax.set_xlim()                                                        # limits of the x-axes
ax.set_title("$LWC - JOYCE$", fontsize=16)
ax.set_xlabel("$time \ [hh:mm]$", fontsize=16)
ax.set_ylabel("$height \ [m]$", fontsize=16)
#plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#plt.plot(time_ICON, CB_array, color='black',label='cloud base')
plt.legend(loc='upper left')
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.ticks=([0,1,2,3])
#cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
cbar.set_label(label='LWC $[kg m^{-3}]$',size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
fig.tight_layout()
#plt.savefig(pathDebugFig+'Hwind_iconlem_joyce_'+date+'.png', format='png')

#%%
# plot hourly mean LWC profiles
NprofilesOut = 24
LWCmean_obs  = np.zeros((NprofilesOut, len(height)))
datetime_60m = [datetime.datetime(int(2013),int(5),int(2),0,0,0) + \
                        datetime.timedelta(minutes=60*x) for x in range(0, 25)]
for indHour in range(len(datetime_60m)-1):
    mask_t = np.where(((time > datetime_60m[indHour]) * (time <= datetime_60m[indHour+1])))[0]
    LWCmean_obs[indHour, :] = np.nanmean(LWC[mask_t,:], axis=0)
    if indHour == 8:
        print(np.shape(np.nanmean(LWC[mask_t,:], axis=0)))
        print(np.nanmin(LWC[mask_t,:]))
        print(np.nanmax(LWC[mask_t,:]))
        print(np.nanmean(LWC[mask_t,:], axis=0))
        
    
#%%
# plotting hourly mean LWC profiles 
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
plt.xlim(0., 0.01)
plt.ylim(0., 3000.)
plt.xlabel('LWC $[kg m^{-3}]$', fontsize=16)
plt.ylabel('height [m]', fontsize=16)
cmap = plt.cm.get_cmap('bwr', 19) 
for ind in range(len(datetime_60m)-1):
    LWC_plot = LWCmean_obs[ind,:] + ind*np.repeat(0.0005, len(height))
    plt.plot(LWC_plot, height)
plt.legend()
plt.tight_layout()
#plt.savefig(pathFig+'pblh_scatterplot_obs_mod.png', format='png')
    