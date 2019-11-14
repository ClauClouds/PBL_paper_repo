#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:32:34 2019
codice prova per leggere e plottare le surface longwave and shortwave fluxes del modello e dei dati da Pyrgeometer
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
from myFunctions import f_resamplingfield

pathFig        = '/work/cacquist/HDCP2_S2/statistics/debug/patch003/'
ICONfile       = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130502.nc'
date           = '20130502'
yy             = date[0:4]
mm             = date[4:6]
dd             = date[6:8]
iconLemData    = Dataset(ICONfile, mode='r')
timeiconLem    = iconLemData.groups['Temp_data'].variables['datetime_ICON'][:].copy()
LWiconLem      = iconLemData.groups['Temp_data'].variables['LWSurfFlux'][:]
SWiconLem      = iconLemData.groups['Temp_data'].variables['SWSurfFlux'][:]
time_ICON      = iconLemData.groups['Temp_data'].variables['datetime_ICON'][:].copy()
datetime_ICON  = nc4.num2date(iconLemData.groups['Temp_data'].variables['datetime_ICON'][:], \
                                  iconLemData.groups['Temp_data'].variables['datetime_ICON'].units) 
height_ICON    = iconLemData.groups['Temp_data'].variables['height2'][:].copy()
cloudMask_ICON = iconLemData.groups['Temp_data'].variables['cloudMask'][:].copy()

meteogram       = '/data/inscape/icon/experiments/juelich/meteo-4nest/METEOGRAM_patch003_20130502_joyce.nc'
meteogramData   = Dataset(meteogram, mode='r')
TempSurficonLem = meteogramData.variables['T_S'][:]
sigma           = 5.67*10**(-8) # W/m2K4
UpLWfluxiconLem = sigma*TempSurficonLem**4



fileLWObs = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/2013/sups_joy_pyrg00_l1_rlds_v00_20130502000000.nc'
LWobsData = Dataset(fileLWObs, mode='r')
LWobs     = LWobsData.variables['rlds'][:]
LWobsErr  = LWobsData.variables['rlds_error'][:]
timeLWObs = nc4.num2date(LWobsData.variables['time'][:], LWobsData.variables['time'].units) 

fileSWObs = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/2013/sups_joy_pyr00_l1_rsd_v00_20130502000000.nc'
SWobsData = Dataset(fileSWObs, mode='r')
SWobs     = SWobsData.variables['rsd'][:]
SWobsErr  = SWobsData.variables['rsd_error'][:]
timeSWObs = nc4.num2date(SWobsData.variables['time'][:], SWobsData.variables['time'].units) 


# resampling obvs to icon resolution
ICON_DF   = pd.DataFrame(cloudMask_ICON, index=datetime_ICON, columns=height_ICON)     
LWobs_res = f_resamplingfield(LWobs, timeLWObs, ICON_DF)
SWobs_res = f_resamplingfield(SWobs, timeSWObs, ICON_DF)
LWobsErr_res = f_resamplingfield(LWobsErr, timeLWObs, ICON_DF)
SWobsErr_res = f_resamplingfield(SWobsErr, timeSWObs, ICON_DF)


#%%
#creating new algorithm to detect cloud base and cloud top for liquid clouds
# idea: 
# - filter all ice clouds above 4000 mt ( they are not PBL clouds)
# - loop in height for every time: check 





#%% 

# calculating half hourly mean values of the longwave and shortwave surface fluxes for obs and model 
SW_DF_obs        = pd.Series(SWobs_res.values[:,0], index=datetime_ICON)
SW_obs_30min     = SW_DF_obs.resample('30min').mean() 
LW_DF_obs        = pd.Series(LWobs_res.values[:,0], index=datetime_ICON)
LW_obs_30min     = LW_DF_obs.resample('30min').mean() 

LW_DF_obs_Err    = pd.Series(LWobsErr_res.values[:,0], index=datetime_ICON)
LW_obs_Err_30min = LW_DF_obs_Err.resample('30min').mean() 
SW_DF_obs_Err    = pd.Series(SWobsErr_res.values[:,0], index=datetime_ICON)
SW_obs_Err_30min = SW_DF_obs_Err.resample('30min').mean() 


LW_mod           = LWiconLem+UpLWfluxiconLem
SW_DF_mod        = pd.Series(SWiconLem, index=datetime_ICON)
SW_mod_30min     = SW_DF_mod.resample('30min').mean() 
LW_DF_mod        = pd.Series(LW_mod, index=datetime_ICON)
LW_mod_30min     = LW_DF_mod.resample('30min').mean() 
datetime_30m     = [datetime.datetime(int(yy),int(mm),int(dd),0,0,0) + \
                    datetime.timedelta(minutes=60*x) for x in range(0, 49)]

data = np.arange(-1000, 2000)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_30m))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(275., 400.)
plt.ylim(275., 400.)
plt.xlabel('Longwave downward flux obs [W/m^2]', fontsize=16)
plt.ylabel('Longwave downward flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_30m)) 

LWF_iconlemPlot = LW_mod_30min.values
LWF_obsPlot     = LW_obs_30min.values
sizeDots         = LW_obs_Err_30min.values
cax = ax.scatter(LWF_obsPlot[:], LWF_iconlemPlot[:], c=colors, cmap=cmap, \
                 s=20*sizeDots)
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_30m[0])[11:16],\
                         str(datetime_30m[8])[11:16],\
                         str(datetime_30m[16])[11:16],\
                         str(datetime_30m[24])[11:16],\
                         str(datetime_30m[32])[11:16],\
                         str(datetime_30m[40])[11:16],\
                         str(datetime_30m[47])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'LW_scatterplot_obs_mod_allDays.png', format='png')



data = np.arange(-1000, 2000)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_30m))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(0., 1000.)
plt.ylim(0., 1000.)
plt.xlabel('Shortwave downward flux obs [W/m^2]', fontsize=16)
plt.ylabel('Shortwave downward flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_30m)) 

SWF_iconlemPlot  = SW_mod_30min.values
SWF_obsPlot      = SW_obs_30min.values
sizeDots         = SW_obs_Err_30min.values
cax = ax.scatter(SWF_obsPlot[:], SWF_iconlemPlot[:], c=colors, cmap=cmap, \
                 s=20*sizeDots)
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_30m[0])[11:16],\
                         str(datetime_30m[8])[11:16],\
                         str(datetime_30m[16])[11:16],\
                         str(datetime_30m[24])[11:16],\
                         str(datetime_30m[32])[11:16],\
                         str(datetime_30m[40])[11:16],\
                         str(datetime_30m[47])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'SW_scatterplot_obs_mod_allDays.png', format='png')

#%%
# plotting short wave downward fluxes 
fig, ax = plt.subplots(2, 1, figsize=(6,10))

ax2 = plt.subplot(2,1,1) 
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False)  
ax2.get_xaxis().tick_bottom()  
ax2.get_yaxis().tick_left()  
#ax2.set_ylim(-1.5, 1.5)                                               # limits of the y-axe
#ax2.set_xlim(0, 0.0006)
#ax2.set_title("ICON forward modeled skewness vs Qc ", fontsize=16)
ax2.set_xlabel("time [hh:mm]", fontsize=14)
ax2.set_ylabel("Longwave downward fluxes [W/m^2]", fontsize=14)
plt.plot(datetime_ICON, LWobs_res.values, color='black', label='obs')
y1        = LWobs_res.values-LWobsErr_res.values
y2        = LWobs_res.values+LWobsErr_res.values
plt.fill_between(datetime_ICON, y1[:,0], y2[:,0], facecolor='black', alpha=0.2)

plt.plot(datetime_ICON, LWiconLem+UpLWfluxiconLem, color='red', label='icon lem')
plt.legend(frameon=False)
#plt.hlines(0., 10**(-7), 0.003, color='black', linestyles=':')
#plt.vlines(0.0005, -1.5, 3.0, color='black', linestyles='dashed ')
myFmt = mdates.DateFormatter('%H:%M')
ax2.xaxis.set_major_formatter(myFmt)
#xticks(datetime_ICON)
ax2.set_xticks([datetime.datetime(2013,5,2,0), \
                datetime.datetime(2013,5,2,6), \
               datetime.datetime(2013,5,2,12),\
               datetime.datetime(2013,5,2,18),\
               datetime.datetime(2013,5,2,23)])
ax2.tick_params(labelsize=14)
plt.tight_layout()
ax3 =plt.subplot(2,1,2) 
ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
ax3.get_xaxis().tick_bottom()  
ax3.get_yaxis().tick_left()  
#ax3.set_ylim(-1.5, 1.5)                                               # limits of the y-axe
#ax3.set_xlim(0, 0.0006)
#ax3.set_title("ICON forward modeled skewness vs Qc+Qr ", fontsize=16)
ax3.set_xlabel("time [hh:mm]", fontsize=14)
ax3.set_ylabel("Shortwave downward fluxes [W/m^2]", fontsize=14)
plt.plot(datetime_ICON, SWobs_res, color='black', label='obs')
y1        = SWobs_res.values-SWobsErr_res.values
y2        = SWobs_res.values+SWobsErr_res.values
plt.fill_between(datetime_ICON, y1[:,0], y2[:,0], facecolor='black', alpha=0.2)


plt.plot(datetime_ICON, SWiconLem, color='red', label='icon lem')
plt.legend(frameon=False)
myFmt = mdates.DateFormatter('%H:%M')
ax3.xaxis.set_major_formatter(myFmt)
#xticks(datetime_ICON)
ax3.set_xticks([datetime.datetime(2013,5,2,0), \
                datetime.datetime(2013,5,2,6), \
               datetime.datetime(2013,5,2,12),\
               datetime.datetime(2013,5,2,18),\
               datetime.datetime(2013,5,2,23)])
ax3.tick_params(labelsize=14)
#plt.hlines(0., 10**(-7), 0.003, color='black', linestyles=':')
#plt.vlines(0.0005, -1.5, 3.0, color='black', linestyles='dashed')
#ax3.set_xticks([0.0005, 0.0015, 0.0025])
plt.savefig(pathFig+'LWF_timeSeries_obs_mod.png', format='png')


#%%

data = np.arange(2000)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_ICON))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(250., 400.)
plt.ylim(250., 400.)
plt.xlabel('LWF obs [W/m^2]', fontsize=16)
plt.ylabel('LWF icon lem [W/m^2]', fontsize=16)
cmap = plt.cm.get_cmap('bwr', len(datetime_ICON)) 
cax = ax.scatter(LWobs_res.values, LWiconLem+UpLWfluxiconLem, c=colors, cmap=cmap, \
                 s=100)
plt.grid(b=True, which='major', color='#666666', linestyle=':')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 1600, 3200, 4800, 6400,8000, 9599])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_ICON[0])[11:16],\
                         str(datetime_ICON[1600])[11:16],\
                         str(datetime_ICON[3200])[11:16],\
                         str(datetime_ICON[4800])[11:16],\
                         str(datetime_ICON[6400])[11:16],\
                         str(datetime_ICON[8000])[11:16],\
                         str(datetime_ICON[9599])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'LWF_scatterplot_obs_mod.png', format='png')




data = np.arange(2000)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_ICON))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(0., 1200.)
plt.ylim(0., 1200.)
plt.xlabel('SWF obs [W/m^2]', fontsize=16)
plt.ylabel('SWF icon lem [W/m^2]', fontsize=16)
cmap = plt.cm.get_cmap('bwr', len(datetime_ICON)) 
cax = ax.scatter(SWobs_res.values, SWiconLem, c=colors, cmap=cmap, \
                 s=100)
plt.grid(b=True, which='major', color='#666666', linestyle=':')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 1600, 3200, 4800, 6400,8000, 9599])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_ICON[0])[11:16],\
                         str(datetime_ICON[1600])[11:16],\
                         str(datetime_ICON[3200])[11:16],\
                         str(datetime_ICON[4800])[11:16],\
                         str(datetime_ICON[6400])[11:16],\
                         str(datetime_ICON[8000])[11:16],\
                         str(datetime_ICON[9599])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'SWF_scatterplot_obs_mod.png', format='png')

#%%