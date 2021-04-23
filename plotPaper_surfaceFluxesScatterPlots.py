#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:16:28 2020

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
date_arr = []
LHSF_mod                        = []
SHSF_mod                        = []
LHSF_obs                        = []
SHSF_obs                        = []
LHSF_err_obs                    = []
SHSF_err_obs                    = []
datetime_30m                    = []


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
    
    LHSF_mod.append(new_dict[10]['LHF_iconlem'])
    SHSF_mod.append(new_dict[10]['SHF_iconlem'])
    LHSF_obs.append(new_dict[10]['LHF_obs'])
    SHSF_obs.append(new_dict[10]['SHF_obs'])
    LHSF_err_obs.append(new_dict[10]['LHF_Err_obs'])
    SHSF_err_obs.append(new_dict[10]['SHF_Err_obs'])     
    datetime_30m.append(new_dict[10]['datetime_30m'])
    
    # =============================================================================
# plotting scatter plots of surface latent and sensible heat flux for every half hour 
# =============================================================================
#%%
##### Latent heat surface flux
datetime_fluxes = datetime_30m[0][:-1]
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
colors = np.arange(0,len(datetime_fluxes))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(0., 500.)
plt.ylim(0., 500.)
plt.xlabel('Latent heat surface flux obs [W/m^2]', fontsize=16)
plt.ylabel('Latent heat surface flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')
cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 

for indFile in range(Nfiles):
    LHSF_iconlemPlot = -LHSF_mod[indFile]
    LHSF_obsPlot = LHSF_obs[indFile]
    sizeDots = LHSF_err_obs[indFile]
    cax = ax.scatter(LHSF_obsPlot[:], LHSF_iconlemPlot[:-1], c=colors, cmap=cmap, \
                 s=10*sizeDots, edgecolors='black')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
                         str(datetime_fluxes[8])[11:16],\
                         str(datetime_fluxes[16])[11:16],\
                         str(datetime_fluxes[24])[11:16],\
                         str(datetime_fluxes[32])[11:16],\
                         str(datetime_fluxes[40])[11:16],\
                         str(datetime_fluxes[47])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'LHSF_scatterplot_obs_mod_allDays.png', format='png')
print('Latent heat plotted')

#%%
###### Sensible heat surface flux
datetime_fluxes = datetime_30m[0][:-1]
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
colors = np.arange(0,len(datetime_fluxes))
plt.plot(data, data, color='black', linestyle=':')
plt.xlim(-100., 450.)
plt.ylim(-100., 450.)
plt.xlabel('Sensible heat surface flux obs [W/m^2]', fontsize=16)
plt.ylabel('Sensible heat surface flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')
cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 

for indFile in range(Nfiles):
    print(indFile)
    SHSF_iconlemPlot = -SHSF_mod[indFile]
    SHSF_obsPlot = SHSF_obs[indFile]
    sizeDots = SHSF_err_obs[indFile]
    cax = ax.scatter(SHSF_obsPlot[:], SHSF_iconlemPlot[:-1], c=colors, cmap=cmap, \
                 s=10*sizeDots, edgecolors='black')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
                         str(datetime_fluxes[8])[11:16],\
                         str(datetime_fluxes[16])[11:16],\
                         str(datetime_fluxes[24])[11:16],\
                         str(datetime_fluxes[32])[11:16],\
                         str(datetime_fluxes[40])[11:16],\
                         str(datetime_fluxes[47])[11:16]], fontsize=14) 
plt.tight_layout()
plt.savefig(pathFig+'SHSF_scatterplot_obs_mod_allDays.png', format='png')

#%%
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8,12))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax[0] = plt.subplot(2,1,1)  
ax[0].spines["top"].set_visible(False)  
ax[0].spines["right"].set_visible(False)  
ax[0].get_xaxis().tick_bottom()  
ax[0].get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_fluxes))
ax[0].plot(data, data, color='black', linestyle=':')
ax[0].set_xlim(-100., 450.)
ax[0].set_ylim(-100., 450.)
ax[0].set_xlabel('Sensible heat surface flux obs [W/m^2]', fontsize=16)
ax[0].set_ylabel('Sensible heat surface flux icon lem [W/m^2]', fontsize=16)
ax[0].grid(b=True, which='major', color='#666666', linestyle=':')
cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
SHSFall_obs = []
SHSFall_mod = []
LHSFall_obs = []
LHSFall_mod = []

for indFile in range(Nfiles):
    print(indFile)
    SHSF_iconlemPlot = -SHSF_mod[indFile]
    SHSF_obsPlot = SHSF_obs[indFile]
    sizeDots = SHSF_err_obs[indFile]
    SHSFall_obs.append(SHSF_obs[indFile])
    SHSFall_mod.append(-SHSF_mod[indFile][:-1])
    cax = ax[0].scatter(SHSF_obsPlot[:], SHSF_iconlemPlot[:-1], c=colors, cmap=cmap, \
                 s=10*sizeDots, edgecolors='black')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
                         str(datetime_fluxes[8])[11:16],\
                         str(datetime_fluxes[16])[11:16],\
                         str(datetime_fluxes[24])[11:16],\
                         str(datetime_fluxes[32])[11:16],\
                         str(datetime_fluxes[40])[11:16],\
                         str(datetime_fluxes[47])[11:16]], fontsize=14) 
x = np.hstack(SHSFall_obs)
y = np.hstack(SHSFall_mod)
xfit = []
yfit = []
for ind in range(len(x)):
    if ~(np.isnan(x[ind])) & ~(np.isnan(y[ind])):
        xfit.append(x[ind])
        yfit.append(y[ind])

m, b= np.polyfit(xfit, yfit, 1)
dataFit = np.arange(0., 500.)
plt.plot(dataFit, b + m*dataFit, color='red')

ax[1] = plt.subplot(2,1,2)  
ax[1].spines["top"].set_visible(False)  
ax[1].spines["right"].set_visible(False)  
ax[1].get_xaxis().tick_bottom()  
ax[1].get_yaxis().tick_left() 
colors = np.arange(0,len(datetime_fluxes))
ax[1].plot(data, data, color='black', linestyle=':')
ax[1].set_xlim(0., 500.)
ax[1].set_ylim(0., 500.)
ax[1].set_xlabel('Latent heat surface flux obs $[Wm^{-2}]$', fontsize=16)
ax[1].set_ylabel('Latent heat surface flux model $[Wm^{-2}]$', fontsize=16)
ax[1].grid(b=True, which='major', color='#666666', linestyle=':')
cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 

for indFile in range(Nfiles):
    LHSF_iconlemPlot = -LHSF_mod[indFile]
    LHSF_obsPlot = LHSF_obs[indFile]
    sizeDots = LHSF_err_obs[indFile]
    LHSFall_obs.append(LHSF_obs[indFile])
    LHSFall_mod.append(-LHSF_mod[indFile][:-1])
    cax = ax[1].scatter(LHSF_obsPlot[:], LHSF_iconlemPlot[:-1], c=colors, cmap=cmap, \
                 s=10*sizeDots, edgecolors='black')
cbar = fig.colorbar(cax, \
                    cmap=cmap, \
                    ticks= [0, 8, 16, 24, 32, 40, 47])
cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
                         str(datetime_fluxes[8])[11:16],\
                         str(datetime_fluxes[16])[11:16],\
                         str(datetime_fluxes[24])[11:16],\
                         str(datetime_fluxes[32])[11:16],\
                         str(datetime_fluxes[40])[11:16],\
                         str(datetime_fluxes[47])[11:16]], fontsize=14) 
x = np.hstack(LHSFall_obs)
y = np.hstack(LHSFall_mod)
xfit = []
yfit = []
for ind in range(len(x)):
    if ~(np.isnan(x[ind])) & ~(np.isnan(y[ind])):
        xfit.append(x[ind])
        yfit.append(y[ind])
mm, bb= np.polyfit(xfit, yfit, 1)
plt.plot(dataFit, bb + mm*dataFit, color='red')

plt.tight_layout()
plt.savefig(pathFig+'LHSF_SHSF_scatterplot_obs_mod_allDays.png', format='png')


#%%
##### Evaporative fraction 
#datetime_fluxes = datetime_30m[0][:-1]
#data = np.arange(-1000, 2000)
#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
#matplotlib.rcParams['savefig.dpi'] = 100
#plt.gcf().subplots_adjust(bottom=0.15)
#fig.tight_layout()
#ax = plt.subplot(1,1,1)  
#ax.spines["top"].set_visible(False)  
#ax.spines["right"].set_visible(False)  
#ax.get_xaxis().tick_bottom()  
#ax.get_yaxis().tick_left() 
#colors = np.arange(0,len(datetime_fluxes))
#plt.plot(data, data, color='black', linestyle=':')
#plt.xlim(-5., 5.)
#plt.ylim(-50., 50.)
#plt.xlabel('Evaporative fraction obs [W/m^2]', fontsize=16)
#plt.ylabel('Evaporative fraction icon lem [W/m^2]', fontsize=16)
#plt.grid(b=True, which='major', color='#666666', linestyle=':')
#
#
## calculating evaporative fraction as ratio of latent/(sens+latent heat)
#EvapFraction_iconlem = []
#EvapFraction_obs = []
#
#for indFile in range(Nfiles):
#    print(indFile)
#    SHSF_iconlemPlot = -SHSF_mod[indFile]
#    SHSF_obsPlot = SHSF_obs[indFile]
#    LHSF_iconlemPlot = LHSF_mod[indFile]
#    LHSF_obsPlot = -LHSF_obs[indFile]
#    sizeDots = 5#SHSF_err_obs[indFile]
#    EvapFractionDay_iconlem = np.zeros((len(SHSF_iconlemPlot)))
#    EvapFractionDay_iconlem.fill(np.nan)
#    EvapFractionDay_obs = np.zeros((len(SHSF_obsPlot)))
#    EvapFractionDay_obs.fill(np.nan)
#    for ind in range(len(SHSF_iconlemPlot)):
#        EvapFractionDay_iconlem[ind] = LHSF_iconlemPlot[ind]/(LHSF_iconlemPlot[ind]+SHSF_iconlemPlot[ind])
#    for ind in range(len(SHSF_obsPlot)):        
#        EvapFractionDay_obs[ind] = LHSF_obsPlot[ind]/(LHSF_obsPlot[ind]+SHSF_obsPlot[ind])
#    EvapFraction_iconlem.append(EvapFractionDay_iconlem)
#    EvapFraction_obs.append(EvapFractionDay_obs)
#    
#cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
#for indFile in range(Nfiles):
#    EvapFractionPlot_iconlem = EvapFraction_iconlem[indFile]
#    EvapFractionPlot_obs =EvapFraction_obs[indFile]
#    
#    cax = ax.scatter(EvapFractionPlot_obs[:], EvapFractionPlot_iconlem[:-1], c=colors, cmap=cmap, \
#                 s=10*sizeDots)
#cbar = fig.colorbar(cax, \
#                    cmap=cmap, \
#                    ticks= [0, 8, 16, 24, 32, 40, 47])
#cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
#cbar.ax.tick_params(labelsize=14)
#cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
#                         str(datetime_fluxes[8])[11:16],\
#                         str(datetime_fluxes[16])[11:16],\
#                         str(datetime_fluxes[24])[11:16],\
#                         str(datetime_fluxes[32])[11:16],\
#                         str(datetime_fluxes[40])[11:16],\
#                         str(datetime_fluxes[47])[11:16]], fontsize=14) 
#plt.tight_layout()
#plt.savefig(pathFig+'EvapFraction_scatterplot_obs_mod_allDays.png', format='png')