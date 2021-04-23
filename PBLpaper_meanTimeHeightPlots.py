

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur Nov 21 14:13 2019
@ author: cacquist
@ contact: cacquist@meteo.uni-koeln.de
@  date : 21 November 2019, modified on 19 Dec 2019
@  goal : elaborate mean time height plots of
- u
- v
- w
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

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
# setting parameters for calculating averaging and domain size of the model:
NprofilesOut  = 24  # hourly means
timeIncrement = 60  # hourly means
patch         = 'patch003'

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/old/'
#pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathFig = pathObs
#filenameObs = 'dataset_PBLcloudPaper_ModObs_20130502.p'
#filenameMod = 'icon_lem_derivedproperties20130502.nc'

fileListObs               = sorted(glob.glob(pathObs+'*.p'))
fileListMod               = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                    = len(fileListMod)    
date_arr  = []
skewW_mod = []
skewW_obs = []
w_mod     = []
w_obs     = []
uWind_obs = []
uWind_mod = []
vWind_obs = []
vWind_mod = []
timeArr   = []
# loop on the ensemble of days (statistics)
for indFile in range(Nfiles):
    
    filenameMod   = fileListMod[indFile]
    filenameObs   = fileListObs[indFile]
    date = fileListMod[indFile][91:99]
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])    
    date_arr.append(date)
    print('processing date '+date)
    # reading time and height from ncdf file (grid of ICON LEM 
    #( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata        = Dataset(filenameMod, mode='r')
    time          = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                     ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height        = ncdata.groups['Temp_data'].variables['height2']
    timeArr.append(time)
    # opening the file containing all the data
    infile        = open(filenameObs,'rb')
    new_dict      = pickle.load(infile, encoding='latin1')
    
    # reading skewness of vertical velocity from model and obs
    skewW_m = ncdata.groups['Temp_data'].variables['skewnessW']
    skewW_o = new_dict[3]['skewnessW']    

    if (skewW_m.shape != (9600,150)):
         fullSkew_m = np.zeros((9600, 150))
         fullSkew_m.fill(np.nan)
         fullSkew_m[0:len(time), :] = skewW_m[:,:]
         fullSkew_o = np.zeros((9600, 150))
         fullSkew_o.fill(np.nan)
         fullSkew_o[0:len(time), :] = skewW_o[:,:]
         skewW_mod.append(fullSkew_m)
         skewW_obs.append(fullSkew_o)
    else:
        skewW_mod.append(skewW_m)
        skewW_obs.append(skewW_o)
        
    # reading vertical velocity from model
    w_m = new_dict[9]['w_iconlem']
    w_o = new_dict[3]['verticalWind']


    
    # reading horizontal wind from model
    uWind_m = new_dict[9]['u_iconlem']
    uWind_o = new_dict[3]['zonalWind']

    if (w_m.shape[0] != 9600):
        if (w_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = w_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = w_o[:,1:]
             w_mod.append(full_m)
             w_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = w_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = w_o[:,:]
             w_mod.append(full_m)
             w_obs.append(full_o)
    else:
        if (w_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = w_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = w_o[:,1:]
             w_mod.append(full_m)
             w_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = uWind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = uWind_o[:,:]
             w_mod.append(full_m)
             w_obs.append(full_o)
             
             

    
    # reading horizontal wind from model
    uWind_m = new_dict[9]['u_iconlem']
    uWind_o = new_dict[3]['zonalWind']

    if (uWind_m.shape[0] != 9600):
        if (uWind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = uWind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = uWind_o[:,1:]
             uWind_mod.append(full_m)
             uWind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = uWind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = uWind_o[:,:]
             uWind_mod.append(full_m)
             uWind_obs.append(full_o)
    else:
        if (uWind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = uWind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = uWind_o[:,1:]
             uWind_mod.append(full_m)
             uWind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = uWind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = uWind_o[:,:]
             uWind_mod.append(full_m)
             uWind_obs.append(full_o)
             
             

    
    # reading horizontal wind from model
    vWind_m = new_dict[9]['v_iconlem']
    vWind_o = new_dict[3]['meridionalWind']

    if (vWind_m.shape[0] != 9600):
        if (vWind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = vWind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = vWind_o[:,1:]
             vWind_mod.append(full_m)
             vWind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = vWind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = vWind_o[:,:]
             vWind_mod.append(full_m)
             vWind_obs.append(full_o)
    else:
        if (vWind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = vWind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = vWind_o[:,1:]
             vWind_mod.append(full_m)
             vWind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = vWind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = vWind_o[:,:]
             vWind_mod.append(full_m)
             vWind_obs.append(full_o)
    
#%%

W_allDays_mod     = np.dstack(w_mod[:])
Uwind_allDays_mod = np.dstack(uWind_mod[:])
skewW_allDays_mod = np.dstack(skewW_mod[:])
Vwind_allDays_mod = np.dstack(vWind_mod[:])

W_allDays_obs     = np.dstack(w_obs[:])
Uwind_allDays_obs = np.dstack(uWind_obs[:])
skewW_allDays_obs = np.dstack(skewW_obs[:])
Vwind_allDays_obs = np.dstack(vWind_obs[:])

# calculating mean values 
meanW_mod = np.nanmean(W_allDays_mod, axis=2)
meanW_obs = np.nanmean(W_allDays_obs, axis=2)

meanUwind_mod = np.nanmean(Uwind_allDays_mod, axis=2)
meanUwind_obs = np.nanmean(Uwind_allDays_obs, axis=2)

meanVwind_mod = np.nanmean(Vwind_allDays_mod, axis=2)
meanVwind_obs = np.nanmean(Vwind_allDays_obs, axis=2)

meanSkewW_mod = np.nanmean(skewW_allDays_mod, axis=2)
meanSkewW_obs = np.nanmean(skewW_allDays_obs, axis=2)


#%%
timeStart= datetime.datetime(2013,4,24,6,0,0)
timeEnd= datetime.datetime(2013,4,24,20,0,0)
Hmax=2000.
fontSizeTitle = 16
fontSizeX     = 10
fontSizeY     = 10
# =============================================================================
#             print('plotting graphs for debugging in debugging mode')         
fig, ax = plt.subplots(4, 2, figsize=(14,10), sharex=True, sharey=True, constrained_layout=True)
ax1 = plt.subplot(4,2,1) 
meanSkewW_obs_plot =  np.ma.array(meanSkewW_obs, mask=np.isnan(meanSkewW_obs))
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False)  
ax1.get_xaxis().tick_bottom()  
ax1.get_yaxis().tick_left()
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax1.xaxis_date()
cax = ax1.pcolormesh(timeArr[0], height, meanSkewW_obs_plot.transpose(), vmin=-2., vmax=2., cmap='PiYG')
ax1.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax1.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax1.set_title("mean skewness W obs", fontsize=16)#
ax1.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax1.set_ylabel("height [m]", fontsize=fontSizeY)
plt.tight_layout()

#cbar = fig.colorbar(cax, orientation='vertical')
#cbar.set_label(label="skewness W",size=16)
#cbar.ax.tick_params(labelsize=14)
#cbar.aspect=80

ax2 = plt.subplot(4,2,2) 
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False)  
ax2.get_xaxis().tick_bottom()  
ax2.get_yaxis().tick_left()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis_date()
cax = ax2.pcolormesh(timeArr[0], height, meanSkewW_mod.transpose(), vmin=-2., vmax=2., cmap='PiYG')
ax2.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax2.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax2.set_title("mean skewness W ICON_LEM", fontsize=16)#
ax2.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax2.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label="$skewness_W$",size=16)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
#plt.savefig(pathDebugFig+'meridionalWind_iconlem_'+date+'.png', format='png')
plt.tight_layout()

ax3 = plt.subplot(4,2,3) 
meanW_obs_plot =  np.ma.array(meanW_obs[:,:-1], mask=np.isnan(meanW_obs[:,:-1]))

ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
ax3.get_xaxis().tick_bottom()  
ax3.get_yaxis().tick_left()
ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis_date()
cax = ax3.pcolormesh(timeArr[0], height, meanW_obs_plot.transpose(), vmin=-10., vmax=10., cmap='seismic')
ax3.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax3.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax3.set_title("mean vertical velocity W ICON_LEM", fontsize=16)#
ax3.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax3.set_ylabel("height [m]", fontsize=fontSizeY)
plt.tight_layout()


ax4 = plt.subplot(4,2,4) 
meanW_mod_plot =  np.ma.array(meanW_mod[:,:-1], mask=np.isnan(meanW_mod[:,:-1]))
ax4.spines["top"].set_visible(False)  
ax4.spines["right"].set_visible(False)  
ax4.get_xaxis().tick_bottom()  
ax4.get_yaxis().tick_left()
ax4.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax4.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax4.xaxis_date()
cax = ax4.pcolormesh(timeArr[0], height, meanW_mod_plot.transpose(), vmin=-10., vmax=10., cmap='seismic')
ax4.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax4.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax4.set_title("mean vertical velocity W ICON_LEM", fontsize=16)#
ax4.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax4.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label="$w \ [ms^{-1}$]",size=16)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
plt.tight_layout()


ax5 = plt.subplot(4,2,5) 
uWind_obs_plot =  np.ma.array(meanUwind_obs[:,:-1], mask=np.isnan(meanUwind_obs[:,:-1]))

ax5.spines["top"].set_visible(False)  
ax5.spines["right"].set_visible(False)  
ax5.get_xaxis().tick_bottom()  
ax5.get_yaxis().tick_left()
ax5.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax5.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax5.xaxis_date()
cax = ax5.pcolormesh(timeArr[0], height, uWind_obs_plot.transpose(), vmin=-5., vmax=5, cmap='YlGnBu')
ax5.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax5.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax5.set_title("mean zonal wind obs [$ms^{-1}$]", fontsize=16)#
ax5.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax5.set_ylabel("height [m]", fontsize=fontSizeY)
plt.tight_layout()


ax6 = plt.subplot(4,2,6) 
uWind_mod_plot =  np.ma.array(meanUwind_mod[:,:-1], mask=np.isnan(meanUwind_mod[:,:-1]))
ax6.spines["top"].set_visible(False)  
ax6.spines["right"].set_visible(False)  
ax6.get_xaxis().tick_bottom()  
ax6.get_yaxis().tick_left()
ax6.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax6.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax6.xaxis_date()
cax = ax6.pcolormesh(timeArr[0], height, uWind_mod_plot.transpose(), vmin=-5., vmax=5., cmap='YlGnBu')
ax6.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax6.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax6.set_title("mean zonal wind ICON_LEM", fontsize=16)#
ax6.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax6.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label=" $u \ [ms^{-1}$]",size=16)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
plt.tight_layout()





ax7 = plt.subplot(4,2,7) 
vWind_obs_plot =  np.ma.array(meanVwind_obs[:,:-1], mask=np.isnan(meanVwind_obs[:,:-1]))

ax7.spines["top"].set_visible(False)  
ax7.spines["right"].set_visible(False)  
ax7.get_xaxis().tick_bottom()  
ax7.get_yaxis().tick_left()
ax7.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax7.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax7.xaxis_date()
cax = ax7.pcolormesh(timeArr[0], height, vWind_obs_plot.transpose(), vmin=-5., vmax=5., cmap='YlGnBu')
ax7.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax7.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax5.set_title("mean meridional wind obs [$ms^{-1}$]", fontsize=16)#
ax7.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax7.set_ylabel("height [m]", fontsize=fontSizeY)
plt.tight_layout()


ax8 = plt.subplot(4,2,8) 
vWind_mod_plot =  np.ma.array(meanVwind_mod[:,:-1], mask=np.isnan(meanVwind_mod[:,:-1]))
ax8.spines["top"].set_visible(False)  
ax8.spines["right"].set_visible(False)  
ax8.get_xaxis().tick_bottom()  
ax8.get_yaxis().tick_left()
ax8.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax8.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax8.xaxis_date()
cax = ax8.pcolormesh(timeArr[0], height, vWind_mod_plot.transpose(), vmin=-5., vmax=5., cmap='YlGnBu')
ax8.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax8.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax8.set_title("mean wind speed ICON_LEM", fontsize=16)#
ax8.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax8.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label="$v \ [ms^{-1}$]", size=16)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
plt.tight_layout()

#fig.subplots_adjust(hspace=0.4, wspace=0.4)


plt.savefig('/work/cacquist/HDCP2_S2/statistics/figs/patch003/stat_mean_field_w_skew_Hwind.png')
# ============================================================================
#%%




from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(14,3),
                  gridspec_kw={"width_ratios":[1,1, 0.05]}, )# constrained_layout = true
fig.subplots_adjust(wspace=0.3)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()

ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False)  
ax2.get_xaxis().tick_bottom()  
ax2.get_yaxis().tick_left()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis_date()
im  = ax.pcolormesh(timeArr[0], height, vWind_mod_plot.transpose(), vmin=-5., vmax=5., cmap='YlGnBu')
im2 = ax2.pcolormesh(timeArr[0], height, vWind_obs_plot.transpose(), vmin=-5., vmax=5., cmap='YlGnBu')
ax.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax5.set_title("mean meridional wind obs [$ms^{-1}$]", fontsize=16)#
ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax.set_ylabel("height [m]", fontsize=fontSizeY)
ax2.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax2.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax5.set_title("mean meridional wind obs [$ms^{-1}$]", fontsize=16)#
ax2.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax2.set_ylabel("height [m]", fontsize=fontSizeY)

ip = InsetPosition(ax2, [1.05,0,0.05,1]) 
cax.set_axes_locator(ip)

cbar = fig.colorbar(im, cax=cax, ax=[ax,ax2])
cbar.set_label(label="$v \ [ms^{-1}$]", size=16)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
fig.tight_layout()
plt.savefig('/work/cacquist/HDCP2_S2/statistics/figs/patch003/stat_mean_field_w_skew_Hwind.png')

































# =============================================================================
#         fig, ax = plt.subplots(figsize=(14,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis_date()
#         cax = ax.pcolormesh(datetime_ICON, height, zonalWind.transpose(), vmin=-8., vmax=20.)
#         ax.set_ylim(107.,2000.)    # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
#         ax.set_xlim()                                 # limits of the x-axes
#         ax.set_title("zonal Wind ICON", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         cbar = fig.colorbar(cax, orientation='vertical')
#         cbar.set_label(label="wind shear",size=16)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         plt.savefig(pathDebugFig+'zonalWind_iconlem_'+date+'.png', format='png')
#         
# =============================================================================