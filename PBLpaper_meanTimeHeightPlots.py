

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur Nov 21 14:13 2019
@ author: cacquist
@ contact: cacquist@meteo.uni-koeln.de
@  date : 21 November 2019
@  goal : elaborate mean time height plots of
- u
- v
- w
- RH 
general information on the content of the extracted pckle file:
# the files .p contain the following array of dictionaries/lists: 
#     outputArray   = [radiosondeList, tower_dict, dictCosmo, dictObsWindLidarMwr, Thermodyn_cosmo, \
#                   Thermodyn_iconlem, Thermodyn_obs, dynamics_iconlem, cloudDict_iconlem, cloudDict_obs]
# Below we provide the definitions of each of them:
# radiosondeList: 
# =============================================================================
# dict_day           = {
#                 'time':DatetimeRadiosonde,
#                 'P':P,
#                 'T':T,
#                 'Td': Td,
#                 'Td_surf': Td[0],
#                 'RH':RH,
#                 'z_lcl':z_lcl,
#                 'z_ccl':z_ccl, 
#                 'T_ccl':T_ground_CCL, 
#                 'PBLheight':PBLheight,
#                 'EISWood':EIS,
#                 'EIS2':EIS2,
#                 'LTS':LTS, 
#                 'theta_v':Theta_v,
#                 'surfaceTemperature':T_surf,
#                 }
# =============================================================================
# tower_dict:     
# =============================================================================
# dictOut={
#          'time':datetime_tower, 
#          'T':T, 
#          'P':P, 
#          'windSpeed':windSpeed,
#          'wDir':wDir,
#      #    'RH':RHArray,
#          'height':height,
#          'Tsurf':Tsurf,
#      #    'RHsurf':RHsurf
#      }
# =============================================================================
# dictCosmo
# =============================================================================
#     dictCosmo = {'pressure':P_cosmo_res.values.transpose(),
#                  'temperature':T_cosmo_res.values.transpose(),
#                  'absoluteHumidity':Q_cosmo_res.values.transpose(),
#                  'cloudnetCategorization':cloudnet_res.data.transpose()
#                  }
# 
# =============================================================================
# dictObsWindLidarMwr
# =============================================================================
#     dictObsWindLidarMwr = {
#             'verticalWind':w_obs_res.values.transpose(),
#             'horizontalWind':Hwind_obs_res.values.transpose(),
#             'skewnessW':skew_obs_res.values.transpose(),
#             'PBLclassObs':PBLclass_obs_res.values.transpose(),
#             'shear':shear_obs_res.values.transpose(),
#             'windDirection':wDir_obs_res.values.transpose(),
#             'absoluteHumidity':qProf_obs_res.values.transpose(),
#             'temperature':tProf_obs_res.values.transpose(),
#             'IWV_mwr':IWV_obs_res, 
#             'LWP_mwr':LWP_obs_res,
#             }
# =============================================================================
# Thermodyn_cosmo, Thermodyn_iconlem, Thermodyn_obs, 
# =============================================================================
# 
#     ThermodynPar={'mixingRatio':r, 
#                   'relativeHumidity':rh, 
#                   'virtualTemperature':tv,
#                   'cclHeight':result_ccl['z_ccl'],
#                   'cclTemperature':result_ccl['T_ccl'],
#                   'lclHeight':lclArray,
#                   'surfaceTemperature':TSurf, 
#                   'virtualPotentialTemperature':Theta_v,
#                   'time': time,
#                   'height':height,
#                   }
# =============================================================================
# dynamics_iconlem
# =============================================================================
#     DynPar={'varianceW':varW, 
#             'PBLHeight':PBLHeightArr, 
#             'windSpeed':windData_ICON['windSpeed'], 
#             'windDirection':windData_ICON['windDirection'], 
#             }
# =============================================================================
# cloudDict_iconlem, cloudDict_obs
# =============================================================================
#     dictOut = {'cloudMask':cloudMask, 
#                'cloudBase':CB_array,
#                'cloudTop':CT_array, 
#                'liquidCloudFraction':mean_CF_liquid,
#                'iceCloudFraction':mean_CF_ice, 
#                'totalCloudFraction':mean_CF_tot, 
#                'datetimeCloudFraction':datetime_CF, 
#                'heightCloudFraction':height,
#                'duration':duration,
#                'cloudLWP':cloudLWP,
#                'chordLength':chordLength, 
#                'massFlux':massFlux, 
#                'Nclouds':Nclouds,
#             }
# =============================================================================
    dict_iconlem_variables = {
            'IWV_iconlem':IWV_iconlem, 
            'LTS_iconlem':LTS_iconlem,
            'PBLheight_iconlem':PBLheight_iconlem,
            'datetime_iconlem':datetime_ICON,
            }
# =============================================================================
    dict_surface_fluxes = {
            'SHF_iconlem':SHFL_30min.values,
            'LHF_iconlem':LHFL_30min.values,
            'datetime_30m':datetime_30m, 
            'SHF_obs':SHF_obs, 
            'LHF_obs':LHF_obs, 
            'SHF_Err_obs':SHF_Err_obs, 
            'LHF_Err_obs':LHF_Err_obs,
            'LW_iconlem':LW_mod_30min.values,
            'SW_iconlem':SW_mod_30min.values,
            'LW_obs':LW_obs_30min.values,
            'SW_obs':SW_obs_30min.values,
            'SW_Err_obs':SW_obs_Err_30min.values,
            'LW_Err_obs':LW_obs_Err_30min.values,            
            }
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

# flags for activating/deactivating plotting routines
flagPlotThetaVglobalProfiles = 1
flagThermodynVarScatterPlots = 1
flagPlotCloudFractionGlobal  = 1
flagPlotCloudProperties      = 1
flagPlotMeanVarianceW        = 1

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'

#filenameObs = 'dataset_PBLcloudPaper_ModObs_20130502.p'
#filenameMod = 'icon_lem_derivedproperties20130502.nc'

fileListObs               = sorted(glob.glob(pathObs+'*.p'))
fileListMod               = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                    = len(fileListMod)    
date_arr = []    
skewW_mod = []
Hwind_mod = []
w_mod = []
w_obs = []
Hwind_obs = []
skewW_obs = []
timeArr = []
# loop on the ensemble of days (statistics)
for indFile in range(Nfiles):
    
    filenameMod   = fileListMod[indFile]
    filenameObs   = fileListObs[indFile]
    date = fileListMod[indFile][87:95]
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
    w_m = ncdata.groups['Temp_data'].variables['vertWind']
    w_o = new_dict[3]['verticalWind']

    if (w_m.shape != (9600,151)):
         full_m = np.zeros((9600, 151))
         full_m.fill(np.nan)
         full_m[0:len(time), :] = w_m[:,:]
         w_mod.append(full_m)
    else:
        w_mod.append(w_m)
        

    if (w_o.shape[0] != 9600):
        if (w_o.shape[1] != 151):
             full_o = np.zeros((9600, 151))
             full_o.fill(np.nan)
             full_o[0:len(time), :-1] = w_o[:,:]
             w_obs.append(full_o)
        else:
             full_o = np.zeros((9600, 151))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = w_m[:,:]
             w_obs.append(full_o)
    else:
        if (w_o.shape[1] != 151):
             full_o = np.zeros((9600, 151))
             full_o.fill(np.nan)
             full_o[0:len(time), :-1] = w_o[:,:]
             w_obs.append(full_o)
        else:
             full_o = np.zeros((9600, 151))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = w_m[:,:]
             w_obs.append(full_o)
                
    
    # reading horizontal wind from model
    Hwind_m = ncdata.groups['Temp_data'].variables['windSpeed']
    Hwind_o = new_dict[3]['horizontalWind']

    if (Hwind_m.shape[0] != 9600):
        if (Hwind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = Hwind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = Hwind_o[:,1:]
             Hwind_mod.append(full_m)
             Hwind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = Hwind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = Hwind_o[:,:]
             Hwind_mod.append(full_m)
             Hwind_obs.append(full_o)
    else:
        if (Hwind_m.shape[1] != 150):
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = Hwind_m[:,1:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = Hwind_o[:,1:]
             Hwind_mod.append(full_m)
             Hwind_obs.append(full_o)
        else:
             full_m = np.zeros((9600, 150))
             full_m.fill(np.nan)
             full_m[0:len(time), :] = Hwind_m[:,:]
             full_o = np.zeros((9600, 150))
             full_o.fill(np.nan)
             full_o[0:len(time), :] = Hwind_o[:,:]
             Hwind_mod.append(full_m)
             Hwind_obs.append(full_o)
    
#%%

W_allDays_mod = np.dstack(w_mod[:])
Hwind_allDays_mod = np.dstack(Hwind_mod[:])
skewW_allDays_mod = np.dstack(skewW_mod[:])
W_allDays_obs = np.dstack(w_obs[:])
Hwind_allDays_obs = np.dstack(Hwind_obs[:])
skewW_allDays_obs = np.dstack(skewW_obs[:])

# calculating mean values 
meanW_mod = np.nanmean(W_allDays_mod, axis=2)
meanW_obs = np.nanmean(W_allDays_obs, axis=2)

meanHwind_mod = np.nanmean(Hwind_allDays_mod, axis=2)
meanHwind_obs = np.nanmean(Hwind_allDays_obs, axis=2)

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
fig, ax = plt.subplots(3, 2, figsize=(14,10), sharex=True, sharey=True)
ax1 = plt.subplot(3,2,1) 
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

ax2 = plt.subplot(3,2,2) 
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
cbar.set_label(label="skewness W",size=16)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
#plt.savefig(pathDebugFig+'meridionalWind_iconlem_'+date+'.png', format='png')
plt.tight_layout()

ax3 = plt.subplot(3,2,3) 
meanW_obs_plot =  np.ma.array(meanW_obs[:,:-1], mask=np.isnan(meanW_obs[:,:-1]))

ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
ax3.get_xaxis().tick_bottom()  
ax3.get_yaxis().tick_left()
ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax3.xaxis_date()
cax = ax3.pcolormesh(timeArr[0], height, meanW_obs_plot.transpose(), vmin=-1.5, vmax=1.5, cmap='seismic')
ax3.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax3.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax3.set_title("mean vertical velocity W ICON_LEM", fontsize=16)#
ax3.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax3.set_ylabel("height [m]", fontsize=fontSizeY)

plt.tight_layout()
ax4 = plt.subplot(3,2,4) 
meanW_mod_plot =  np.ma.array(meanW_mod[:,:-1], mask=np.isnan(meanW_mod[:,:-1]))
ax4.spines["top"].set_visible(False)  
ax4.spines["right"].set_visible(False)  
ax4.get_xaxis().tick_bottom()  
ax4.get_yaxis().tick_left()
ax4.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax4.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax4.xaxis_date()
cax = ax4.pcolormesh(timeArr[0], height, meanW_mod_plot.transpose(), vmin=-1.5, vmax=1.5, cmap='seismic')
ax4.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax4.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax4.set_title("mean vertical velocity W ICON_LEM", fontsize=16)#
ax4.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax4.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label="w [m/s]",size=1)
cbar.ax.tick_params(labelsize=10)
cbar.aspect=80
plt.tight_layout()


ax5 = plt.subplot(3,2,5) 
Hwind_obs_plot =  np.ma.array(meanHwind_obs[:,:-1], mask=np.isnan(meanHwind_obs[:,:-1]))

ax5.spines["top"].set_visible(False)  
ax5.spines["right"].set_visible(False)  
ax5.get_xaxis().tick_bottom()  
ax5.get_yaxis().tick_left()
ax5.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax5.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax5.xaxis_date()
cax = ax5.pcolormesh(timeArr[0], height, Hwind_obs_plot.transpose(), vmin=0., vmax=10, cmap='YlGnBu')
ax5.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax5.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
#ax5.set_title("mean wind speed obs", fontsize=16)#
ax5.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax5.set_ylabel("height [m]", fontsize=fontSizeY)

plt.tight_layout()
ax6 = plt.subplot(3,2,6) 
Hwind_mod_plot =  np.ma.array(meanHwind_mod[:,:-1], mask=np.isnan(meanHwind_mod[:,:-1]))
ax6.spines["top"].set_visible(False)  
ax6.spines["right"].set_visible(False)  
ax6.get_xaxis().tick_bottom()  
ax6.get_yaxis().tick_left()
ax6.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax6.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax6.xaxis_date()
cax = ax6.pcolormesh(timeArr[0], height, Hwind_mod_plot.transpose(), vmin=0., vmax=10, cmap='YlGnBu')
ax6.set_ylim(107.,Hmax)                                               # limits of the y-axesn  cmap=plt.cm.get_cmap("viridis", 256)
ax6.set_xlim(timeStart, timeEnd)                                 # limits of the x-axes
ax6.set_title("mean wind speed ICON_LEM", fontsize=16)#
ax6.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
ax6.set_ylabel("height [m]", fontsize=fontSizeY)
cbar = fig.colorbar(cax, orientation='vertical')
cbar.set_label(label="wind speed [m/s]",size=16)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
plt.tight_layout()
plt.savefig('/work/cacquist/HDCP2_S2/statistics/figs/patch003/stat_mean_field_w_skew_Hwind.png')
# ============================================================================
#%%
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