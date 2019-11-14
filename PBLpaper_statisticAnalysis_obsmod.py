#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:20:56 2019
@ author: cacquist
@ contact: cacquist@meteo.uni-koeln.de
@  date : 02 august 2019
@  goal : elaborate statistics of the dataset produced with the code 
PBLpaper_prepare_dataset_model_obs_stat.py. It reads as input the files 
dataset_PBLcloudPaper_ModObs_date.p and produces plots of quantities of interest 
listed below.

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
dataset_mean_variance_obs = []
dataset_mean_variance_mod = []
duration_obs              = []
chordLength_obs           = []
massFlux_obs              = []
cloudLWP_obs              = []
LWC_obs                   = []
meanheightFromCB_obs      = []
meanheightFromCB_mod      = []
cloudMeanCB_mod           = []
cloudMeanCT_mod           = []
cloudMeanCB_obs           = []
cloudMeanCT_obs           = []
timeCloudStart_obs        = []
timeCloudEnd_obs          = []
duration_mod              = []
chordLength_mod           = []
massFlux_mod              = []
LWC_mod                   = []
timeCloudStart_mod        = []
timeCloudEnd_mod          = []
dataset_cloudFraction_obs = []
dataset_cloudFraction_mod = []
dataset_cloudFractionLiquid_obs = []
dataset_cloudFractionLiquid_mod = []
dataset_cloudFractionIce_obs    = []
dataset_cloudFractionIce_mod    = []

thetaV_radiosObs = []
thetaV_iconlem   = []
         
P_radiosondes = []
T_radiosondes = []
theta_v_radiosondes = []
height_radiosondes = []     
time_radiosondes = []
theta_v_mod = []
dateArr = []
time_mod = []
height_mod = []
lcl_radiosondes = []
ccl_radiosondes = []
lts_radiosondes =[]
pblHeight_radiosondes = []
lcl_mod = []
ccl_mod = []
lts_mod = []
pblHeight_mod = []
IWV_obs = []
IWV_mod = []
datetime_ICON_arr = []
date_arr = []
LHSF_mod = []
SHSF_mod = []
LHSF_obs = []
SHSF_obs = []
LHSF_err_obs = []
SHSF_err_obs = []
datetime_30m = []
cloudLWP_mod = []
LWF_mod = []
SWF_mod = []
LWF_obs = []
SWF_obs = []
LWF_err_obs = []
SWF_err_obs = []
# ----------------------------------------------------------------------------------------
# ------- Analysis of the mean hourly profiles of variance of vertical velocity for obs. ICON-LEM
# ----------------------------------------------------------------------------------------
# loop on the ensemble of days (statistics)
for indFile in range(Nfiles):
    
    
     print(indFile)
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
     w_mod         = ncdata.groups['Temp_data'].variables['vertWind']
     varW_mod      = ncdata.groups['Temp_data'].variables['varianceW']
     
     
     # opening the file containing all the data
     infile        = open(filenameObs,'rb')
     new_dict      = pickle.load(infile, encoding='latin1')
     W_obs         = new_dict[3]['verticalWind']
 
     timeWindow    = 200 #10 #200 for 9 seconds time window corresponding 
     #to 30 min considering that PBL data have time resolution of 3 minutes
     varianceW_obs = f_calcWvariance(W_obs,time,height,timeWindow)
     print(len(new_dict[8]['cloudLWC']))
     # reading cloud properties data: 
     # duration, chord length, cloud LWP, massflux, cloud fraction
     duration_obs.append(np.asarray(new_dict[9]['duration']))
     duration_mod.append(np.asarray(new_dict[8]['duration']))
     chordLength_obs.append(np.asarray(new_dict[9]['chordLength']))
     chordLength_mod.append(np.asarray(new_dict[8]['chordLength']))
     cloudLWP_obs.append(np.asarray(new_dict[9]['cloudLWP']))
     cloudLWP_mod.append(np.asarray(new_dict[8]['cloudLWP']))    
     massFlux_obs.append(np.asarray(new_dict[9]['massFlux']))
     massFlux_mod.append(np.asarray(new_dict[8]['massFlux']))
     meanheightFromCB_mod.append(np.asarray(new_dict[8]['meanheightFromCB']))
     meanheightFromCB_obs.append(np.asarray(new_dict[9]['meanheightFromCB']))
     cloudMeanCB_mod.append(np.asarray(new_dict[8]['cloudMeanCB']))
     cloudMeanCT_mod.append(np.asarray(new_dict[8]['cloudMeanCT']))
     cloudMeanCB_obs.append(np.asarray(new_dict[9]['cloudMeanCB']))
     cloudMeanCT_obs.append(np.asarray(new_dict[9]['cloudMeanCT']))
     LWC_obs.append(np.asarray(new_dict[9]['cloudLWC']))
     LWC_mod.append(np.asarray(new_dict[8]['cloudLWC']))
     timeCloudStart_obs.append(np.asarray(new_dict[9]['timeCloudStart']))
     timeCloudStart_mod.append(np.asarray(new_dict[8]['timeCloudStart']))
     timeCloudEnd_obs.append(np.asarray(new_dict[9]['timeCloudEnd']))
     timeCloudEnd_mod.append(np.asarray(new_dict[8]['timeCloudEnd']))
     dataset_cloudFraction_obs.append(new_dict[9]['totalCloudFraction'])
     dataset_cloudFraction_mod.append(new_dict[8]['totalCloudFraction'])
     dataset_cloudFractionLiquid_obs.append(new_dict[9]['liquidCloudFraction'])
     dataset_cloudFractionLiquid_mod.append(new_dict[8]['liquidCloudFraction'])
     dataset_cloudFractionIce_obs.append(new_dict[9]['iceCloudFraction'])
     dataset_cloudFractionIce_mod.append(new_dict[8]['iceCloudFraction'])
     dateArr.append(date)
     IWV_obs.append(new_dict[3]['IWV_mwr'])
     IWV_mod.append(new_dict[10]['IWV_iconlem'])
     datetime_ICON_arr.append(new_dict[10]['datetime_iconlem'])
     

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
     
     # reading model data from lem
     theta_v_mod.append(new_dict[5]['virtualPotentialTemperature'])
     time_mod.append(new_dict[5]['time'])
     height_mod.append(new_dict[5]['height'])
     lcl_mod.append(new_dict[5]['lclHeight'])
     ccl_mod.append(new_dict[5]['cclHeight'])
     lts_mod.append(new_dict[10]['LTS_iconlem'])
     LHSF_mod.append(new_dict[11]['LHF_iconlem'])
     SHSF_mod.append(new_dict[11]['SHF_iconlem'])
     LHSF_obs.append(new_dict[11]['LHF_obs'])
     SHSF_obs.append(new_dict[11]['SHF_obs'])
     LHSF_err_obs.append(new_dict[11]['LHF_Err_obs'])
     SHSF_err_obs.append(new_dict[11]['SHF_Err_obs'])     
     datetime_30m.append(new_dict[11]['datetime_30m'])
     LWF_mod.append(new_dict[11]['LW_iconlem'])
     SWF_mod.append(new_dict[11]['SW_iconlem'])
     LWF_obs.append(new_dict[11]['LW_obs'])
     SWF_obs.append(new_dict[11]['SW_obs'])
     LWF_err_obs.append(new_dict[11]['LW_Err_obs'])
     SWF_err_obs.append(new_dict[11]['SW_Err_obs'])
     pblHeight_mod.append(new_dict[7]['PBLHeight'])
     # ----------------------------------------------------------------------------------------
     # ------- Analysis of the mean hourly profiles of variance of vertical velocity for obs. ICON-LEM, ICON-INSCAPE
     # ----------------------------------------------------------------------------------------
     #---- calculating mean variance and standard deviation profiles for each hour of the day for obs and model
     print('calculating mean variance and standard deviation profiles for each hour of the day for obs and models')
     varianceWmean_obs = f_calcMeanStdVarProfiles(varianceW_obs, time[:], height[:],\
                                     date, yy, mm, dd, NprofilesOut, timeIncrement) 
     varianceWmean_mod = f_calcMeanStdVarProfiles(varW_mod[:,:], time[:], height[:], \
                                     date, yy, mm, dd, NprofilesOut, timeIncrement) 
 
     dataset_mean_variance_obs.append(varianceWmean_obs)
     dataset_mean_variance_mod.append(varianceWmean_mod)
     
     #print(LHSF_mod)
     #print(new_dict[11]['LHF_iconlem'])

#%%
# processing of LWC profile: goal is to derive hourly mean LWC profile from the collected clouds for model and obs
LWC_All_obs              = np.concatenate(LWC_obs, axis=0)
LWC_All_mod              = np.concatenate(LWC_mod, axis=0)
duration_all_obs         = np.concatenate(duration_obs, axis=0)
duration_all_mod         = np.concatenate(duration_mod, axis=0)
CB_all_obs               = np.concatenate(LWC_obs, axis=0)
meanheightFromCB_All_obs = np.concatenate(meanheightFromCB_obs, axis=0)
meanheightFromCB_All_mod = np.concatenate(meanheightFromCB_mod, axis=0)
cloudMeanCT_All_obs      = np.concatenate(cloudMeanCB_obs, axis=0)
cloudMeanCB_All_obs      = np.concatenate(cloudMeanCB_obs, axis=0)
cloudMeanCT_All_mod      = np.concatenate(cloudMeanCT_mod, axis=0)
cloudMeanCB_All_mod      = np.concatenate(cloudMeanCB_mod, axis=0)
#duration_intervals = [0, 400, 800, 1600, 3200]
duration_intervals       = [0, 200, 800, 1600]

LWC_obs_DF              = pd.DataFrame(LWC_All_obs, index=duration_all_obs, columns=height)
LWC_mod_DF              = pd.DataFrame(LWC_All_mod, index=duration_all_mod, columns=height)
meanHeightFromCB_obs_DF = pd.DataFrame(meanheightFromCB_All_obs, index=duration_all_obs, columns=height)
meanHeightFromCB_mod_DF = pd.DataFrame(meanheightFromCB_All_mod, index=duration_all_mod, columns=height)
CB_obs_DF               = pd.Series(cloudMeanCB_All_obs, index=duration_all_obs)
CB_mod_DF               = pd.Series(cloudMeanCB_All_mod, index=duration_all_mod)
CT_obs_DF               = pd.Series(cloudMeanCT_All_obs, index=duration_all_obs)
CT_mod_DF               = pd.Series(cloudMeanCT_All_mod, index=duration_all_mod)


#%%
# calculating mean profiles: concept is to select all profiles within the \
# specific duration interval. Then, we select all heights between 0 and 1 and we loop on them:
# for every height 
height_grid = np.linspace(0., 1., num=5)
Nprofiles_durationIntervalsObs = []
Nprofiles_durationIntervalsMod = []

LWC_mean_rescaled_obs = []
LWC_mean_rescaled_mod = []
CB_obs_distr_arr = []
CB_mod_distr_arr = []
CT_obs_distr_arr = []
CT_mod_distr_arr = []
for indDuration in range(0, len(duration_intervals)):
    #print(indDuration)
    #print(indDuration-len(duration_intervals))
    if ((indDuration-len(duration_intervals)+1) == 0.):
        mask_duration_obs    = LWC_obs_DF.index > duration_intervals[indDuration]
        mask_duration_mod    = LWC_mod_DF.index > duration_intervals[indDuration]
        #print('interval sampled > '+str(duration_intervals[indDuration])+' s')
    else:
        mask_duration_obs    = (LWC_obs_DF.index > duration_intervals[indDuration]) * (LWC_obs_DF.index <= duration_intervals[indDuration+1])
        mask_duration_mod    = (LWC_mod_DF.index > duration_intervals[indDuration]) * (LWC_mod_DF.index <= duration_intervals[indDuration+1])
        #print('interval sampled'+str(duration_intervals[indDuration])+'-'+str(duration_intervals[indDuration+1])+' s')
    
    # selecting the relative heights and the LWC corresponding profiles of the selected profiles     
    relHeightCBCT_obs_durInt = meanHeightFromCB_obs_DF.loc[mask_duration_obs,:].values
    LWC_obs_durInt           = LWC_obs_DF.loc[mask_duration_obs,:].values
    relHeightCBCT_mod_durInt = meanHeightFromCB_mod_DF.loc[mask_duration_mod,:].values
    LWC_mod_durInt           = LWC_mod_DF.loc[mask_duration_mod,:].values    
    # reading the number of profiles corresponding to the selected duration 
    # interval and storing the value
    N_obs                    = np.shape(relHeightCBCT_obs_durInt)[0]
    Nprofiles_durationIntervalsObs.append(N_obs)
    N_mod                    = np.shape(relHeightCBCT_mod_durInt)[0]
    Nprofiles_durationIntervalsMod.append(N_mod)    
    
    # building array to use as index given by the duration selected in the picked interval
    durationSelected_obs     = duration_all_obs[mask_duration_obs]
    durationSelected_mod     = duration_all_mod[mask_duration_mod]
    
    # defining pandas structures for filtering and averaging the data
    RescaledLWCMatrixObs     = np.zeros((N_obs,len(height_grid)))
    RescaledLWCMatrixMod     = np.zeros((N_mod,len(height_grid)))
    
    # loop on heights: for every height and every profile, we check the heights 
    # between the two height grid points and we average them and store them in the rescaled matrix
    for indHeight in range(len(height_grid)-1):
        Hmin = height_grid[indHeight]
        Hmax = height_grid[indHeight+1]
        
        for indTimeProf in range(len(durationSelected_obs)):
            ProfSel_obs                                 = pd.Series(LWC_obs_durInt[indTimeProf,:], \
                                                index=relHeightCBCT_obs_durInt[indTimeProf,:])
            mask_h                                      = (ProfSel_obs.index >= Hmin) * (ProfSel_obs.index <= Hmax)
            RescaledLWCMatrixObs[indTimeProf,indHeight] = np.nanmean(ProfSel_obs.loc[mask_h])

        for indTimeProf in range(len(durationSelected_mod)):
            ProfSel_mod                                 = pd.Series(LWC_mod_durInt[indTimeProf,:], \
                                                index=relHeightCBCT_mod_durInt[indTimeProf,:])
            mask_h                                      = (ProfSel_mod.index >= Hmin) * (ProfSel_mod.index <= Hmax)
            RescaledLWCMatrixMod[indTimeProf,indHeight] = np.nanmean(ProfSel_mod.loc[mask_h])
    
    
    LWC_mean_rescaled_obs.append(RescaledLWCMatrixObs)
    LWC_mean_rescaled_mod.append(RescaledLWCMatrixMod)
    
    
    # calculating CB/CT distributions for each of the duration types
    CB_obs_distr_arr.append(CB_obs_DF.loc[mask_duration_obs])
    CB_mod_distr_arr.append(CB_mod_DF.loc[mask_duration_mod])
    CT_obs_distr_arr.append(CT_obs_DF.loc[mask_duration_obs])
    CT_mod_distr_arr.append(CT_mod_DF.loc[mask_duration_mod])
    
# calculating mean profiles for observations on the rescaled grid
LWC_mean_prof_obs = np.zeros((len(duration_intervals),len(height_grid)))
LWC_std_prof_obs  = np.zeros((len(duration_intervals),len(height_grid)))

for indDurInt in range(len(LWC_mean_rescaled_obs)):
    LWC_mean_prof_obs[indDurInt,:] = np.nanmean(LWC_mean_rescaled_obs[indDurInt], axis=0)
    LWC_std_prof_obs[indDurInt,:]  = np.nanstd(LWC_mean_rescaled_obs[indDurInt], axis=0)
    
LWC_mean_prof_mod = np.zeros((len(duration_intervals),len(height_grid)))
LWC_std_prof_mod  = np.zeros((len(duration_intervals),len(height_grid)))

for indDurInt in range(len(LWC_mean_rescaled_mod)):
    LWC_mean_prof_mod[indDurInt,:] = np.nanmean(LWC_mean_rescaled_mod[indDurInt], axis=0)
    LWC_std_prof_mod[indDurInt,:]  = np.nanstd(LWC_mean_rescaled_mod[indDurInt], axis=0)   
    
#%%
   
# plotting CB distributions
fig = plt.figure(figsize=(10,14))
plt.gcf().subplots_adjust(bottom=0.15)
stringplot = ['a)','b)','c)','d)']
stringmidHours = ['0-200 s', '200 - 800 s', '800-1600s', '>1600 s']
Rangearr_CB = [[600., 3000.], [600., 3000.],[600., 3000.], [600., 3000.]]
Rangearr_CT = [[1300., 5000.], [1300., 5000.],[1300.,5000.], [1300., 5000.]]
Nbins_CB = 10
Nbins_CT = 10
ymaxArr = [0.0025, 0.0009, 0.0015, 0.025]
xlabelArr = [1400., 5000., 1000., 250.]
ylabelArr = [0.007, 0.0022, 0.0037, 0.056]
colorArr = ['black', 'black', 'black', 'black']
indObs = [1,3,5,7]
indMod = [2,4,6,8]
for indPlot in range(0,4):

    ax = plt.subplot(4,2,indObs[indPlot])  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('cloud base height [m]') 
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')

    plt.hist(CB_obs_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    plt.hist(CB_mod_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    ax = plt.subplot(4,2,indMod[indPlot])  
    ax.set_title(stringmidHours[indPlot])
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('cloud top height [m]')   
    plt.hist(CT_obs_distr_arr[indPlot], \
                      bins=10, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    plt.hist(CT_mod_distr_arr[indPlot], \
                      bins=10, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)  
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')

plt.tight_layout()

#plt.savefig(pathDebugFig+'histograms_cloudProperties_'+date+'.png', format='png')
#%%
# plotting all profiles for each duration interval after regridding
fig, ax       = plt.subplots(nrows=2, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.002, 0.002, 0.0004, 0.001]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr_obs = [1,2,3,4]
indPlotArr_mod = [5,6,7,8]

for ind in range(0, len(LWC_mean_rescaled_obs)):
     xmax          = xmaxArr[ind]
     ax        = plt.subplot(2,len(duration_intervals),indPlotArr_obs[ind])  
     ax.spines["top"].set_visible(False)  
     ax.spines["right"].set_visible(False)  
     ax.get_xaxis().tick_bottom()  
     ax.get_yaxis().tick_left() 
     ax.set_ylim(ymin, ymax)
     ax.set_xlim(xmin, xmax)
     ax.set_xlabel('LWC $Kg m^{-3}$')
     ax.set_ylabel('rel. dist. from cloud base')
     LWC_plot_obs = LWC_mean_rescaled_obs[ind]
     for indProfile in range(np.shape(LWC_mean_rescaled_obs[ind])[0]):
         plt.plot(LWC_plot_obs[indProfile,:], height_grid, color='black')
         
         
     ax1        = plt.subplot(2,len(duration_intervals),indPlotArr_mod[ind])  
     ax1.spines["top"].set_visible(False)  
     ax1.spines["right"].set_visible(False)  
     ax1.get_xaxis().tick_bottom()  
     ax1.get_yaxis().tick_left() 
     ax1.set_ylim(ymin, ymax)
     ax1.set_xlim(xmin, xmax)
     ax1.set_xlabel('LWC $Kg m^{-3}$')
     ax1.set_ylabel('rel. dist. from cloud base')
     LWC_plot_mod = LWC_mean_rescaled_mod[ind]
     for indProfile in range(np.shape(LWC_mean_rescaled_mod[ind])[0]):
         plt.plot(LWC_plot_mod[indProfile,:], height_grid, color='red')        
plt.tight_layout()
#%%
# =============================================================================
fig, ax       = plt.subplots(nrows=1, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.002, 0.002, 0.0004, 0.001]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr = [1,2,3,4]
for ind in range(0, len(LWC_mean_rescaled_obs)):
     xmax          = xmaxArr[ind]
     ax        = plt.subplot(1,len(duration_intervals),indPlotArr[ind])  
     ax.spines["top"].set_visible(False)  
     ax.spines["right"].set_visible(False)  
     ax.get_xaxis().tick_bottom()  
     ax.get_yaxis().tick_left() 
     ax.set_ylim(ymin, ymax)
     ax.set_xlim(xmin, xmax)
     ax.set_xlabel('LWC $Kg m^{-3}$')
     ax.set_ylabel('rel. dist. from cloud base')
     plt.plot(LWC_mean_prof_obs[ind,:], height_grid, color='black', label='obs')
     y1        = LWC_mean_prof_obs[ind,:]-LWC_std_prof_obs[ind,:]
     y2        = LWC_mean_prof_obs[ind,:]+LWC_std_prof_obs[ind,:]
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
     
     plt.plot(LWC_mean_prof_mod[ind,:], height_grid, color='red', label='icon-lem')
     y1        = LWC_mean_prof_mod[ind,:]-LWC_std_prof_mod[ind,:]
     y2        = LWC_mean_prof_mod[ind,:]+LWC_std_prof_mod[ind,:]
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
# =============================================================================

    
#%%

# definitions for plotting the variance of the single day and patch as a reference
# =============================================================================
#      varWmean_obs = varianceWmean_obs['meanProfiles']
#      varWmean_mod = varianceWmean_mod['meanProfiles']
#      varWstd_obs  = varianceWmean_obs['stdProfiles']
#      varWstd_mod  = varianceWmean_mod['stdProfiles']
#      timeMean     = varianceWmean_obs['meanTime']
#      heightMean   = varianceWmean_obs['height']
#      indHourPlotStart = 8
#      plotDoneDay  = f_plotVarianceWSingleDays(date,varWmean_obs,varWmean_mod, \
#                             varWstd_obs,varWstd_mod,indHourPlotStart, height, pathFig)
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
#%%
# =============================================================================
# plotting scatter plots of surface latent and sensible heat flux for every half hour 
# =============================================================================

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
plt.xlim(-350., 50.)
plt.ylim(-350., 50.)
plt.xlabel('Latent heat surface flux obs [W/m^2]', fontsize=16)
plt.ylabel('Latent heat surface flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
for indFile in range(Nfiles):
    print(indFile)
    LHSF_iconlemPlot = LHSF_mod[indFile]
    LHSF_obsPlot = -LHSF_obs[indFile]
    sizeDots = LHSF_err_obs[indFile]
    cax = ax.scatter(LHSF_obsPlot[:], LHSF_iconlemPlot[:-1], c=colors, cmap=cmap, \
                 s=10*sizeDots)
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

#%%

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
                 s=10*sizeDots)
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

# =============================================================================
# plotting scatter plots of longwave ans shortwave radiative fluxes for every half hour 
# =============================================================================
datetime_30m     = [datetime.datetime(int(yy),int(mm),int(dd),0,0,0) + \
                        datetime.timedelta(minutes=30*x) for x in range(0, 49)] 
datetime_fluxes = datetime_30m
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
plt.xlim(225., 400.)
plt.ylim(225., 400.)
plt.xlabel('Longwave downward flux obs [W/m^2]', fontsize=16)
plt.ylabel('Longwave downward flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
for indFile in range(Nfiles):
    print(indFile)
    LWF_iconlemPlot = LWF_mod[indFile]
    LWF_obsPlot = LWF_obs[indFile]
    sizeDots = LWF_err_obs[indFile]
    cax = ax.scatter(LWF_obsPlot[:], LWF_iconlemPlot[:], c=colors, cmap=cmap, \
                 s=10*sizeDots)
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
plt.savefig(pathFig+'LWF_scatterplot_obs_mod_allDays.png', format='png')

#%%
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
plt.xlim(-20., 950.)
plt.ylim(-150., 50.)
plt.xlabel('Shortwave downward flux obs [W/m^2]', fontsize=16)
plt.ylabel('Shortwave downward flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
for indFile in range(Nfiles):
    print(indFile)
    SWF_iconlemPlot = SWF_mod[indFile]
    SWF_obsPlot = SWF_obs[indFile]
    sizeDots = SWF_err_obs[indFile]
    cax = ax.scatter(SWF_obsPlot[:], SWF_iconlemPlot[:], c=colors, cmap=cmap, \
                 s=10*sizeDots)
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
plt.savefig(pathFig+'SWF_scatterplot_obs_mod_allDays.png', format='png')


#%%

# =============================================================================
# calculating and plotting IWV distributions for each hour of the day
# =============================================================================
IWV_mod_nd = np.stack(IWV_mod).T
IWV_obs_nd = np.stack(IWV_obs).T
IWV_obs_nd = np.vstack([IWV_obs_nd,np.nan*np.ones(6)]) 
hours = [0, 6, 9, 12, 15, 18, 24]

def hour2idx(hour, dtime=9):
    return int(hour*3600/dtime)

nbins = 10
ymax  = 0.6
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
xmin = 10.
xmax=25.
indplot = 1
for indHour in range(len(hours)-1):
    hourInf = hours[indHour]
    hourSup = hours[indHour+1]
    hinf_idx = hour2idx(hourInf)
    hsup_idx = hour2idx(hourSup)
    
    mod = IWV_mod_nd[hinf_idx:hsup_idx,:]
    obs = IWV_obs_nd[hinf_idx:hsup_idx,:]
    percentiles_mod = np.nanpercentile(mod, [25, 50, 75, 90])
    percentiles_obs = np.nanpercentile(obs, [25, 50, 75, 90])
    
    #plt.figure()
    #plt.hist(mod.flatten(),range=(10,30), normed=True, alpha=0.5)
    #lt.hist(obs.flatten(),range=(10,30), normed=True, alpha=0.5)
    ax = plt.subplot(2,3,indplot)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()     
    matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
    matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('norm occurrences')
    ax.set_xlabel('IWV [Kg/m^2]')
    #ax.ylim(ymax)
    plt.ylim(0.,ymax)
    plt.xlim(xmin, xmax)
    plt.text(19., ymax-1.5*ymax/10., 'median mod = '+str(round(percentiles_mod[1], 1)))
    plt.text(19., ymax-2.*ymax/10., 'median obs = '+str(round(percentiles_obs[1], 1)))  
    plt.grid(b=True, which='major', color='#666666', linestyle=':')
    plt.hist(mod.flatten(), bins=nbins, normed=True, color='red', cumulative=False, range=[xmin, xmax], alpha=0.5, label='icon-lem')       
    plt.hist(obs.flatten(), bins=nbins, normed=True, color='black', cumulative=False, range=[xmin, xmax], alpha=0.5, label='obs')       
    plt.legend(loc='upper left', fontsize=12, frameon=False)
    ax.set_title(str(hourInf)+' - '+str(hourSup)+' UTC')
    indplot= indplot+1
fig.tight_layout()
plt.savefig(pathFig+'IWV_distrib_mod_obs_stat_global.png', format='png')    





# =============================================================================
# calculating and plotting potential temperature profiles 
# =============================================================================
from myFunctions import f_calculateMeanThetaVModelProfiles
theta_v_dict_obs_mod_arr = f_calculateMeanThetaVModelProfiles(time_radiosondes, \
                                                              theta_v_radiosondes,\
                                                              height_radiosondes, \
                                                              lcl_radiosondes, \
                                                              lts_radiosondes, \
                                                              pblHeight_radiosondes, \
                                                              theta_v_mod, \
                                                              time_mod, \
                                                              height_mod, \
                                                              lcl_mod, \
                                                              lts_mod, \
                                                              pblHeight_mod)


from myFunctions import f_calculateMeanProfilesPlotThetaVRadiosondes
result = f_calculateMeanProfilesPlotThetaVRadiosondes(theta_v_dict_obs_mod_arr, height_mod)
MatrixHourMeanProfileThetaRad = result[0]
MatrixHourStdProfileThetaRad = result[1]
listHourDict = result[2]
gridHeight = height_mod[0]

if flagPlotThetaVglobalProfiles == 1:
    Ncols = 5
    Nrows = 2
    Nplots = 11
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 3000.
    ymin          = height_mod[0][-1]
    xmin          = 280.
    xmax          = 300.
    fontSizeTitle = 16
    fontSizeX     = 12
    fontSizeY     = 12
    #timeTitles = [']
    
    for indPlot in range(1, Nplots):
        ax        = plt.subplot(2,5,indPlot)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left() 
        #ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
        prof_mod = listHourDict[indPlot-1]['meanProfile_mod']
        prof_obs = MatrixHourMeanProfileThetaRad[:, indPlot-1]
        std_mod = listHourDict[indPlot-1]['stdProfileMod']
        labelHour = listHourDict[indPlot-1]['hour']
        std_obs = MatrixHourStdProfileThetaRad[:, indPlot-1]
        plt.plot(prof_obs, gridHeight, label='obs '+str(labelHour)+' UTC',  color='black')
        plt.plot(prof_mod, height_mod[0], label='icon-lem',  color='red')
        y1        = prof_obs-std_obs
        y2        = prof_obs+std_obs
        plt.fill_betweenx(gridHeight, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = prof_mod-std_mod
        y2        = prof_mod+std_mod
        plt.fill_betweenx(height_mod[0], y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        plt.legend(loc='upper right', fontsize=14, frameon=False)
        plt.ylim(ymin,ymax)
        plt.xlim(xmin,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        plt.xlabel('${\Theta_v}$[K]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
    plt.savefig(pathFig+'theta_v_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')


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
for ind in range(len(listHourDict)):
    hour = listHourDict[ind]['hour']
    dim = listHourDict[ind]['n_lcl_hour']
    hourList.append(np.repeat(hour, dim))
    LCLobs.append(listHourDict[ind]['lcl_rad_hour'])
    LCLmod.append(listHourDict[ind]['lcl_mod_hour'])
    PBLobs.append(listHourDict[ind]['pblHeight_rad_hour'])
    PBLmod.append(listHourDict[ind]['pblHeight_mod_hour'])
    LTSobs.append(listHourDict[ind]['lts_rad_hour'])
    LTSmod.append(listHourDict[ind]['lts_mod_hour'])
    
    
hourList = np.concatenate(hourList)
LCLobs = np.concatenate(LCLobs)
LCLmod = np.concatenate(LCLmod)
LTSobs = np.concatenate(LTSobs)
LTSmod = np.concatenate(LTSmod)
PBLobs = np.concatenate(PBLobs)
PBLmod = np.concatenate(PBLmod)


if flagThermodynVarScatterPlots == 1:
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
    colors = hourList
    plt.plot(data, data, color='black', linestyle=':')
    plt.xlim(0., 2000.)
    plt.ylim(0., 2000.)
    plt.xlabel('LCL obs [m]', fontsize=16)
    plt.ylabel('LCL mod [m]', fontsize=16)
    cmap = plt.cm.get_cmap('bwr', 19) 
    cax = ax.scatter(LCLobs, LCLmod, c=colors, cmap=cmap, s=100, vmin=5, vmax=23.)
    plt.grid(b=True, which='major', color='#666666', linestyle=':')
    fig.colorbar(cax, label='time [hh]', cmap=cmap)
    plt.tight_layout()
    plt.savefig(pathFig+'lcl_scatterplot_obs_mod.png', format='png')
    
    
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
    plt.xlim(0., 2000.)
    plt.ylim(0., 2000.)
    plt.xlabel('PBL height obs [m]', fontsize=16)
    plt.ylabel('PBL height mod [m]', fontsize=16)
    cmap = plt.cm.get_cmap('bwr', 19) 
    cax = ax.scatter(PBLobs, PBLmod, c=colors, cmap=cmap, s=100, vmin=5, vmax=23.)
    fig.colorbar(cax, label='time [hh]', cmap=cmap)
    plt.tight_layout()
    plt.savefig(pathFig+'pblh_scatterplot_obs_mod.png', format='png')
    


# =============================================================================
# calculating and plotting global cloud fraction
# =============================================================================

print('calculating mean profiles of cloud fraction every 30 minutes for the whole dataset')
CFmean_obs = np.zeros((np.shape(dataset_cloudFraction_obs)[1], len(height)))
CFmean_mod = np.zeros((np.shape(dataset_cloudFraction_mod)[1], len(height)))
CFmeanLiquid_obs = np.zeros((np.shape(dataset_cloudFraction_obs)[1], len(height)))
CFmeanLiquid_mod = np.zeros((np.shape(dataset_cloudFraction_mod)[1], len(height)))
CFmeanIce_obs = np.zeros((np.shape(dataset_cloudFraction_obs)[1], len(height)))
CFmeanIce_mod = np.zeros((np.shape(dataset_cloudFraction_mod)[1], len(height)))

for indHour in range(np.shape(dataset_cloudFraction_mod)[1]):
    CFmean_obs[indHour,:] = np.mean(np.asarray(dataset_cloudFraction_obs)[:,indHour,:], axis=0)
    CFmean_mod[indHour,:] = np.mean(np.asarray(dataset_cloudFraction_mod)[:,indHour,:], axis=0)
    CFmeanLiquid_obs[indHour,:] = np.mean(np.asarray(dataset_cloudFractionLiquid_obs)[:,indHour,:], axis=0)
    CFmeanLiquid_mod[indHour,:] = np.mean(np.asarray(dataset_cloudFractionLiquid_mod)[:,indHour,:], axis=0)
    CFmeanIce_obs[indHour,:] = np.mean(np.asarray(dataset_cloudFractionIce_obs)[:,indHour,:], axis=0)
    CFmeanIce_mod[indHour,:] = np.mean(np.asarray(dataset_cloudFractionIce_mod)[:,indHour,:], axis=0)    
    
datetime_CF = new_dict[8]['datetimeCloudFraction']

if flagPlotCloudFractionGlobal == 1:
    from myFunctions import f_plotCloudFraction
    cloudFractionPlotsDone = f_plotCloudFraction(datetime_CF, height, pathFig, \
                                                  CFmean_mod, \
                                                  CFmean_obs, \
                                                  CFmeanLiquid_mod, \
                                                  CFmeanLiquid_obs, \
                                                  CFmeanIce_mod, \
                                                  CFmeanIce_obs)
    
    
    
# =============================================================================
# calculating and plotting cloud properties
# =============================================================================


if flagPlotCloudProperties == 1: 
    nbins      = 20
    ymax       = 0.2
    fig = plt.figure(figsize=(15,7))
    plt.gcf().subplots_adjust(bottom=0.15)
    #fig.tight_layout()
    stringplot = ['a)','b)','c)','d)']
    stringmidHours = ['cloud duration [s]', 'chord length [m]', 'Mass flux [Kg/ms]', 'mean LWP [g/m^2]']
    HistQuantities_obs = [np.concatenate(duration_obs), np.concatenate(chordLength_obs), \
                          np.concatenate(massFlux_obs), np.concatenate(np.asarray(cloudLWP_obs))*1000.]
    HistQuantities_iconlem = [np.concatenate(duration_mod), np.concatenate(chordLength_mod), \
                          np.concatenate(massFlux_mod), np.concatenate(np.asarray(cloudLWP_mod))*1000.]
    
    Rangearr = [[0., 1500.], [0., 6000.], [-1500., 1500.], [0., 300.]]
    ymaxArr = [0.0025, 0.0009, 0.0015, 0.025]
    xlabelArr = [1400., 5000., 1000., 250.]
    ylabelArr = [0.007, 0.0022, 0.0037, 0.056]
    colorArr = ['blue', 'red', 'green', 'purple']
    for indPlot in range(1,5):
        ax = plt.subplot(2,2,indPlot)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()    
        ax.get_yaxis().tick_left() 
        #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
        #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
        #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
        ax.set_ylabel('occurrences', fontsize=15 )
        ax.set_xlabel(stringmidHours[indPlot-1], fontsize=15)   
        plt.hist(HistQuantities_obs[indPlot-1], \
                 bins=nbins, \
                 normed=True, \
                 color='black', \
                 range=Rangearr[indPlot-1], \
                 cumulative=False, \
                 alpha=0.5, label='obs')    
        plt.hist(HistQuantities_iconlem[indPlot-1], \
                 bins=nbins, \
                 normed=True, \
                 color=colorArr[indPlot-1], \
                 range=Rangearr[indPlot-1], \
                 cumulative=False, \
                 alpha=0.5, label='icon-lem')  
        plt.legend(frameon=False)
        plt.tight_layout()
    plt.savefig(pathFig+'histograms_cloudProperties_obs_mod.png', format='png')

print('Plotting hourly profiles of variance of vertical velocity during the day for obs, ICON, ICON_INSCAPE')
 


#=============================================================================
# calculating and plotting global mean profiles and stds of variance of vertical velocity
#=============================================================================

matrixVar_obs    = np.zeros((len(height),len(dataset_mean_variance_obs), NprofilesOut))
matrixVar_mod    = np.zeros((len(height),len(dataset_mean_variance_mod), NprofilesOut))

meanVariance_obs = np.zeros((len(height),NprofilesOut))
for indHour in range(NprofilesOut):
    for indDay in range(len(dataset_mean_variance_obs)):
        matrixVar_obs[:,indDay,indHour] = dataset_mean_variance_obs[indDay]['meanProfiles'][:,indHour]
        matrixVar_mod[:,indDay,indHour] = dataset_mean_variance_mod[indDay]['meanProfiles'][:,indHour]
meanVariance_obs = np.nanmean(matrixVar_obs, axis=1)
stdVariance_obs  = np.nanstd(matrixVar_obs, axis=1)
meanVariance_mod = np.nanmean(matrixVar_mod, axis=1)
stdVariance_mod  = np.nanstd(matrixVar_mod, axis=1)

if flagPlotCloudProperties == 1:
    Ncols = 5
    Nrows = 2
    Nplots = 11
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 3000.
    ymin          = 107.
    xmax          = 1.5
    fontSizeTitle = 16
    fontSizeX     = 10
    fontSizeY     = 10
    #timeTitles = [']
    indHourPlotStart = 8 # corresponding to starting at 8 UTC 
    
    for indPlot in range(1, Nplots):
        ax        = plt.subplot(2,5,indPlot)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left() 
        #ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        plt.plot(meanVariance_obs[:,indHourPlotStart], height, label='obs',  color='black')
        plt.plot(meanVariance_mod[:,indHourPlotStart], height, label='icon-lem',  color='red')
        y1        = meanVariance_obs[:,indHourPlotStart]-stdVariance_obs[:,indHourPlotStart]
        y2        = meanVariance_obs[:,indHourPlotStart]+stdVariance_obs[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = meanVariance_mod[:,indHourPlotStart]-stdVariance_mod[:,indHourPlotStart]
        y2        = meanVariance_mod[:,indHourPlotStart]+stdVariance_mod[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        plt.legend(loc='upper right', fontsize=14, frameon=False)
        plt.ylim(ymin,ymax)
        plt.xlim(0.,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        plt.xlabel('${\sigma}$ [m/s]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart+1
    
    plt.savefig(pathFig+'varW_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')
