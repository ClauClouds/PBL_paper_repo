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
    outputArray   = [radiosondeList, \           0
                     tower_dict, \               1
                     dictCosmo, \                2
                     dictObsWindLidarMwr, \      3
                     Thermodyn_cosmo, \          4
                     Thermodyn_iconlem, \        5
                     Thermodyn_obs, \            6
                     cloudDict_iconlem, \        7 
                     cloudDict_obs, \            8
                     dict_iconlem_variables, \   9 
                     dict_surface_fluxes]        10
# Below we provide the definitions of each of them:
# =============================================================================
# radiosondeList: index 0
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
# tower_dict: index 1
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
# dictCosmo: index 2
# =============================================================================
#     dictCosmo = {'pressure':P_cosmo_res.values.transpose(),
#                  'temperature':T_cosmo_res.values.transpose(),
#                  'absoluteHumidity':Q_cosmo_res.values.transpose(),
#                  'cloudnetCategorization':cloudnet_res.data.transpose()
#                  }
# 
# =============================================================================
# dictObsWindLidarMwr: index 3
# =============================================================================
#    dictObsWindLidarMwr = {
#            'datetime_obs':datetime_lwp_iwv_joyce,
#            'verticalWind':w_obs_res.values.transpose(),
#            'horizontalWind':Hwind_obs_res.values.transpose(),
#            'zonalWind':u_obs_res.values.transpose(),
#            'meridionalWind':v_obs_res.values.transpose(), 
#            'skewnessW':skew_obs_res.values.transpose(),
#            'PBLclassObs':PBLclass_obs_res.values.transpose(),
#            'shear':shear_obs_res.values.transpose(),
#            'windDirection':wDir_obs_res.values.transpose(),
#            'mixingLayerHeight_w_obs': zmlaw_obs_res.values.transpose(),
#            'mixingLayerHeightErr_w_obs': zmlawErr_obs_res.values.transpose(),
#            'absoluteHumidity':qProf_obs_res.values.transpose(),
#            'temperature':tProf_obs_res.values.transpose(),
#            'IWV_mwr':IWV_obs_res, 
#            'LWP_mwr':LWP_obs_res,
#            'LWC_cloudnet':LWC_obs_res,
#            'height':height_ICON
#            }

# =============================================================================
# Thermodyn_cosmo, Thermodyn_iconlem, Thermodyn_obs:  index 4, 5, 6
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
#                   'potentialTemperature':Theta,
#                   'time': time,
#                   'height':height,
#                   'LTS':LTS,
#                   }
# =============================================================================
# cloudDict_iconlem, cloudDict_obs:  index 7, 8
# =============================================================================
#     dictOut = {'cloudMask':cloudMask, 
               'cloudBase':CB_array,
               'cloudTop':CT_array, 
               'Ncloudlayers':NcloudLayers,
               'liquidCloudFraction':mean_CF_liquid,
               'iceCloudFraction':mean_CF_ice, 
               'totalCloudFraction':mean_CF_tot, 
               'datetimeCloudFraction':datetime_CF, 
               'heightCloudFraction':height,
               'duration':duration,
               'cloudLWP':cloudLWP,
               'cloudLWC':cloudLWC,
               'meanheightFromCB':meanheightFromCB,
               'cloudMeanCB':meanCB,
               'cloudMeanCT':meanCT,              
               'cloudTimeEnd':cloudTimeEnd,
               'chordLength':chordLength, 
               'massFlux':massFlux, 
               'timeCloudStart':cloudTimeStart, 
               'timeCloudEnd':cloudTimeEnd,
               'Nclouds':Nclouds,
               'LWPall':LWP,
               'cloudThicknessAll':cloudThickness, 
               'UpdraftCBAll':UpdraftCB, 
               'cloudMaturity':cloudMaturity
#             }
# =============================================================================
dict_iconlem_variables: index 9
# =============================================================================
#     dict_out = {
            'IWV_iconlem':IWV_iconlem, 
            'LTS_iconlem':LTS_iconlem,
            'PBLheightRN_iconlem':PBLheightRN_iconlem,
            'PBLheightTW_iconlem':PBLheightTW_iconlem,
            'datetime_iconlem':datetime_ICON,
            'T_iconlem':Tmatrix_iconlem, 
            'u_iconlem':u_iconlem,
            'v_iconlem':v_iconlem,
            'w_iconlem':w_iconlem,
            'varianceW':varianceW_icon_lem, 
            'stdW':stdW_icon_lem,
            'PBLHeightRN':PBLheightRN_iconlem,
            'PBLHeightTW':PBLheightTW_iconlem,
            'windSpeed':windSpeed_icon_lem,#windData_ICON['windSpeed'], 
            'windDirection':windDirection_icon_lem, #windData_ICON['windDirection'], 
            }
# =============================================================================
    dict_surface_fluxes: index 10
# =============================================================================
#      dict_out = {
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
from myFunctions import f_closest



try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle


# setting parameters for calculating averaging and domain size of the model:
NprofilesOut                    = 24  # hourly means
timeIncrement                   = 60  # hourly means
patch                           = 'patch003'

# flags for activating/deactivating plotting routines
flagPlotThetaVglobalProfiles    = 1
flagThermodynVarScatterPlots    = 1
flagPlotCloudFractionGlobal     = 1
flagPlotCloudProperties         = 1
flagPlotMeanVarianceW           = 1

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'


#filenameObs = 'dataset_PBLcloudPaper_ModObs_20130502.p'
#filenameMod = 'icon_lem_derivedproperties20130502.nc'

fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)
dataset_mean_variance_obs       = []
dataset_mean_variance_mod       = []


# cloud variables (mean of cloud units)
cloudUpdraftCB_mod              = []
cloudUpdraftCB_obs              = []
cloudMaturity_mod               = []
cloudMaturity_obs               = []
datetime_ICON_arr               = []
duration_obs                    = []
chordLength_obs                 = []
cloudMassFlux_obs               = []
LWC_obs                         = []
duration_mod                    = []
chordLength_mod                 = []
cloudMassFlux_mod               = []
LWC_mod                         = []
timeCloudStart_mod              = []
timeCloudEnd_mod                = []
cloudLWP_obs                    = []
cloudLWP_mod                    = []
meanheightFromCB_obs            = []
meanheightFromCB_mod            = []
cloudMeanCB_mod                 = []
cloudMeanCT_mod                 = []
cloudMeanCB_obs                 = []
cloudMeanCT_obs                 = []
cloudThickness_mod              = []
cloudThickness_obs              = []
cloudUpdraftCB_mod              = []
cloudUpdraftCB_obs              = []
timeCloudStart_obs              = []
timeCloudEnd_obs                = []


# cloud fraction variables
dataset_cloudFraction_obs       = []
dataset_cloudFraction_mod       = []
dataset_cloudFractionLiquid_obs = []
dataset_cloudFractionLiquid_mod = []
dataset_cloudFractionIce_obs    = []
dataset_cloudFractionIce_mod    = []

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


# global variables
IWV_obs                         = []
IWV_mod                         = []
LWP_all_obs                     = []
LWP_all_mod                     = []
cloudThicknessAll_obs           = []
cloudThicknessAll_mod           = []
updraftCBall_mod                = []
updraftCBall_obs                = []
LHSF_mod                        = []
SHSF_mod                        = []
LHSF_obs                        = []
SHSF_obs                        = []
LHSF_err_obs                    = []
SHSF_err_obs                    = []
datetime_30m                    = []
LWF_mod                         = []
SWF_mod                         = []
LWF_obs                         = []
SWF_obs                         = []
LWF_err_obs                     = []
SWF_err_obs                     = []


# dinamical properties
varianceW_obs                   = []
varianceW_mod                   = []
stdW_obs                        = []
stdW_mod                        = []
skewnessW_obs                   = []
skewnessW_mod                   = []

# ----------------------------------------------------------------------------------------
# loop on the ensemble of days (statistics)
# ----------------------------------------------------------------------------------------

for indFile in range(Nfiles):
    
    
     print(indFile)
     filenameMod   = fileListMod[indFile]
     filenameObs   = fileListObs[indFile]
     date          = fileListMod[indFile][87:95]
     yy            = int(date[0:4])
     mm            = int(date[4:6])
     dd            = int(date[6:8])
     
     date_arr.append(date)
     print('processing date '+date)
     
     # reading time and height from ncdf file (grid of ICON LEM 
     #( ( sec resolution, 9600 time stamps, and 150 height levels)))
     ncdata        = Dataset(filenameMod, mode='r')
     time          = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                     ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
     height        = ncdata.groups['Temp_data'].variables['height2'][:]
     #w_mod         = ncdata.groups['Temp_data'].variables['vertWind']
     #varW_mod      = ncdata.groups['Temp_data'].variables['varianceW']
     
     
     
     # opening the file containing all the data
     infile        = open(filenameObs,'rb')
     new_dict      = pickle.load(infile, encoding='latin1')
          
     # dinamic properties
     varianceW_obs.append(new_dict[3]['varianceW'])
     varianceW_mod.append(new_dict[9]['varianceW'])
     stdW_obs.append(new_dict[3]['stdW'])
     stdW_mod.append(new_dict[9]['stdW'])
     skewnessW_obs.append(new_dict[3]['skewnessW'])
     skewnessW_mod.append(ncdata.groups['Temp_data'].variables['skewnessW'][:])
          
     # duration, chord length, cloud LWP, massflux, cloud fraction
     duration_obs.append(np.asarray(new_dict[8]['duration']))
     duration_mod.append(np.asarray(new_dict[7]['duration']))
     chordLength_obs.append(np.asarray(new_dict[8]['chordLength']))
     chordLength_mod.append(np.asarray(new_dict[7]['chordLength']))
     cloudLWP_obs.append(np.asarray(new_dict[8]['cloudLWP']))
     cloudLWP_mod.append(np.asarray(new_dict[7]['cloudLWP']))    
     cloudMassFlux_obs.append(np.asarray(new_dict[8]['massFlux']))
     cloudMassFlux_mod.append(np.asarray(new_dict[7]['massFlux']))
     meanheightFromCB_mod.append(np.asarray(new_dict[7]['meanheightFromCB']))
     meanheightFromCB_obs.append(np.asarray(new_dict[8]['meanheightFromCB']))
     cloudMeanCB_mod.append(np.asarray(new_dict[7]['cloudMeanCB']))
     cloudMeanCT_mod.append(np.asarray(new_dict[7]['cloudMeanCT']))
     cloudMeanCB_obs.append(np.asarray(new_dict[8]['cloudMeanCB']))
     cloudMeanCT_obs.append(np.asarray(new_dict[8]['cloudMeanCT']))
     LWC_obs.append(np.asarray(new_dict[8]['cloudLWC']))
     LWC_mod.append(np.asarray(new_dict[7]['cloudLWC']))
     timeCloudStart_obs.append(np.asarray(new_dict[8]['timeCloudStart']))
     timeCloudStart_mod.append(np.asarray(new_dict[7]['timeCloudStart']))
     timeCloudEnd_obs.append(np.asarray(new_dict[8]['timeCloudEnd']))
     timeCloudEnd_mod.append(np.asarray(new_dict[7]['timeCloudEnd']))
     
     # reading 30 min averaged cloud fractions for ice, liquid and both
     dataset_cloudFraction_obs.append(new_dict[8]['totalCloudFraction'])
     dataset_cloudFraction_mod.append(new_dict[7]['totalCloudFraction'])
     dataset_cloudFractionLiquid_obs.append(new_dict[8]['liquidCloudFraction'])
     dataset_cloudFractionLiquid_mod.append(new_dict[7]['liquidCloudFraction'])
     dataset_cloudFractionIce_obs.append(new_dict[8]['iceCloudFraction'])
     dataset_cloudFractionIce_mod.append(new_dict[7]['iceCloudFraction'])
     
     
     # reading time series on icon time resolution.
     if (len(new_dict[9]['IWV_iconlem']) == 9600):
         IWV_mod.append(new_dict[9]['IWV_iconlem'])
         LWP_all_mod.append(new_dict[7]['LWPall'])
         cloudThicknessAll_mod.append(new_dict[7]['cloudThicknessAll'])
         updraftCBall_mod.append(new_dict[7]['UpdraftCBAll'])
     else:
         IWVarraymod = np.zeros((9600))
         IWVarraymod.fill(np.nan)
         LWParraymod = np.zeros((9600))
         LWParraymod.fill(np.nan)
         thicknessarraymod = np.zeros((9600))
         thicknessarraymod.fill(np.nan)
         maturityarraymod = np.zeros((9600))
         maturityarraymod.fill(np.nan)
         updraftCBarraymod = np.zeros((9600))
         updraftCBarraymod.fill(np.nan)
         
         IWVarraymod[0:len(new_dict[9]['IWV_iconlem'])] = new_dict[9]['IWV_iconlem']
         LWParraymod[0:len(new_dict[7]['LWPall'])] = new_dict[7]['LWPall']
         thicknessarraymod[0:len(new_dict[7]['cloudThicknessAll'])] = new_dict[7]['cloudThicknessAll']
         updraftCBarraymod[0:len(new_dict[7]['UpdraftCBAll'])] = new_dict[7]['UpdraftCBAll']
         
         IWV_mod.append(IWVarraymod)
         LWP_all_mod.append(LWParraymod)
         cloudThicknessAll_mod.append(thicknessarraymod)
         updraftCBall_mod.append(updraftCBarraymod)
         
     if (len(new_dict[3]['IWV_mwr']) == 9600):
         IWV_obs.append(new_dict[3]['IWV_mwr'])
         LWP_all_obs.append(new_dict[3]['LWP_mwr'])
         cloudThicknessAll_obs.append(new_dict[8]['cloudThicknessAll'])
         updraftCBall_obs.append(new_dict[8]['UpdraftCBAll'])

     else:
         IWVarrayobs = np.zeros((9600))
         IWVarrayobs.fill(np.nan)
         LWParrayobs = np.zeros((9600))
         LWParrayobs.fill(np.nan)
         thicknessarrayobs = np.zeros((9600))
         thicknessarrayobs.fill(np.nan)
         maturityarrayobs = np.zeros((9600))
         maturityarrayobs.fill(np.nan)
         updraftCBarrayobs = np.zeros((9600))
         updraftCBarrayobs.fill(np.nan)
         
         IWVarrayobs[0:len(new_dict[3]['IWV_mwr'])] = new_dict[3]['IWV_mwr']
         LWParrayobs[0:len(new_dict[3]['LWP_mwr'])] = new_dict[3]['LWP_mwr']
         thicknessarrayobs[0:len(new_dict[8]['cloudThicknessAll'])] = new_dict[8]['cloudThicknessAll']
         updraftCBarrayobs[0:len(new_dict[8]['UpdraftCBAll'])] = new_dict[8]['UpdraftCBAll']
         
         
         IWV_obs.append(IWVarrayobs)
         LWP_all_obs.append(LWParrayobs)
         cloudThicknessAll_obs.append(thicknessarrayobs)
         updraftCBall_obs.append(updraftCBarrayobs)
         
         
     datetime_ICON_arr.append(new_dict[9]['datetime_iconlem'])
     

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
     LHSF_mod.append(new_dict[10]['LHF_iconlem'])
     SHSF_mod.append(new_dict[10]['SHF_iconlem'])
     LHSF_obs.append(new_dict[10]['LHF_obs'])
     SHSF_obs.append(new_dict[10]['SHF_obs'])
     LHSF_err_obs.append(new_dict[10]['LHF_Err_obs'])
     SHSF_err_obs.append(new_dict[10]['SHF_Err_obs'])     
     datetime_30m.append(new_dict[10]['datetime_30m'])
     LWF_mod.append(new_dict[10]['LW_iconlem'])
     SWF_mod.append(new_dict[10]['SW_iconlem'])
     LWF_obs.append(new_dict[10]['LW_obs'])
     SWF_obs.append(new_dict[10]['SW_obs'])
     LWF_err_obs.append(new_dict[10]['LW_Err_obs'])
     SWF_err_obs.append(new_dict[10]['SW_Err_obs'])
     
     pblHeightTW_mod.append(new_dict[9]['PBLHeightTW'])
     pblHeightRN_mod.append(new_dict[9]['PBLHeightRN'])
     rh_mod.append(new_dict[5]['relativeHumidity'])
     T_mod.append(new_dict[9]['T_iconlem'])
     
     cloudThickness_mod.append(new_dict[7]['cloudThickness'])
     cloudThickness_obs.append(new_dict[8]['cloudThickness'])
     cloudUpdraftCB_mod.append(new_dict[7]['cloudUpdraftCB'])
     cloudUpdraftCB_obs.append(new_dict[8]['cloudUpdraftCB'])

     # ----------------------------------------------------------------------------------------
     # ------- Analysis of the mean hourly profiles of variance of vertical velocity for obs. ICON-LEM, ICON-INSCAPE
     # ----------------------------------------------------------------------------------------
     #---- calculating mean variance and standard deviation profiles for each hour of the day for obs and model
     print('calculating mean variance and standard deviation profiles for each \
           hour of the day for obs and models')
     varianceWmean_obs = f_calcMeanStdVarProfiles(new_dict[3]['varianceW'], time[:], height[:],\
                                     date, yy, mm, dd, NprofilesOut, timeIncrement) 
     varianceWmean_mod = f_calcMeanStdVarProfiles(new_dict[9]['varianceW'], time[:], height[:], \
                                     date, yy, mm, dd, NprofilesOut, timeIncrement) 
 
     dataset_mean_variance_obs.append(varianceWmean_obs)
     dataset_mean_variance_mod.append(varianceWmean_mod)
     
     #print(LHSF_mod)
     #print(new_dict[11]['LHF_iconlem'])



#%%
# drop days with no clouds in model data (real indeces are 0, 7, 9) Since every loop one index is removed, 
# 
LengthFull = len(LWC_mod)
for ind in range(LengthFull):
    if (ind < len(LWC_mod)):
        print(LWC_mod[ind].shape)
        if (LWC_mod[ind].shape == (0,)):
            print('deleted')
            print(ind)
            del LWC_mod[ind]
            del duration_mod[ind]
            del meanheightFromCB_mod[ind]
            del cloudMeanCT_mod[ind]
            del cloudMeanCB_mod[ind]
            del cloudMaturity_mod[ind]
            del cloudUpdraftCB_mod[ind]
            del cloudThickness_mod[ind]


Ncloud_mod = len(np.concatenate(duration_mod, axis=0))
print('Number of clouds found in model')
print(len(np.concatenate(duration_mod, axis=0)))
print('Number of clouds found in obs')
print(len(np.concatenate(duration_obs, axis=0)))

# processing of LWC profile: goal is to derive hourly mean LWC profile from the collected clouds for model and obs
LWC_All_obs              = np.concatenate(LWC_obs, axis=0)
LWC_All_mod              = np.concatenate(LWC_mod, axis=0)
duration_all_obs         = np.concatenate(duration_obs, axis=0)
duration_all_mod         = np.concatenate(duration_mod, axis=0)
cloudMeanLWP_all_obs     = np.concatenate(LWC_obs, axis=0)
meanheightFromCB_All_obs = np.concatenate(meanheightFromCB_obs, axis=0)
meanheightFromCB_All_mod = np.concatenate(meanheightFromCB_mod, axis=0)
cloudMeanCT_All_obs      = np.concatenate(cloudMeanCT_obs, axis=0)
cloudMeanCB_All_obs      = np.concatenate(cloudMeanCB_obs, axis=0)
cloudMeanCT_All_mod      = np.concatenate(cloudMeanCT_mod, axis=0)
cloudMeanCB_All_mod      = np.concatenate(cloudMeanCB_mod, axis=0)
cloudThickness_All_obs   = np.concatenate(cloudThickness_obs, axis=0)
cloudThickness_All_mod   = np.concatenate(cloudThickness_mod, axis=0)
cloudUpdraftCB_All_obs   = np.concatenate(cloudUpdraftCB_obs, axis=0)
cloudUpdraftCB_All_mod   = np.concatenate(cloudUpdraftCB_mod, axis=0)
cloudLWP_All_mod         = np.concatenate(cloudLWP_mod, axis=0)
cloudLWP_All_obs         = np.concatenate(cloudLWP_obs, axis=0)
cloudMassFlux_All_mod    = np.concatenate(cloudMassFlux_mod, axis=0)
cloudMassFlux_All_obs    = np.concatenate(cloudMassFlux_obs, axis=0)


#duration_intervals = [0, 400, 800, 1600, 3200]
#duration_intervals       = [0, 200, 800, 1600]
duration_intervals       = [0, 900]
label_durations          = ['0-15 min','> 15 min']
LWC_obs_DF               = pd.DataFrame(LWC_All_obs, index=duration_all_obs, columns=height)
LWC_mod_DF               = pd.DataFrame(LWC_All_mod, index=duration_all_mod, columns=height)
meanHeightFromCB_obs_DF  = pd.DataFrame(meanheightFromCB_All_obs, index=duration_all_obs, columns=height)
meanHeightFromCB_mod_DF  = pd.DataFrame(meanheightFromCB_All_mod, index=duration_all_mod, columns=height)
CB_obs_DF                = pd.Series(cloudMeanCB_All_obs, index=duration_all_obs)
CB_mod_DF                = pd.Series(cloudMeanCB_All_mod, index=duration_all_mod)
CT_obs_DF                = pd.Series(cloudMeanCT_All_obs, index=duration_all_obs)
CT_mod_DF                = pd.Series(cloudMeanCT_All_mod, index=duration_all_mod)
cloudUpdraftCB_obs_DF    = pd.Series(cloudUpdraftCB_All_obs, index=duration_all_obs)
cloudUpdraftCB_mod_DF    = pd.Series(cloudUpdraftCB_All_mod, index=duration_all_mod)
cloudThickness_obs_DF    = pd.Series(cloudThickness_All_obs, index=duration_all_obs)
cloudThickness_mod_DF    = pd.Series(cloudThickness_All_mod, index=duration_all_mod)
cloudLWP_obs_DF          = pd.Series(cloudLWP_All_obs, index=duration_all_obs)
cloudLWP_mod_DF          = pd.Series(cloudLWP_All_mod, index=duration_all_mod)
cloudMassFlux_obs_DF     = pd.Series(cloudMassFlux_All_obs, index=duration_all_obs)
cloudMassFlux_mod_DF     = pd.Series(cloudMassFlux_All_mod, index=duration_all_mod)
#%%


# =============================================================================
# calculating and plotting IWV /LWP distributions and sigma(IWV)/sigma(LWP) 30 min std for all data
# =============================================================================
IWV_mod_nd = np.stack(IWV_mod).T
IWV_obs_nd = np.stack(IWV_obs).T
IWV_mod_nd = pd.DataFrame(IWV_mod_nd, index=datetime_ICON_arr[0], columns=np.arange(0,13))
IWV_obs_nd = pd.DataFrame(IWV_obs_nd, index=datetime_ICON_arr[0], columns=np.arange(0,13))

LWP_mod_nd = np.stack(LWP_all_mod).T
LWP_obs_nd = np.stack(LWP_all_obs).T
LWP_mod_nd = pd.DataFrame(LWP_mod_nd, index=datetime_ICON_arr[0], columns=np.arange(0,13))
LWP_obs_nd = pd.DataFrame(LWP_obs_nd, index=datetime_ICON_arr[0], columns=np.arange(0,13))

IWV_mod_nd.values[IWV_mod_nd.values < 10**(-7)] = np.nan
LWP_mod_nd.values[LWP_mod_nd.values < 10**(-7)] = np.nan
LWP_obs_nd.values[LWP_obs_nd.values < 0.] = np.nan
std_matrix_mod = IWV_mod_nd.resample('30min').std()   # calculating variance every 30 min for model
std_matrix_obs = IWV_obs_nd.resample('30min').std()   # calculating variance every 30 min for obs
std_LWP_matrix_mod = LWP_mod_nd.resample('30min').std()   # calculating variance every 30 min for model
std_LWP_matrix_obs = LWP_obs_nd.resample('30min').std()   # calculating variance every 30 min for obs

Nbins_IWV = 10
Nbins_LWP = 20
range_IWV = [5., 35.]
range_sigma_IWV = [0., 0.9]
range_sigma_LWP = [0., 0.3]
range_LWP = [0., 0.8]

fig, ax       = plt.subplots(nrows=2, ncols=2, figsize=(10,8))
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax                 = plt.subplot(2,2,1)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()  

# calculating totals
IWVDistr_mod =  IWV_mod_nd.values.flatten()
NtotIWV_mod = len(IWVDistr_mod[np.where(~np.isnan(IWVDistr_mod))])

IWVDistr_obs =  IWV_obs_nd.values.flatten()
NtotIWV_obs = len(IWVDistr_obs[np.where(~np.isnan(IWVDistr_obs))])

LWPDistr_mod =  LWP_mod_nd.values.flatten()
NtotLWP_mod = len(LWPDistr_mod[np.where(~np.isnan(LWPDistr_mod))])

LWPDistr_obs =  LWP_obs_nd.values.flatten()
NtotLWP_obs = len(LWPDistr_obs[np.where(~np.isnan(LWPDistr_obs))])

percentiles_LWPAll_mod = np.nanpercentile(LWPDistr_mod, [25, 50, 75, 90])
percentiles_LWPAll_obs = np.nanpercentile(LWPDistr_obs, [25, 50, 75, 90])
percentiles_IWVAll_mod = np.nanpercentile(IWVDistr_mod, [25, 50, 75, 90])
percentiles_IWVAll_obs = np.nanpercentile(IWVDistr_obs, [25, 50, 75, 90])

    
stringplot_IWV_obs     = 'median obs = '+str(round(percentiles_IWVAll_obs[1],3))+' $kg m^{-2}$'
stringplot_IWV_mod     = 'median mod = '+str(round(percentiles_IWVAll_mod[1],3))+' $kg m^{-2}$'
stringplot_LWP_obs     = 'median obs = '+str(round(percentiles_LWPAll_obs[1],3))+' $kg m^{-2}$'
stringplot_LWP_mod     = 'median mod = '+str(round(percentiles_LWPAll_mod[1],3))+' $kg m^{-2}$'


#ax.text(1800., ylabelArrCB[indPlot], stringplot_obs, fontsize=10)
#ax.text(1800., ylabelArrCB[indPlot]-0.0001, stringplot_mod, fontsize=10)
ax.set_ylim(0., 0.15)
plt.hist(IWV_obs_nd.values.flatten(), \
                      bins=Nbins_IWV, \
                      normed=True, \
                      color='black', \
                      range=range_IWV, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(IWV_obs_nd.values.flatten(), \
                      bins=Nbins_IWV, \
                      normed=True, \
                      color='black', \
                      range=range_IWV, \
                      histtype='step', \
                      cumulative=False)
plt.hist(IWV_mod_nd.values.flatten(), \
                      bins=Nbins_IWV, \
                      normed=True, \
                      color='red', \
                      range=range_IWV, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(IWV_mod_nd.values.flatten(), \
                      bins=Nbins_IWV, \
                      normed=True, \
                      color='red', \
                      range=range_IWV, \
                      histtype='step', \
                      cumulative=False)
ax.set_xlabel('IWV [Kg/m2]')
ax.set_ylabel('norm occ')
ax.text(21., 0.12, 'N model = '+str(NtotIWV_mod), fontsize=10)
ax.text(21., 0.11, 'N obs     = '+str(NtotIWV_obs), fontsize=10)
ax.text(21., 0.10, stringplot_IWV_obs, fontsize=10)
ax.text(21., 0.09, stringplot_IWV_mod, fontsize=10)

ax                 = plt.subplot(2,2,2)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()  
#ax.text(1800., ylabelArrCB[indPlot], stringplot_obs, fontsize=10)
#ax.text(1800., ylabelArrCB[indPlot]-0.0001, stringplot_mod, fontsize=10)
ax.set_ylim(0., 15.)
plt.hist(LWP_obs_nd.values.flatten(), \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='black', \
                      range=range_LWP, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(LWP_obs_nd.values.flatten(), \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='black', \
                      range=range_LWP, \
                      histtype='step', \
                      cumulative=False)
plt.hist(LWP_mod_nd.values.flatten(), \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='red', \
                      range=range_LWP, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(LWP_mod_nd.values.flatten(), \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='red', \
                      range=range_LWP, \
                      histtype='step', \
                      cumulative=False)
ax.set_xlabel('LWP [Kg/m2]')
ax.set_ylabel('norm occ')
ax.text(0.5, 12., 'N model = '+str(NtotLWP_mod), fontsize=10)
ax.text(0.5, 11., 'N obs     = '+str(NtotLWP_obs), fontsize=10)
ax.text(0.5, 10., stringplot_LWP_obs, fontsize=10)
ax.text(0.5, 9., stringplot_LWP_mod, fontsize=10)


ax                 = plt.subplot(2,2,3)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()  
#ax.text(1800., ylabelArrCB[indPlot], stringplot_obs, fontsize=10)
#ax.text(1800., ylabelArrCB[indPlot]-0.0001, stringplot_mod, fontsize=10)
ax.set_ylim(0., 5.)
Nbins_sigma_IWV = 10
plt.hist(std_matrix_obs.values.flatten(), \
                      bins=Nbins_sigma_IWV, \
                      normed=True, \
                      color='black', \
                      range=range_sigma_IWV, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(std_matrix_obs.values.flatten(), \
                      bins=Nbins_sigma_IWV, \
                      normed=True, \
                      color='black', \
                      range=range_sigma_IWV, \
                      histtype='step', \
                      cumulative=False)
plt.hist(std_matrix_mod.values.flatten(), \
                      bins=Nbins_sigma_IWV, \
                      normed=True, \
                      color='red', \
                      range=range_sigma_IWV, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(std_matrix_mod.values.flatten(), \
                      bins=Nbins_sigma_IWV, \
                      normed=True, \
                      color='red', \
                      range=range_sigma_IWV, \
                      histtype='step', \
                      cumulative=False)
ax.set_xlabel('${\sigma}$(IWV) ')
ax.set_ylabel('norm occ')

ax                 = plt.subplot(2,2,4)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()  
#ax.text(1800., ylabelArrCB[indPlot], stringplot_obs, fontsize=10)
#ax.text(1800., ylabelArrCB[indPlot]-0.0001, stringplot_mod, fontsize=10)
ax.set_ylim(0., 50.)
Nbins_sigma_LWP = 20
plt.hist(std_LWP_matrix_obs.values.flatten(), \
                      bins=Nbins_sigma_LWP, \
                      normed=True, \
                      color='black', \
                      range=range_sigma_LWP, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(std_LWP_matrix_obs.values.flatten(), \
                      bins=Nbins_sigma_LWP, \
                      normed=True, \
                      color='black', \
                      range=range_sigma_LWP, \
                      histtype='step', \
                      cumulative=False)
plt.hist(std_LWP_matrix_mod.values.flatten(), \
                      bins=Nbins_sigma_LWP, \
                      normed=True, \
                      color='red', \
                      range=range_sigma_LWP, \
                      cumulative=False, \
                      alpha=0.2)    
plt.hist(std_LWP_matrix_mod.values.flatten(), \
                      bins=Nbins_sigma_LWP, \
                      normed=True, \
                      color='red', \
                      range=range_sigma_LWP, \
                      histtype='step', \
                      cumulative=False)
ax.set_xlabel('${\sigma}$(LWP) ')
ax.set_ylabel('norm occ')
fig.tight_layout()

fig.savefig(pathFig+'LWP_IWV_alldata.png', format='png') 
#percentiles_mod = np.nanpercentile(std_matrix_mod.values, [25, 50, 75, 90])
#percentiles_obs = np.nanpercentile(std_matrix_obs.values, [25, 50, 75, 90])
#percentiles_LWP_mod = np.nanpercentile(std_matrix_mod.values, [25, 50, 75, 90])
#percentiles_LWP_obs = np.nanpercentile(std_matrix_obs.values, [25, 50, 75, 90])
#xlabelArr = ['${\sigma}$(IWV) [Kg/m^2]', '${\sigma}$(LWP) [Kg/m^2]']
#variablesObsArr = [std_matrix_obs.values.flatten(), std_LWP_matrix_obs.values.flatten()]
#variablesModArr = [std_matrix_mod.values.flatten(), std_LWP_matrix_mod.values.flatten()]
#titleArr = [' Distributions of ${\sigma}$(IWV) over 30 min' , ' Distributions of ${\sigma}$(LWP) over 30 min' ]
#perc_mod_arr = [percentiles_mod[1], percentiles_LWP_mod[1]]
#perc_obs_arr = [percentiles_obs[1], percentiles_LWP_obs[1]]

#%%
"""
This part of the code has the goal to characterize the macroscopic cloud properties
of clouds with duration smaller than 15 min (short lived) 
and clouds with duration larger than 15 min (long lived)
Here we: 
1_ calculate mean profiles of LWC for short and long lived clouds : 
    concept is to select all profiles within the specific duration interval. 
    Then, we select all heights between 0 and 1 and we loop on them: for every height
2_ calculate distributions of LWP,geometrical thickness, maturity, updrafts speeds at cloud base
    for short and long lived clouds
"""
height_grid                    = np.linspace(0., 1., num=10)
Nprofiles_durationIntervalsObs = []
Nprofiles_durationIntervalsMod = []

LWC_mean_rescaled_obs   = []
LWC_mean_rescaled_mod   = []
CB_obs_distr_arr        = []
CB_mod_distr_arr        = []
CT_obs_distr_arr        = []
CT_mod_distr_arr        = []
thickness_obs_distr_arr = []
thickness_mod_distr_arr = []
updraftCB_obs_distr_arr = []
updraftCB_mod_distr_arr = []
LWP_obs_distr_arr       = []
LWP_mod_distr_arr       = []
massFlux_obs_distr_arr  = []
massFlux_mod_distr_arr  = []

for indDuration in range(0, len(duration_intervals)):

    # case for duration > last element of duration_intervals array
    if ((indDuration-len(duration_intervals)+1) == 0.):
        mask_duration_obs    = LWC_obs_DF.index > duration_intervals[indDuration]
        mask_duration_mod    = LWC_mod_DF.index > duration_intervals[indDuration]
        #print('interval sampled > '+str(duration_intervals[indDuration])+' s')
    else:
        # case for all intermediate duration intervals
        mask_duration_obs    = (LWC_obs_DF.index > duration_intervals[indDuration]) \
        * (LWC_obs_DF.index <= duration_intervals[indDuration+1])
        mask_duration_mod    = (LWC_mod_DF.index > duration_intervals[indDuration]) \
        * (LWC_mod_DF.index <= duration_intervals[indDuration+1])
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
    thickness_obs_distr_arr.append(cloudThickness_obs_DF.loc[mask_duration_obs])
    thickness_mod_distr_arr.append(cloudThickness_mod_DF.loc[mask_duration_mod])
    updraftCB_obs_distr_arr.append(cloudUpdraftCB_obs_DF.loc[mask_duration_obs])
    updraftCB_mod_distr_arr.append(cloudUpdraftCB_mod_DF.loc[mask_duration_mod])
    LWP_obs_distr_arr.append(cloudLWP_obs_DF.loc[mask_duration_obs])
    LWP_mod_distr_arr.append(cloudLWP_mod_DF.loc[mask_duration_mod])
    massFlux_obs_distr_arr.append(cloudMassFlux_obs_DF.loc[mask_duration_obs])
    massFlux_mod_distr_arr.append(cloudMassFlux_mod_DF.loc[mask_duration_mod])    


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

# plotting cloud LWP, updrafts at cloud base and mass flux distributions for short and long lived clouds
fig = plt.figure(figsize=(12,8))
plt.gcf().subplots_adjust(bottom=0.15)
stringmidHours = ['0-15 min', '>15 min']

Rangearr_LWP = [[0., 0.7], [0., 0.7]]#,[600., 3000.], [600., 3000.]]
Rangearr_updraft_CB = [[-3., 3], [-3., 3]]#,[1300.,5000.], [1300., 5000.]]
Rangearr_MassFlux = [[-4000.,4000.], [-4000.,4000.]]
Nbins_LWP = 15 
Nbins_updrafts = 15
Nbins_massFlux = 15
ylabelArrLWP = [15., 15.]
ylabelArrUpdraft = [1.5, 1.5]
ylabelArrMassFLux = [1.5, 1.5]

indLWP = [1,4]
indUpdraft  = [2,5]
indMassflux = [3,6]
ymaxArrLWP = [20., 20.]
ymaxArrUpdraft = [2., 2.]
ymaxArrMassFlux = [2., 2.]


# loop on duration intervals (2)
for indPlot in range(0,2):
    print(indPlot)
    LWP_cloud_obs_distr = LWP_obs_distr_arr[indPlot]
    LWP_cloud_obs_distr[LWP_cloud_obs_distr < 0.] = np.nan
    LWP_cloud_mod_distr = LWP_mod_distr_arr[indPlot]
    LWP_cloud_mod_distr[LWP_cloud_mod_distr < 0.] = np.nan
    
    percentiles_LWP_mod = np.nanpercentile(LWP_cloud_mod_distr, [25, 50, 75, 90])
    percentiles_LWP_obs = np.nanpercentile(LWP_cloud_obs_distr, [25, 50, 75, 90])
    percentiles_updraft_mod = np.nanpercentile(updraftCB_mod_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_updraft_obs = np.nanpercentile(updraftCB_obs_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_massFlux_mod = np.nanpercentile(massFlux_mod_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_massFlux_obs = np.nanpercentile(massFlux_obs_distr_arr[indPlot].values, [25, 50, 75, 90])
    print(percentiles_LWP_mod)
    print(percentiles_LWP_obs)
    print(percentiles_updraft_mod)
    print(percentiles_updraft_obs)
    print(percentiles_massFlux_mod)
    print(percentiles_massFlux_obs)
    print('*********************')
    

    stringplot_obs     = 'median obs = '+str(round(percentiles_LWP_obs[1],3))+' $kg m^{-2}$'
    stringplot_mod     = 'median mod = '+str(round(percentiles_LWP_mod[1],3))+' $kg m^{-2}$'
    ax                 = plt.subplot(2,3,indLWP[indPlot])  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    ax.text(0.1, ylabelArrLWP[indPlot], stringplot_obs, fontsize=10)
    ax.text(0.1, ylabelArrLWP[indPlot]-1., stringplot_mod, fontsize=10)
    ax.set_ylabel('normalized occurrences')
    ax.set_xlabel('LWP [$kg m^{-2}$]') 
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')
    ax.set_ylim(0., ymaxArrLWP[indPlot])
    plt.hist(LWP_cloud_obs_distr, \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_LWP[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)    
    plt.hist(LWP_cloud_obs_distr, \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_LWP[indPlot-1], \
                      histtype='step', \
                      cumulative=False)
    plt.hist(LWP_cloud_mod_distr, \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_LWP[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(LWP_cloud_mod_distr, \
                      bins=Nbins_LWP, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_LWP[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    
    
    ax                 = plt.subplot(2,3,indUpdraft[indPlot])  
    ax.set_title(stringmidHours[indPlot])
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    stringplot_obs = 'median obs = '+str(round(percentiles_updraft_obs[1],3))+' $ ms^{-1}$'
    #plt.text(0.8, ymax-2.*ymax/10., 'median mod = '+str(round(percenties_mod[1], 2)))
    stringplot_mod = 'median mod = '+str(round(percentiles_updraft_mod[1],3))+' $ ms^{-1}$'
    ax.legend(loc='upper right', fontsize=12, frameon=False)
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('normalized occurrences')
    ax.set_xlabel('Updraft at cloud base [$ ms^{-1}$]')   
    ax.text(0., ylabelArrUpdraft[indPlot], stringplot_obs, fontsize=10)
    ax.text(0., ylabelArrUpdraft[indPlot]-0.1, stringplot_mod, fontsize=10)
    ax.set_ylim(0., ymaxArrUpdraft[indPlot])
    plt.hist(updraftCB_obs_distr_arr[indPlot], \
                      bins=Nbins_updrafts, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_updraft_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(updraftCB_obs_distr_arr[indPlot], \
                      bins=Nbins_updrafts, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_updraft_CB[indPlot-1], \
                      cumulative=False, \
                      histtype='step') 
    plt.hist(updraftCB_mod_distr_arr[indPlot], \
                      bins=Nbins_updrafts, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_updraft_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(updraftCB_mod_distr_arr[indPlot], \
                      bins=Nbins_updrafts, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_updraft_CB[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')


    stringplot_obs     = 'median obs = '+str(round(percentiles_massFlux_obs[1],3))+' $kg ms^{-1}$'
    stringplot_mod     = 'median mod = '+str(round(percentiles_massFlux_mod[1],3))+' $kg ms^{-1}$'
    
    
    ax                 = plt.subplot(2,3,indMassflux[indPlot])  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    ax.text(0.1, ylabelArrMassFLux[indPlot], stringplot_obs, fontsize=10)
    ax.text(0.1, ylabelArrMassFLux[indPlot]-0.5, stringplot_mod, fontsize=10)
    ax.set_ylabel('normalized occurrences')
    ax.set_xlabel('Mass Flux [$kg (ms)^{-1}$]') 
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')
    ax.set_ylim(0., ymaxArrMassFlux[indPlot])
    plt.hist(massFlux_obs_distr_arr[indPlot], \
                      bins=Nbins_massFlux, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_MassFlux[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)    
    plt.hist(massFlux_obs_distr_arr[indPlot], \
                      bins=Nbins_massFlux, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_MassFlux[indPlot-1], \
                      histtype='step', \
                      cumulative=False)
    plt.hist(massFlux_mod_distr_arr[indPlot], \
                      bins=Nbins_massFlux, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_MassFlux[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(massFlux_mod_distr_arr[indPlot], \
                      bins=Nbins_massFlux, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_MassFlux[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    

plt.tight_layout()
plt.savefig(pathFig+'LWP_updrafts_massFlux_distributions_vs_duration.png', format='png') 






#%%   

# plotting CB , CT, geom thickness distributions for short and long lived clouds
fig = plt.figure(figsize=(12,8))
plt.gcf().subplots_adjust(bottom=0.15)
stringmidHours = ['0-15 min', '>15 min']

Rangearr_CB = [[600., 3000.], [600., 3000.]]#,[600., 3000.], [600., 3000.]]
Rangearr_CT = [[800., 5000.], [800., 5000.]]#,[1300.,5000.], [1300., 5000.]]
Rangearr_Thickness =  [[0., 1500.], [0., 2000.]]#, [0., 600.], [0., 600.]]
Nbins_CB = 10
Nbins_CT = 10
Nbins_thickness = 10
ymaxArr = [0.0017, 0.0017, 0.03]#, 0.025]
#xlabelArr = [1400., 5000., 1000., 250.]
ylabelArrCB = [0.0013, 0.0013]
ylabelArrCT = [0.0013, 0.0013]
ylabelArrCloudThickness = [0.003, 0.003]

indCloudBase = [1,4]
indCloudTop  = [2,5]
indThickness = [3,6]
ymaxArrCB = [0.0017, 0.0017]
ymaxArrCT = [0.0017, 0.0017]
ymaxArrThickness = [0.004, 0.004]



# loop on duration intervals (2)
for indPlot in range(0,2):
    print(indPlot)
    percentiles_CB_mod = np.nanpercentile(CB_mod_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_CB_obs = np.nanpercentile(CB_obs_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_CT_mod = np.nanpercentile(CT_mod_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_CT_obs = np.nanpercentile(CT_obs_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_thickness_mod = np.nanpercentile(thickness_mod_distr_arr[indPlot].values, [25, 50, 75, 90])
    percentiles_thickness_obs = np.nanpercentile(thickness_obs_distr_arr[indPlot].values, [25, 50, 75, 90])
    print(percentiles_CB_mod)
    print(percentiles_CB_obs)
    print(percentiles_CT_mod)
    print(percentiles_CT_obs)
    print(percentiles_thickness_mod)
    print(percentiles_thickness_obs)
    print('*********************')
    
    stringplot_obs     = 'median obs = '+str(round(percentiles_CB_obs[1],1))+' m'
    stringplot_mod     = 'median mod = '+str(round(percentiles_CB_mod[1],1))+' m'
    ax                 = plt.subplot(2,3,indCloudBase[indPlot])  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
    matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots

    ax.text(1800., ylabelArrCB[indPlot], stringplot_obs, fontsize=10)
    ax.text(1800., ylabelArrCB[indPlot]-0.0001, stringplot_mod, fontsize=10)
    ax.set_ylabel('occurrences')
    ax.set_xlabel('Cloud base height [m]') 
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')
    ax.set_ylim(0., ymaxArrCB[indPlot])
    plt.hist(CB_obs_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)    
    plt.hist(CB_obs_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CB[indPlot-1], \
                      histtype='step', \
                      cumulative=False)
    plt.hist(CB_mod_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(CB_mod_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    ax                 = plt.subplot(2,3,indCloudTop[indPlot])  
    ax.set_title(stringmidHours[indPlot])
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    stringplot_obs = 'median obs = '+str(round(percentiles_CT_obs[1],1))+' m'
    #plt.text(0.8, ymax-2.*ymax/10., 'median mod = '+str(round(percenties_mod[1], 2)))
    stringplot_mod = 'median mod = '+str(round(percentiles_CT_mod[1],1))+' m'
    ax.legend(loc='upper right', fontsize=12, frameon=False)
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('Cloud top height [m]')   
    ax.text(1200., ylabelArrCT[indPlot], stringplot_obs, fontsize=10)
    ax.text(1200., ylabelArrCT[indPlot]-0.0001, stringplot_mod, fontsize=10)
    ax.set_ylim(0., ymaxArrCT[indPlot])
    plt.hist(CT_obs_distr_arr[indPlot], \
                      bins=Nbins_CT, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(CT_obs_distr_arr[indPlot], \
                      bins=Nbins_CT, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      histtype='step') 
    plt.hist(CT_mod_distr_arr[indPlot], \
                      bins=Nbins_CT, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(CT_mod_distr_arr[indPlot], \
                      bins=Nbins_CT, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    ax.xaxis.set_ticks(np.arange(500.,5500, 1000))
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')


    ax                 = plt.subplot(2,3,indThickness[indPlot])  
    ax.set_title(stringmidHours[indPlot])
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
    stringplot_obs = 'median obs = '+str(round(percentiles_thickness_obs[1],1))+' m'
    #plt.text(0.8, ymax-2.*ymax/10., 'median mod = '+str(round(percenties_mod[1], 2)))
    stringplot_mod = 'median mod = '+str(round(percentiles_thickness_mod[1],1))+' m'
    ax.legend(loc='upper right', fontsize=12, frameon=False)
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('Geometrical thickness [m]')   
    ax.text(500., ylabelArrCloudThickness[indPlot], stringplot_obs, fontsize=10)
    ax.text(500., ylabelArrCloudThickness[indPlot]-0.0002, stringplot_mod, fontsize=10)
    ax.set_ylim(0., ymaxArrThickness[indPlot])
    plt.hist(thickness_obs_distr_arr[indPlot], \
                      bins=Nbins_thickness, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_Thickness[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(thickness_obs_distr_arr[indPlot], \
                      bins=Nbins_thickness, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_Thickness[indPlot-1], \
                      cumulative=False, \
                      histtype='step') 
    plt.hist(thickness_mod_distr_arr[indPlot], \
                      bins=Nbins_thickness, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_Thickness[indPlot-1], \
                      cumulative=False, \
                      alpha=0.2)
    plt.hist(thickness_mod_distr_arr[indPlot], \
                      bins=Nbins_thickness, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_Thickness[indPlot-1], \
                      cumulative=False, \
                      histtype='step')
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')
    ax.xaxis.set_ticks(np.arange(0.,2000., 500))

plt.tight_layout()
plt.savefig(pathFig+'CB_CT_Thickness_distributions_vs_duration.png', format='png') 

#%%
# plotting all profiles for each duration interval after regridding
fig, ax       = plt.subplots(nrows=2, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.005, 0.003]#, 0.004, 0.01]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr_obs = [1,2]#,3,4]
indPlotArr_mod = [3,4]#5,6,7,8]

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
plt.savefig(pathFig+'LWC_all_global_profiles.png', format='png') 

    
#%%
fig, ax       = plt.subplots(nrows=1, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 0.88888889#1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.002, 0.002]#, 0.001, 0.001]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr = [1,2]#,3,4]
for ind in range(0, len(LWC_mean_rescaled_obs)):
     xmax          = xmaxArr[ind]
     ax        = plt.subplot(1,len(duration_intervals),indPlotArr[ind])  
     ax.spines["top"].set_visible(False)  
     ax.spines["right"].set_visible(False)  
     ax.get_xaxis().tick_bottom()  
     ax.get_yaxis().tick_left() 
     ax.set_ylim(ymin, ymax)
     ax.set_xlim(xmin, xmax)
     ax.set_xlabel('LWC [$kg m^{-3}$]')
     ax.set_ylabel('rel. dist. from cloud base')
     ax.set_title(stringmidHours[ind])
     plt.plot(LWC_mean_prof_obs[ind,:], height_grid, color='black', label='obs')

     y1        = LWC_mean_prof_obs[ind,:]-LWC_std_prof_obs[ind,:]
     y2        = LWC_mean_prof_obs[ind,:]+LWC_std_prof_obs[ind,:]
     print(y1)
     print(y2)
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
     
     plt.plot(LWC_mean_prof_mod[ind,:], height_grid, color='red', label='icon-lem')
     y1        = LWC_mean_prof_mod[ind,:]-LWC_std_prof_mod[ind,:]
     y2        = LWC_mean_prof_mod[ind,:]+LWC_std_prof_mod[ind,:]
     print(y1)
     print(y2)
     plt.legend(loc='upper right', fontsize=12, frameon=False)
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
plt.tight_layout()

plt.savefig(pathFig+'LWC_mean_global_profiles.png', format='png') 


#%%
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
                                                              lts_radiosondes, \
                                                              pblHeight_radiosondes, \
                                                              theta_v_mod, \
                                                              T_mod, \
                                                              rh_mod, \
                                                              time_mod, \
                                                              height_mod, \
                                                              lcl_mod, \
                                                              lts_mod, \
                                                              pblHeightRN_mod)


from myFunctions import f_calculateMeanProfilesPlotThetaVRadiosondes
result = f_calculateMeanProfilesPlotThetaVRadiosondes(theta_v_dict_obs_mod_arr, height_mod)

def f_calculateMeanThetaVModelProfiles(time_radiosondes, \
                                       theta_v_radiosondes,\
                                       T_radiosondes, \
                                       rh_radiosObs, \
                                       height_radiosondes, \
                                       lcl_radiosondes, \
                                       lts_radiosondes, \
                                       pblHeight_radiosondes, \
                                       theta_v_mod, \
                                       T_mod, \
                                       rh_mod, \
                                       time_mod, \
                                       height_mod, \
                                       lcl_mod, \
                                       lts_mod, \
                                       pblHeight_mod):
MatrixHourMeanProfileThetaRad = result[0]
MatrixHourStdProfileThetaRad  = result[1]
listHourDict                  = result[2]
MatrixHourMeanProfileTRad     = result[3]
MatrixHourStdProfileTRad      = result[4]
MatrixHourMeanProfileRHRad    = result[5]
MatrixHourStdProfileRHRad     = result[6]

gridHeight                    = height_mod[0]

# hours with radiosondes: 5,7,8,9,11,13,15,17,19,20,21,23
# indeces               : 0,1,2,3, 4, 5, 6, 7, 8, 9,10,11
#%%
indexPlot = [1,3,4,5,6,7,11]

if flagPlotThetaVglobalProfiles == 1:
    Ncols = 4
    Nrows = 2
    Nplots = 11
    
    
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = height_mod[0][-1]
    xmin          = 283.
    xmax          = 298.
    fontSizeTitle = 16
    fontSizeX     = 12
    fontSizeY     = 12
    #timeTitles = [']
    indxArr = [0,0,0,0,1,1,1]
    indyArr = [0,1,2,3,0,1,2]
    lclx_obs_Arr = [289.2, 289.8, 290.8, 292., 292., 292.2, 288.]
    lclx_mod_Arr = [286., 287., 288.1, 290., 291.2, 293., 289.5]

    for indPlot in range(0,len(indxArr)):
        
        # reading number of profiles, lcl and pbl heights for the selected hour
        Nprofiles = listHourDict[indexPlot[indPlot]]['Nprofiles']
        lcl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_rad_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_mod_hour'])
        lcl_mod_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])

        pbl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_rad_hour'])
        pbl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_mod_hour'])
        
        # assigning indeces for subplot positions 
        indx      = indxArr[indPlot]
        indy      = indyArr[indPlot]
        
        #removing subplot box top and right lines
        ax[indx,indy].spines["top"].set_visible(False)  
        ax[indx,indy].spines["right"].set_visible(False)  
        ax[indx,indy].get_xaxis().tick_bottom()  
        ax[indx,indy].get_yaxis().tick_left() 
        ax[indx,indy].text(284, 2000., 'N = '+str(Nprofiles), fontsize=10)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        prof_mod  = listHourDict[indexPlot[indPlot]]['meanProfile_mod']
        prof_obs  = MatrixHourMeanProfileThetaRad[:, indexPlot[indPlot]]
        std_mod   = listHourDict[indexPlot[indPlot]]['stdProfileMod']
        labelHour = listHourDict[indexPlot[indPlot]]['hour']
        std_obs   = MatrixHourStdProfileThetaRad[:, indexPlot[indPlot]]
        ax[indx,indy].plot(prof_obs, gridHeight, label='obs '+str(labelHour)+' UTC',  color='black')
        ax[indx,indy].plot(prof_mod, height_mod[0], label='icon-lem',  color='red')
        y1        = prof_obs-std_obs
        y2        = prof_obs+std_obs
        ax[indx,indy].fill_betweenx(gridHeight, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = prof_mod-std_mod
        y2        = prof_mod+std_mod
        ax[indx,indy].fill_betweenx(height_mod[0], y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        ax[indx,indy].errorbar(lclx_obs_Arr[indPlot], lcl_rad_plot, fmt='ko', xerr=0, yerr=lcl_rad_err, ecolor='black')
        ax[indx,indy].errorbar(lclx_mod_Arr[indPlot], lcl_mod_plot, fmt='ro', xerr=0, yerr=lcl_mod_err, ecolor='red')
        ax[indx,indy].legend(loc='upper left', fontsize=12, frameon=False)
        ax[indx,indy].set_ylim(ymin,ymax)
        ax[indx,indy].set_xlim(xmin,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        ax[indx,indy].set_xlabel('${\Theta_v}$[K]', fontsize=fontSizeX)
        ax[indx,indy].set_ylabel('height [m]', fontsize=fontSizeY)
    fig.subplots_adjust(hspace=0.15, bottom=0.1,left=0.05)   
    ax[1,3].set_visible(False) # to remove last plot
    fig.savefig(pathFig+'theta_v_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')






    
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = height_mod[0][-1]
    xmin          = -5.
    xmax          = 15.
    fontSizeTitle = 16
    fontSizeX     = 12
    fontSizeY     = 12
    #timeTitles = [']
    indxArr = [0,0,0,0,1,1,1]
    indyArr = [0,1,2,3,0,1,2]

    for indPlot in range(0,len(indxArr)):
        
        # reading number of profiles, lcl and pbl heights for the selected hour
        Nprofiles = listHourDict[indexPlot[indPlot]]['Nprofiles']
        
        # assigning indeces for subplot positions 
        indx      = indxArr[indPlot]
        indy      = indyArr[indPlot]
        
        #removing subplot box top and right lines
        ax[indx,indy].spines["top"].set_visible(False)  
        ax[indx,indy].spines["right"].set_visible(False)  
        ax[indx,indy].get_xaxis().tick_bottom()  
        ax[indx,indy].get_yaxis().tick_left() 
        ax[indx,indy].text(284, 2000., 'N = '+str(Nprofiles), fontsize=10)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        prof_mod  = listHourDict[indexPlot[indPlot]]['meanProfile_mod']
        prof_obs  = MatrixHourMeanProfileThetaRad[:, indexPlot[indPlot]]
        std_mod   = listHourDict[indexPlot[indPlot]]['stdProfileMod']
        labelHour = listHourDict[indexPlot[indPlot]]['hour']
        std_obs   = MatrixHourStdProfileThetaRad[:, indexPlot[indPlot]]
        
        # reading values at 2500 for obs and mod
        thetaV2500    = prof_obs[f_closest(gridHeight, 2500.)]
        thetaV2500Arr = np.repeat(thetaV2500, len(gridHeight))
        thetaVdiffObs = thetaV2500Arr - prof_obs
        thetaVdiffMod = thetaV2500Arr - prof_mod
        ax[indx,indy].plot(thetaVdiffObs, gridHeight, label='obs '+str(labelHour)+' UTC',  color='black')
        ax[indx,indy].plot(thetaVdiffMod, gridHeight, label='icon-lem',  color='red')
        y1        = thetaV2500Arr-prof_obs-std_obs
        y2        = thetaV2500Arr-prof_obs+std_obs
        ax[indx,indy].fill_betweenx(gridHeight, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = thetaV2500Arr-prof_mod-std_mod
        y2        = thetaV2500Arr-prof_mod+std_mod
        ax[indx,indy].fill_betweenx(height_mod[0], y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        ax[indx,indy].legend(loc='upper right', fontsize=12, frameon=False)
        ax[indx,indy].set_ylim(ymin,ymax)
        ax[indx,indy].set_xlim(xmin,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        ax[indx,indy].set_xlabel('${\Delta\Theta_v}$[K]', fontsize=fontSizeX)
        ax[indx,indy].set_ylabel('height [m]', fontsize=fontSizeY)
    fig.subplots_adjust(hspace=0.15, bottom=0.1, left=0.05)   
    ax[1,3].set_visible(False) # to remove last plot
    fig.savefig(pathFig+'theta_v_globalMeanDataset_difference_diurnal_cycle_obs_mod.png', format='png')


#%%
#### PLOT for Temperature profiles
flagPlotTglobalProfiles = 1
if flagPlotTglobalProfiles == 1:
    Ncols = 4
    Nrows = 2
    Nplots = 11
    
    
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(12,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = height_mod[0][-1]
    xmin          = 260.
    xmax          = 295.
    fontSizeTitle = 16
    fontSizeX     = 12
    fontSizeY     = 12
    #timeTitles = [']
    indxArr = [0,0,0,0,1,1,1]
    indyArr = [0,1,2,3,0,1,2]
    lclx_obs_Arr = [289.2, 289.8, 290.8, 292., 292., 292.2, 288.]
    lclx_mod_Arr = [286., 287., 288.1, 290., 291.2, 293., 289.5]

    for indPlot in range(0,len(indxArr)):
        
        # reading number of profiles, lcl and pbl heights for the selected hour
        Nprofiles = listHourDict[indexPlot[indPlot]]['Nprofiles']
        lcl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_rad_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_mod_hour'])
        lcl_mod_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])

        pbl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_rad_hour'])
        pbl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_mod_hour'])
        
        # assigning indeces for subplot positions 
        indx      = indxArr[indPlot]
        indy      = indyArr[indPlot]
        
        #removing subplot box top and right lines
        ax[indx,indy].spines["top"].set_visible(False)  
        ax[indx,indy].spines["right"].set_visible(False)  
        ax[indx,indy].get_xaxis().tick_bottom()  
        ax[indx,indy].get_yaxis().tick_left() 
        ax[indx,indy].text(284, 2000., 'N = '+str(Nprofiles), fontsize=10)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        prof_mod  = listHourDict[indexPlot[indPlot]]['meanProfile_T_mod']
        prof_obs  = MatrixHourMeanProfileTRad[:, indexPlot[indPlot]]
        std_mod   = listHourDict[indexPlot[indPlot]]['stdProfile_T_Mod']
        labelHour = listHourDict[indexPlot[indPlot]]['hour']
        std_obs   = MatrixHourStdProfileTRad[:, indexPlot[indPlot]]
        ax[indx,indy].plot(prof_obs, gridHeight, label='obs '+str(labelHour)+' UTC',  color='black')
        ax[indx,indy].plot(prof_mod, height_mod[0], label='icon-lem',  color='red')
        y1        = prof_obs-std_obs
        y2        = prof_obs+std_obs
        ax[indx,indy].fill_betweenx(gridHeight, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = prof_mod-std_mod
        y2        = prof_mod+std_mod
        ax[indx,indy].fill_betweenx(height_mod[0], y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        ax[indx,indy].legend(loc='upper left', fontsize=12, frameon=False)
        ax[indx,indy].set_ylim(ymin,ymax)
        ax[indx,indy].set_xlim(xmin,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        ax[indx,indy].set_xlabel('${T}$[K]', fontsize=fontSizeX)
        ax[indx,indy].set_ylabel('height [m]', fontsize=fontSizeY)
    fig.subplots_adjust(hspace=0.15, bottom=0.1,)   
    ax[1,3].set_visible(False) # to remove last plot
    fig.savefig(pathFig+'T_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')


#%%

# Plot for RH profiles 
flagPlotRHglobalProfiles = 1
if flagPlotRHglobalProfiles == 1:
    Ncols = 4
    Nrows = 2    
    
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(12,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = height_mod[0][-1]
    xmin          = 0.
    xmax          = 100.
    fontSizeTitle = 16
    fontSizeX     = 12
    fontSizeY     = 12
    #timeTitles = [']
    indxArr = [0,0,0,0,1,1,1]
    indyArr = [0,1,2,3,0,1,2]
    lclx_obs_Arr = [289.2, 289.8, 290.8, 292., 292., 292.2, 288.]
    lclx_mod_Arr = [286., 287., 288.1, 290., 291.2, 293., 289.5]

    for indPlot in range(0,len(indxArr)):
        
        # reading number of profiles, lcl and pbl heights for the selected hour
        Nprofiles = listHourDict[indexPlot[indPlot]]['Nprofiles']
        lcl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_rad_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])
        lcl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['lcl_mod_hour'])
        lcl_mod_err    = np.nanstd(listHourDict[indexPlot[indPlot]]['lcl_rad_hour'])

        pbl_rad_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_rad_hour'])
        pbl_mod_plot   = np.nanmedian(listHourDict[indexPlot[indPlot]]['pblHeight_mod_hour'])
        
        # assigning indeces for subplot positions 
        indx      = indxArr[indPlot]
        indy      = indyArr[indPlot]
        
        #removing subplot box top and right lines
        ax[indx,indy].spines["top"].set_visible(False)  
        ax[indx,indy].spines["right"].set_visible(False)  
        ax[indx,indy].get_xaxis().tick_bottom()  
        ax[indx,indy].get_yaxis().tick_left() 
        ax[indx,indy].text(284, 2000., 'N = '+str(Nprofiles), fontsize=10)
        matplotlib.rc('xtick', labelsize=10)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)                        # sets dimension of ticks in the plots
        prof_mod  = listHourDict[indexPlot[indPlot]]['meanProfile_RH_mod']
        prof_obs  = MatrixHourMeanProfileRHRad[:, indexPlot[indPlot]]
        std_mod   = listHourDict[indexPlot[indPlot]]['stdProfile_RH_Mod']
        labelHour = listHourDict[indexPlot[indPlot]]['hour']
        std_obs   = MatrixHourStdProfileRHRad[:, indexPlot[indPlot]]
        ax[indx,indy].plot(prof_obs, gridHeight, label='obs '+str(labelHour)+' UTC',  color='black')
        ax[indx,indy].plot(prof_mod, height_mod[0], label='icon-lem',  color='red')
        y1        = prof_obs-std_obs
        y2        = prof_obs+std_obs
        ax[indx,indy].fill_betweenx(gridHeight, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = prof_mod-std_mod
        y2        = prof_mod+std_mod
        ax[indx,indy].fill_betweenx(height_mod[0], y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        ax[indx,indy].errorbar(lclx_obs_Arr[indPlot], lcl_rad_plot, fmt='ko', xerr=0, yerr=lcl_rad_err, ecolor='black')
        ax[indx,indy].errorbar(lclx_mod_Arr[indPlot], lcl_mod_plot, fmt='ro', xerr=0, yerr=lcl_mod_err, ecolor='red')
        ax[indx,indy].legend(loc='upper left', fontsize=12, frameon=False)
        ax[indx,indy].set_ylim(ymin,ymax)
        ax[indx,indy].set_xlim(xmin,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        ax[indx,indy].set_xlabel('${RH}$[%]', fontsize=fontSizeX)
        ax[indx,indy].set_ylabel('height [m]', fontsize=fontSizeY)
    fig.subplots_adjust(hspace=0.15, bottom=0.1,)   
    ax[1,3].set_visible(False) # to remove last plot
    fig.savefig(pathFig+'RH_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')


#%%
    

#%%

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

# =============================================================================
# 


#%% 
# std(IWV) distribution per intervals of the day
# =============================================================================
# hours = [0, 6, 9, 12, 15, 18, 23]
# 
# nbins   = 20
# ymax    = 10.
# xmin    = 0.
# xmax    = 1.
# indplot = 1
# fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
# matplotlib.rcParams['savefig.dpi'] = 100
# plt.gcf().subplots_adjust(bottom=0.15)
# 
# for indHour in range(len(hours)-1):
#     hourInf = hours[indHour]
#     hourSup = hours[indHour+1]
#     hinf_idx = datetime.datetime(2013,4,24,hourInf,0,0)
#     hsup_idx = datetime.datetime(2013,4,24,hourSup,0,0)
#     maskT_mod = (std_matrix_mod.index > hinf_idx) * (std_matrix_mod.index < hsup_idx)
#     maskT_obs = (std_matrix_obs.index > hinf_idx) * (std_matrix_obs.index < hsup_idx)
#     mod = std_matrix_mod.loc[hinf_idx:hsup_idx,:]
#     obs = std_matrix_obs.loc[hinf_idx:hsup_idx,:]
#     percentiles_mod = np.nanpercentile(mod, [25, 50, 75, 90])
#     percentiles_obs = np.nanpercentile(obs, [25, 50, 75, 90])
#     
#     #plt.figure()
#     #plt.hist(mod.flatten(),range=(10,30), normed=True, alpha=0.5)
#     #lt.hist(obs.flatten(),range=(10,30), normed=True, alpha=0.5)
#     ax = plt.subplot(2,3,indplot)  
#     ax.spines["top"].set_visible(False)  
#     ax.spines["right"].set_visible(False)  
#     ax.get_xaxis().tick_bottom()  
#     ax.get_yaxis().tick_left()     
#     matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
#     matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
#     ax.set_ylabel('norm occurrences')
#     ax.set_xlabel('${\sigma}$(IWV) [Kg/m^2]')
#     #ax.ylim(ymax)
#     plt.ylim(0.,ymax)
#     plt.xlim(xmin, xmax)
#     plt.text(0.6, ymax-1.5*ymax/10., 'median mod = '+str(round(percentiles_mod[1], 2)))
#     plt.text(0.6, ymax-2.*ymax/10., 'median obs = '+str(round(percentiles_obs[1], 2)))  
#     plt.grid(b=True, which='major', color='#666666', linestyle=':')
#     plt.hist(mod.values.flatten(), bins=nbins, normed=True, color='red', \
#              cumulative=False, range=[xmin, xmax], alpha=0.1, label='icon-lem')       
#     plt.hist(obs.values.flatten(), bins=nbins, normed=True, color='black',\
#              cumulative=False, range=[xmin, xmax], alpha=0.1, label='obs')       
#     plt.hist(mod.values.flatten(), bins=nbins, normed=True, color='red', \
#              cumulative=False, range=[xmin, xmax], histtype='step')
#     plt.hist(obs.values.flatten(), bins=nbins, normed=True, color='black', \
#              cumulative=False, range=[xmin, xmax], histtype='step')      
#     plt.legend(loc='upper left', fontsize=12, frameon=False)
#     ax.set_title(str(hourInf)+' - '+str(hourSup)+' UTC')
#     indplot= indplot+1
# fig.tight_layout()
# plt.savefig(pathFig+'sigmaIWV2_distrib_hours_stat_global.png', format='png')    
# 
# 
# =============================================================================


#%%
# IWV distributions per hour of the day
#IWV_obs_nd = np.vstack([IWV_obs_nd,np.nan*np.ones(6)]) 
hours = [0, 6, 9, 12, 15, 18, 24]

def hour2idx(hour, dtime=9):
    return int(hour*3600/dtime)

nbins = 20
ymax  = 0.2
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
xmin    = 5.
xmax    = 30.
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
    plt.hist(mod.flatten(), bins=nbins, normed=True, color='red', cumulative=False, range=[xmin, xmax], alpha=0.1, label='icon-lem')       
    plt.hist(obs.flatten(), bins=nbins, normed=True, color='black', cumulative=False, range=[xmin, xmax], alpha=0.1, label='obs')       
    plt.hist(mod.flatten(), bins=nbins, normed=True, color='red', cumulative=False, range=[xmin, xmax], histtype='step')       
    plt.hist(obs.flatten(), bins=nbins, normed=True, color='black', cumulative=False, range=[xmin, xmax], histtype='step')      
    plt.legend(loc='upper left', fontsize=12, frameon=False)
    ax.set_title(str(hourInf)+' - '+str(hourSup)+' UTC')
    indplot= indplot+1
fig.tight_layout()
plt.savefig(pathFig+'IWV_distrib_mod_hours_obs_stat_global.png', format='png')    

#=============================================================================

#%%
# =============================================================================
# plotting scatter plots of surface latent and sensible heat flux for every half hour 
# =============================================================================
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
plt.xlim(-350., 50.)
plt.ylim(-350., 50.)
plt.xlabel('Latent heat surface flux obs [W/m^2]', fontsize=16)
plt.ylabel('Latent heat surface flux icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')

cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
for indFile in range(Nfiles):
    print(len(LHSF_mod[indFile]))
    print(indFile)
    print(len(LHSF_obs[indFile]))
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
##### Evaporative fraction 
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
plt.xlim(-5., 5.)
plt.ylim(-50., 50.)
plt.xlabel('Evaporative fraction obs [W/m^2]', fontsize=16)
plt.ylabel('Evaporative fraction icon lem [W/m^2]', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')


# calculating evaporative fraction as ratio of latent/(sens+latent heat)
EvapFraction_iconlem = []
EvapFraction_obs = []

for indFile in range(Nfiles):
    print(indFile)
    SHSF_iconlemPlot = -SHSF_mod[indFile]
    SHSF_obsPlot = SHSF_obs[indFile]
    LHSF_iconlemPlot = LHSF_mod[indFile]
    LHSF_obsPlot = -LHSF_obs[indFile]
    sizeDots = 5#SHSF_err_obs[indFile]
    EvapFractionDay_iconlem = np.zeros((len(SHSF_iconlemPlot)))
    EvapFractionDay_iconlem.fill(np.nan)
    EvapFractionDay_obs = np.zeros((len(SHSF_obsPlot)))
    EvapFractionDay_obs.fill(np.nan)
    for ind in range(len(SHSF_iconlemPlot)):
        EvapFractionDay_iconlem[ind] = LHSF_iconlemPlot[ind]/(LHSF_iconlemPlot[ind]+SHSF_iconlemPlot[ind])
    for ind in range(len(SHSF_obsPlot)):        
        EvapFractionDay_obs[ind] = LHSF_obsPlot[ind]/(LHSF_obsPlot[ind]+SHSF_obsPlot[ind])
    EvapFraction_iconlem.append(EvapFractionDay_iconlem)
    EvapFraction_obs.append(EvapFractionDay_obs)
    
cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
for indFile in range(Nfiles):
    EvapFractionPlot_iconlem = EvapFraction_iconlem[indFile]
    EvapFractionPlot_obs =EvapFraction_obs[indFile]
    
    cax = ax.scatter(EvapFractionPlot_obs[:], EvapFractionPlot_iconlem[:-1], c=colors, cmap=cmap, \
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
plt.savefig(pathFig+'EvapFraction_scatterplot_obs_mod_allDays.png', format='png')
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
colors = np.arange(0,len(LWF_mod[0]))
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
    if len(LWF_iconlemPlot) < 49:
        NumberNans = 49 - len(LWF_iconlemPlot)
        LWF_iconlemPlot = np.append(np.asarray(LWF_iconlemPlot), np.repeat(np.nan, float(NumberNans)))
    if len(LWF_obsPlot) < 49:
        NumberNans = 49 - len(LWF_obsPlot)
        LWF_obsPlot = np.append(np.asarray(LWF_obsPlot), np.repeat(np.nan, float(NumberNans)))
    
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
# =============================================================================
# data = np.arange(-1000, 2000)
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
# matplotlib.rcParams['savefig.dpi'] = 100
# plt.gcf().subplots_adjust(bottom=0.15)
# fig.tight_layout()
# ax = plt.subplot(1,1,1)  
# ax.spines["top"].set_visible(False)  
# ax.spines["right"].set_visible(False)  
# ax.get_xaxis().tick_bottom()  
# ax.get_yaxis().tick_left() 
# colors = np.arange(0,len(datetime_fluxes))
# plt.plot(data, data, color='black', linestyle=':')
# plt.xlim(-20., 950.)
# plt.ylim(-150., 50.)
# plt.xlabel('Shortwave downward flux obs [W/m^2]', fontsize=16)
# plt.ylabel('Shortwave downward flux icon lem [W/m^2]', fontsize=16)
# plt.grid(b=True, which='major', color='#666666', linestyle=':')
# 
# cmap = plt.cm.get_cmap('jet', len(datetime_fluxes)) 
# for indFile in range(Nfiles):
#     print(indFile)
#     SWF_iconlemPlot = SWF_mod[indFile]
#     SWF_obsPlot = SWF_obs[indFile]
#     sizeDots = SWF_err_obs[indFile]
#     cax = ax.scatter(SWF_obsPlot[:], SWF_iconlemPlot[:], c=colors, cmap=cmap, \
#                  s=10*sizeDots)
# cbar = fig.colorbar(cax, \
#                     cmap=cmap, \
#                     ticks= [0, 8, 16, 24, 32, 40, 47])
# cbar.set_label(label='time [hh:mm]',size=15, family='helvetica')
# cbar.ax.tick_params(labelsize=14)
# cbar.ax.set_yticklabels([str(datetime_fluxes[0])[11:16],\
#                          str(datetime_fluxes[8])[11:16],\
#                          str(datetime_fluxes[16])[11:16],\
#                          str(datetime_fluxes[24])[11:16],\
#                          str(datetime_fluxes[32])[11:16],\
#                          str(datetime_fluxes[40])[11:16],\
#                          str(datetime_fluxes[47])[11:16]], fontsize=14) 
# plt.tight_layout()
# plt.savefig(pathFig+'SWF_scatterplot_obs_mod_allDays.png', format='png')
# 
# =============================================================================

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

#%%
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
        plt.xlabel('${\sigma^{2}} [m^{2}s^{-2}$]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart+1
    
    plt.savefig(pathFig+'varW_globalMeanDataset_diurnal_cycle_obs_mod.png', format='png')
