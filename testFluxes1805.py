#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:33:37 2019

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

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

from myFunctions2 import f_readingTowerData
from f_processModelOutput import f_processModelOutput
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

# ----- user defined parameters for model data reading and processing
station       = 'joyce'
debuggingFlag = 1
verboseFlag   = 1
reprocICON    = 0        # flag for reprocessing the icon output data
timeSpinOver  = 0.0      # ending time of the spin up of the model (corresponds to the time at which we start to calculate cloud fraction
intTime       = 400.     # time interval over which calculate cloud fraction and coupling parameters [seconds] corresponding to minutes with a resolution of model output of 9 sec (5 min = 50), (15 min = 150)
QcThreshold   = 10**(-7) # threshold in Qc to identify liquid cloud presence
QiThreshold   = 10**(-7) # threshold in Qc to identify ice cloud presence   
SigmaWThres   = 0.2      # threshold in the variance of the vertical velocity, used as a threshold for identifying turbulent regions.
nTimeMean     = 200      # number of time stamps to be collected to average over 30 min=1800 s
timeStep      = 33       # time step for running mean
timewindow    = 200      # time window for calculating running mean of variance corresponding to 30 min with 9 sec resolution (200*9=1800 sec = 30 min)
gradWindThr   = 0.01     # Threshold for wind gradient to determine wind shear presence in the PBLclas
timeWindowSk  = 33       # number of time stamps corresponding to 5 minutes in the ICON file
runningWindow = 200*4    # number of time stamps corresponding to 30 minutes running window for the average

# ----- creating dictionary of input parameters to process icon lem model output
modelInputParameters = {'timeSpinOverVar':timeSpinOver, 'intTimeVar':intTime, 'QcThresholdVar':QcThreshold, \
                  'QiThresholdVar':QiThreshold, 'SigmaWThresVar':SigmaWThres, 'nTimeMeanVar':nTimeMean, \
                  'timeStepVar':timeStep, 'timewindowVar':timewindow, 'gradWindThrVar':gradWindThr, \
                  'timeWindowSkVar':timeWindowSk, 'runningWindowVar':runningWindow }


# ----- define list of days to be processed 
dayList           = ['20130518']#'20130501','20130502','20130425', '20130424']#, '20130429']#]#,, '20130427',
Ndays             = len(dayList)

# ----- defining input directories for data
path_tower        = '/data/hatpro/jue/hdcp2/meteo_data/'#'/data/TR32/D2/data/juelich/met_mast/'
path_cloudnet_cat = '/data/hatpro/jue/cloudnet/juelich/processed/categorize/2013/'
path_icon         = '/data/inscape/icon/experiments/juelich/meteo-4nest/'
path_bl_class     = '/data/hatpro/jue/cloudnet/juelich/products/bl-classification/2013/'
path_mwr_joyce    = '/data/hatpro/hps/data/level2/'
path_mwr_kit      = '/data/hatpro/hpk/data/level2/'
path_mwr_lacross  = '/data/hatpro/hpl/data/level2/'
path_gps          = '/data/hatpro/jue/data/gps/'
path_radiation    = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/'
path_radiosondes  = '/home/juelich/rs_hope/KIT/'
path_LH_SH_fluxes = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
path_LWC          = '/data/hatpro/jue/cloudnet/juelich/products/lwc-scaled-adiabatic/'
patch             = 'patch003' # patch002, patch003, patch004
domSel            = 'DOM03'
pathOut           = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
# ----- defining data flags
dataFlagArr       = np.repeat(1, 14)
dataFlagLabel     = ['tower', 'cloudnetClass', 'PBLclass', 'MWR_joyce','radiosoundings', 'LWC_Cloudnet_prod']


# ----- loop on the number of days
for iDay in range(Ndays):
    
    # set day, month and year string and output path for plots for the selected day
    date              = dayList[iDay]
    yy                = dayList[iDay][0:4]
    mm                = dayList[iDay][4:6]
    dd                = dayList[iDay][6:8]
    pathDebugFig      = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'+date+'/'
    # assigning a string name to the file just created for the next steps of the code
    iconLemData           = Dataset('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/icon_lem_derivedproperties'+date+'.nc', mode='r')
    print('reprocessing or processing the day for the first time')

SHFL_iconlem      = iconLemData.groups['Temp_data'].variables['SHFL'][:].copy()
LHFL_iconlem      = iconLemData.groups['Temp_data'].variables['LHFL'][:].copy()
datetime_ICON     = nc4.num2date(iconLemData.groups['Temp_data'].variables['datetime_ICON'][:], \
                                  iconLemData.groups['Temp_data'].variables['datetime_ICON'].units) 

# resampling ICON lem variables on the same time resolution of the observations (mean half hourly)
SHFL_DF           = pd.Series(SHFL_iconlem, index=datetime_ICON)
SHFL_30min        = SHFL_DF.resample('30min', how='mean')
#SHFL_30min        = SHFL_DF.resample('30min').nanmean()    
LHFL_DF           = pd.Series(LHFL_iconlem, index=datetime_ICON)
LHFL_30min        = LHFL_DF.resample('30min').mean()        
datetime_30m      = [datetime.datetime(int(yy),int(mm),int(dd),0,0,0) + \
                    datetime.timedelta(minutes=30*x) for x in range(0, 49)]


if (len(SHFL_30min) < 48):
    NumberNans = 48 - len(SHFL_30min)
    outSerieSHFL = np.append(np.asarray(SHFL_30min.values), np.repeat(np.nan, float(NumberNans)))
    SHFL_30min = pd.Series(outSerieSHFL, index = datetime_30m[:-1])
    

if (len(LHFL_30min) < 48):
    NumberNans = 48 - len(LHFL_30min)
    outSerieLHFL = np.append(np.asarray(LHFL_30min.values), np.repeat(np.nan, float(NumberNans)))
    LHFL_30min = pd.Series(outSerieLHFL, index = datetime_30m[:-1])
    
#%%
fig, ax = plt.subplots(figsize=(10,4))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
ax.xaxis_date()
plt.plot(datetime_ICON, SHFL_iconlem)
plt.plot(datetime_30m[:-1], SHFL_30min)

#%%