#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:32:40 2019

@author: cacquist
@goal  : plot wind mean profiles from model and obs

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
pathObs   = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/old/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
#pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathMod = pathObs
#pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathFig = pathObs

#filenameObs = 'dataset_PBLcloudPaper_ModObs_20130502.p'
#filenameMod = 'icon_lem_derivedproperties20130502.nc'

fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)
dataset_mean_variance_obs       = []
dataset_mean_variance_mod       = []
uWind_obs                       = []
uWind_mod                       = []
vWind_obs                       = []
vWind_mod                       = []
wWind_obs                       = []
wWind_mod                       = []

date_arr = []

for indFile in range(Nfiles):
    
    
     print(indFile)
     filenameMod   = fileListMod[indFile]
     filenameObs   = fileListObs[indFile]
     date          = fileListMod[indFile][91:99]
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
     
     uWind_obs.append(np.asarray(new_dict[3]['zonalWind']))
     uWind_mod.append(np.asarray(new_dict[9]['u_iconlem']))
     vWind_mod.append(np.asarray(new_dict[9]['v_iconlem']))
     vWind_obs.append(np.asarray(new_dict[3]['meridionalWind']))
     wWind_obs.append(np.asarray(new_dict[3]['verticalWind']))
     wWind_mod.append(np.asarray(new_dict[9]['w_iconlem']))

u_All_obs         = np.concatenate(uWind_obs, axis=0)
u_All_mod         = np.concatenate(uWind_mod, axis=0)
v_all_obs         = np.concatenate(vWind_obs, axis=0)
v_all_mod         = np.concatenate(vWind_mod, axis=0)
w_all_obs         = np.concatenate(wWind_obs, axis=0)
w_all_mod         = np.concatenate(wWind_mod, axis=0)



