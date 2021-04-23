#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:57:03 2020

@author: cacquist
@goal : plot PBL cloud properties. the code reads the PBLcloud files produced with the code PBLclouds_id.py and derives 
plots of cloud properties
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
import xarray as xr
from myFunctions import f_calcMeanStdVarProfiles
from myFunctions import f_plotVarianceWSingleDays
from myFunctions import f_calcWvariance
from myFunctions import f_convertPressureToHeight
from myFunctions import f_closest
from myFunctions import f_resampleArrays2StandardData
from myFunctions import f_resample2StandardData
from myFunctions import f_plotTest
from myFunctions import lcl
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']


# directories where data are stored
patch = 'patch003'
Nlevels_example = np.arange(8)
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathOut                         = pathObs
fileListMod                     = sorted(glob.glob(pathObs+'PBLClouds_iconlem_*.nc'))
fileListObs                     = sorted(glob.glob(pathObs+'PBLClouds_Obs_*.nc'))
print(fileListMod)
print(fileListObs)

