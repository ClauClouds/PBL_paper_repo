#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:20:56 2019

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

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
    
filename = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/dataset_PBLcloudPaper_ModObs_20130502.p'
infile = open(filename,'rb')
new_dict = pickle.load(infile, encoding='latin1')