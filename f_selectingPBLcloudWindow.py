#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:41:35 2019
functions to define,for each day of the dataset, which is the interval of time 
to consider to select typical boundary layer clouds and exclude other cloud types
which can mess up the statistics. the filter is based on human selection based 
on visual inspection of cloudnet target classification and model cloud mask.
As general rule: we exclude all clouds with cloud base above 2000 mt and vertical 
extension larger than 1000 m associated with a mixed / ice phase.

'20130501']#'20130502']#, '20130501','20130424', '20130425', '20130427', '20130429'
@author: cacquist
"""
import datetime
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

