
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
from myFunctions import f_reshapeRadiosondes
from myFunctions import f_selectingPBLcloudWindow


try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
patch                           = 'patch003'

path_radiosondes    = '/home/juelich/rs_hope/KIT/'
dayList = ['20130414', '20130420', \
 '20130424', '20130425', '20130426', '20130427', '20130428', '20130429', '20130430', \
 '20130501', '20130502', '20130503', '20130504', '20130505', '20130506', '20130509', \
 '20130510', '20130518', '20130524', '20130525', '20130527', '20130528']
# new days : = ['20130414','20130420', '20130426','20130428', '20130430','20130524','20130525','20130527', '20130528']
#dayList = ['20130414']
Ndays             = len(dayList)
height_radiosondes              = []
time_radiosondes                = []
theta_v_mod                     = []
date_arr                        = []
time_mod                        = []
height_mod                      = []
lcl_radiosondes                 = []
ccl_radiosondes                 = []

# ----- loop on the number of days
for iDay in range(Ndays):
    # set day, month and year string and output path for plots for the selected day
    date = dayList[iDay]
    yy = dayList[iDay][0:4]
    mm = dayList[iDay][4:6]
    dd = dayList[iDay][6:8]
    pathDebugFig = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_' + patch + '/' + date + '/'

    # set, based on human inspection, the time interval in which to select clouds
    humanInfo = f_selectingPBLcloudWindow(date)

    print('PROCESSING DAY:' + date)
    # -----------------------------------------------------------------------------------
    # ---- radiosoundings
    # -----------------------------------------------------------------------------------   
    pathIn = path_radiosondes + yy + mm + dd + '/'
    fileList = glob.glob(pathIn + '*.txt')
    print('number of radiosondes for the day :',len(fileList))
    for indRemove in range(len(fileList)):
        if fileList[indRemove] == 'KIT_HOPE_2013042013.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042411.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042509.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042511.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042616.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013043011.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013050311.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013050515.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013052513.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013052715.txt':
            del fileList[indRemove]

    print('number of radiosondes for the day :',len(fileList))

    if os.listdir(pathIn):
        # folder full, process radiosondes that are in
        from myFunctions import f_processRadiosondesDay
        radiosondeList = f_processRadiosondesDay(fileList, yy, mm, dd)


    new_dict =  [radiosondeList]
    # reading data radiosoundings
    from myFunctions import f_reshapeRadiosondes
    RadiosondeFormatted = f_reshapeRadiosondes(new_dict)

    height_radiosondes.append(RadiosondeFormatted['height'])
    time_radiosondes.append(RadiosondeFormatted['time'])
    lcl_radiosondes.append(RadiosondeFormatted['lcl'])
    ccl_radiosondes.append(RadiosondeFormatted['ccl'])
    print((RadiosondeFormatted['ccl']))

ccl_values = np.hstack(ccl_radiosondes)
print(ccl_values)
nbins = 30
labelsizeaxes = 10
fontSizeTitle = 12
fontSizeX = 10
fontSizeY = 10
fs = 10
fontSizeCbar = 10
cbarAspect = 10
Ymax = 11000
Ymin = 107.
# histograms parameters
histtypekeyword = 'step'
densitykeyword = True
nbins = 10
nbinsPBL = 6
bins_geo_all = 20
bins_geo_PBL = 30

fig, axes = plt.subplots(figsize=(15, 10))
axes = plt.subplot(1, 1, 1)
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)
axes.get_xaxis().tick_bottom()
axes.get_yaxis().tick_left()
matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
#axes[0,0].grid(True, which="both")
counts, yBins = np.histogram(ccl_values, range=[0., 5000.], bins=nbins, density=True)
yCenterBins = np.diff(yBins) + yBins[:-1]
axes.plot(counts, yCenterBins, color='red', label='mod')
plt.legend(frameon=False)
axes.set_ylabel('Z_ccl [m]', fontsize=fs)
axes.set_xlabel('occurrences [%]', fontsize=fs)
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/patch003/figures_JAMES/'

plt.savefig(pathFig+'CCL_histogram_obs.png')