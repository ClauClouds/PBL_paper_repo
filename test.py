#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:53:07 2019
date : 20.05.2020
author: Claudia Acquistapace
goal : correct / debug function to calculate CCL for raadiosondes and model output

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
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
patch                           = 'patch003'


def f_CCL_new(T, P, RH, height, time, date):
    """
    function to calculate convective condensation level (CCL). For more info on definitions of this level, read pp.250
    of Petty : A first course in atmospheric thermodynamics
    input: T: temperature , to be provided in K
              relative humidity, in % (es: 70.14)
              pressure, in Kpa
              device, string for "model" or "obs"
    procedure:
        step 1: calculate dew point T
        step 2: calculate saturation mixing ratio m0 at t=Td, P=Psurf
        step 3: calculate, for every value of P, Td(m=m0, P)
        step 4: check, for every level of P, if there's a level i  of P for which  T(P)i-1 < Td(m0,P)i < T(P)i+1.
            If the level is found, assign T_star = T(P)i and Z_ccl as the height corresponding to that pressure height.
        step 5: calculate Tc using adiabatic lapse rate to come back at the height of the surface.
    output: T_ccl, z_ccl

    """
    pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/' + patch + '/figures_JAMES/debugging/'

    print('calculating CCL height and T_CCL')
    # defining constants
    cost_rvl = np.power(5423, -1.)  # K
    E0       = 0.611  # Kpa
    T0       = 273.  # K
    Rv       = 461  # J K^-1 Kg^-1
    L        = 5.6 * 10**6 # J/Kg
    epsilon  = 0.622
    Ad_rate  = -9.8  # K/Km

    # assigning dimensions:
    dimHeight = len(height)
    dimTime   = len(time)

    # step 1: calculating due point temperature profile for each time (The dew point is \
    # the temperature to which air must be cooled to become saturated with water vapor. )
    # substituting RH = 0. to RH = nan to avoid log(0) cases
    RH[RH == 0.] = np.nan
    T[T == 0.]   = np.nan
    # calculating Td
    Td           = np.power(np.power(T, -1.) - cost_rvl * np.log(RH / 100.), -1.)
    
    # step 2: calculating mixing ratio at the surface for T = Td and P=Psurf
    # finding index of height corresponding to lowest level in height
    indHmin = np.nanargmin((height))
    # reading values of P, T, RH at the corresponding height
    Td_surf = Td[:, indHmin]
    P_surf  = P[:, indHmin]
    RH_surf = RH[:, indHmin]
    m0 = epsilon*E0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))) / (P_surf - E0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))))

    # step 3: calculating Td(m=m0, P) for every P value
    z_ccl  = np.zeros((dimTime))
    T_cclTop  = np.zeros((dimTime))
    z_ccl.fill(np.nan)
    #indPlotCount = 0
    for indTime in range(dimTime):
        Tdm0_profile = np.zeros((dimHeight))
        Tdm0_profile.fill(np.nan)
        indCCLprofile = []
        for indHeight in range(dimHeight-1):
            Tdm0_profile[indHeight] = 1 / ( (1/T0)  - ( (1/L) * Rv * np.log((m0[indTime] * P[indTime, indHeight])/(E0 * epsilon))))
            if Tdm0_profile[indHmin] > T[indTime, indHmin]:
                if (T[indTime, indHeight] < Tdm0_profile[indHeight]) and (T[indTime, indHeight+1] > Tdm0_profile[indHeight]):
                    indCCLprofile.append(indHeight)
            else:
                if (T[indTime, indHeight] > Tdm0_profile[indHeight]) and (T[indTime, indHeight+1] < Tdm0_profile[indHeight]):
                    indCCLprofile.append(indHeight)
                ##fig, ax = plt.subplots(figsize=(12, 5))
                #plt.plot(Tdm0_profile, height, label='TDm0')
                #plt.plot(T[indTime, :], height, label='T')
                #plt.legend()
                # plt.plot(time, z_ccl3)
                # plt.plot(np.repeat(M0[5000],len(height)), height)
                #plt.ylim(0, 6000)
                #plt.savefig(pathFig + str(indPlotCount) + 'Tm0_profile_Check.png', format='png')
                #indPlotCount = indPlotCount +1
        #print(len(indCCLprofile))
        if len(indCCLprofile) == 0:
            z_ccl[indTime] = np.nan
            T_cclTop[indTime] = np.nan
        else:
            z_ccl[indTime] = np.nanmin(height[indCCLprofile])
            T_cclTop[indTime] = np.nanmin(T[indTime,np.nanargmin(height[indCCLprofile])])



    #fig, ax = plt.subplots(figsize=(12,5))
    #plt.plot(time, z_ccl)
    #plt.ylim(0,6000)
    #plt.savefig(pathFig+date+'_z_ccl_mod.png', format='png')

    # ---- finding z(CCL) using the dry adiabatic lapse rate
    T_ground_CCL = np.zeros((dimTime))
    for indTime in range(dimTime):
        T_ground_CCL[indTime] = (T_cclTop[indTime] - Ad_rate * z_ccl[indTime] * 10. ** (-3))

    # providing output as standardized xarray output format
    DatasetOut = xr.Dataset(
        data_vars={'z_ccl'    : (('time'), z_ccl),
                   't_ccltop' : (('time'), T_cclTop),
                   't_ccl'    : (('time'), T_ground_CCL),
                   'T_dew'    : (('time', 'height'), Td)},
           coords={'time'     : time,
                   'height'   : height})

    return (DatasetOut)





# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*icon_lem*.nc'))
Nfiles                          = len(fileListObs)
date_arr = []
CCL_list = []
height_radiosondes = []
ccl_radiosondes = []
time_radiosondes = []
for indFile in range(1):
    # for indFile in range(1):
    print(indFile)

    date = fileListObs[indFile][81:89]
    print('processing date ' + date)
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)
    timeReferenceFormat = pd.date_range(date, periods=9600, freq='9s')

    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(pathMod + 'icon_lem_derivedproperties' + date + '.nc', mode='r')
    time = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)

    height = ncdata.groups['Temp_data'].variables['height2'][:]
    P = ncdata.groups['Temp_data'].variables['P'][:]
    T = ncdata.groups['Temp_data'].variables['T'][:]
    RH = ncdata.groups['Temp_data'].variables['RH'][:]

    # calling function to calculate CCL height for the day
    #CCLDataset = f_CCL_new(T, P, RH, height, time, date)

    # reformatting data to the standard datetime output
    timeStandard = pd.date_range(start=datetime.datetime(yy, mm, dd, 6, 0, 0), \
                                 end=datetime.datetime(yy, mm, dd, 23, 59, 59), freq='9s')
    #CCLstandard = CCLDataset.reindex({'time': timeStandard})

    #CCL_list.append(CCLstandard)

    # opening the file containing all the data
    infile = open(pathObs + 'dictionaries_ModObs_' + date + '.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')
    RadiosondeFormatted = f_reshapeRadiosondes(new_dict)
    height_radiosondes.append(RadiosondeFormatted['height'])
    time_radiosondes.append(RadiosondeFormatted['time'])
    ccl_radiosondes.append(RadiosondeFormatted['ccl'])


# concatenating cloud base and top datasets along the time dimension
#CCL_dataset = xr.concat(CCL_list, dim='time')
print(ccl_radiosondes)
print(time_radiosondes)

strasuka
Z_cclAll = CCL_dataset.z_ccl.values.flatten()
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
counts, yBins = np.histogram(Z_cclAll, range=[0., 5000.], bins=nbins, density=True)
yCenterBins = np.diff(yBins) + yBins[:-1]
axes.plot(counts, yCenterBins, color='red', label='mod')
plt.legend(frameon=False)
axes.set_ylabel('Z_ccl [m]', fontsize=fs)
axes.set_xlabel('occurrences [%]', fontsize=fs)

plt.savefig(pathFig+'CCL_histogram_model.png')