"""
date : wednesday 8 april 2020
author: Claudia Acquistapace
goal: code extracted from PBLpaper_statisticAnalysis_obsmod.py to produce the figure 1 of the paper, hourly profiles of
variance of vertical velocity and skewness for model and observations from radiosondes.
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
flagPlot                        = 1

# directories where data are stored
#path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
#'/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod                         = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
pathFig                         = '/work/cacquist/HDCP2_S2/statistics/figs/'+patch+'/figures_JAMES/'
#path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs                     = sorted(glob.glob(pathObs+'*.p'))
fileListMod                     = sorted(glob.glob(pathMod+'*.nc'))
Nfiles                          = len(fileListMod)
date_arr = []

# surface fluxes
SHF_mod                      = np.zeros((Nfiles,NprofilesOut))
SHF_obs                      = np.zeros((Nfiles,NprofilesOut))
LHF_mod                      = np.zeros((Nfiles,NprofilesOut))
LHF_obs                      = np.zeros((Nfiles,NprofilesOut))
minSHF_mod                       = np.zeros((Nfiles,NprofilesOut))
maxSHF_mod                       = np.zeros((Nfiles,NprofilesOut))
minSHF_obs                       = np.zeros((Nfiles,NprofilesOut))
maxSHF_obs                       = np.zeros((Nfiles,NprofilesOut))
minLHF_mod                       = np.zeros((Nfiles,NprofilesOut))
maxLHF_mod                       = np.zeros((Nfiles,NprofilesOut))
minLHF_obs                       = np.zeros((Nfiles,NprofilesOut))
maxLHF_obs                       = np.zeros((Nfiles,NprofilesOut))
MeanSHF_mod  = np.zeros((Nfiles,NprofilesOut))
MeanSHF_obs  = np.zeros((Nfiles,NprofilesOut))
MeanLHF_mod  = np.zeros((Nfiles,NprofilesOut))
MeanLHF_obs  = np.zeros((Nfiles,NprofilesOut))
MeanRH_obs = np.zeros((Nfiles,NprofilesOut))
MeanRH_mod = np.zeros((Nfiles,NprofilesOut))
minRH_mod                       = np.zeros((Nfiles,NprofilesOut))
maxRH_mod                       = np.zeros((Nfiles,NprofilesOut))
minRH_obs                       = np.zeros((Nfiles,NprofilesOut))
maxRH_obs                       = np.zeros((Nfiles,NprofilesOut))

# cloud macroscopic properties
# cloud base
MeanCB_mod                      = np.zeros((Nfiles,NprofilesOut))
minCB_mod                       = np.zeros((Nfiles,NprofilesOut))
maxCB_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanCB_obs                      = np.zeros((Nfiles,NprofilesOut))
minCB_obs                       = np.zeros((Nfiles,NprofilesOut))
maxCB_obs                       = np.zeros((Nfiles,NprofilesOut))

# cloud top
MeanCT_mod                      = np.zeros((Nfiles,NprofilesOut))
minCT_mod                       = np.zeros((Nfiles,NprofilesOut))
maxCT_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanCT_obs                      = np.zeros((Nfiles,NprofilesOut))
minCT_obs                       = np.zeros((Nfiles,NprofilesOut))
maxCT_obs                       = np.zeros((Nfiles,NprofilesOut))

#chord length
MeanCL_mod                      = np.zeros((Nfiles,NprofilesOut))
minCL_mod                       = np.zeros((Nfiles,NprofilesOut))
maxCL_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanCL_obs                      = np.zeros((Nfiles,NprofilesOut))
minCL_obs                       = np.zeros((Nfiles,NprofilesOut))
maxCL_obs                       = np.zeros((Nfiles,NprofilesOut))

# cloud thickness
MeanThickness_mod                      = np.zeros((Nfiles,NprofilesOut))
minThickness_mod                       = np.zeros((Nfiles,NprofilesOut))
maxThickness_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanThickness_obs                      = np.zeros((Nfiles,NprofilesOut))
minThickness_obs                       = np.zeros((Nfiles,NprofilesOut))
maxThickness_obs                       = np.zeros((Nfiles,NprofilesOut))
datetime_out                    = []
start = datetime.datetime(2000, 1, 1)
datetime_out = np.array([start + datetime.timedelta(hours=i) for i in range(24)])
print(datetime_out)

#print(fileListMod)
hours_arr = []
for indFile in range(Nfiles):
#for indFile in range(1):

    print(indFile)

    filenameMod = fileListMod[indFile]
    filenameObs = fileListObs[indFile]
    #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
    #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'
    date = fileListMod[indFile][87:95]
    print(date)
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])

    date_arr.append(date)
    print('processing date ' + date)

    # reading time and height from ncdf file (grid of ICON LEM
    # ( ( sec resolution, 9600 time stamps, and 150 height levels)))
    ncdata = Dataset(filenameMod, mode='r')
    time   = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:],\
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)

    height = ncdata.groups['Temp_data'].variables['height2'][:]
    # w_mod         = ncdata.groups['Temp_data'].variables['vertWind']
    # varW_mod      = ncdata.groups['Temp_data'].variables['varianceW']

    # opening the file containing all the data
    infile = open(filenameObs, 'rb')
    new_dict = pickle.load(infile, encoding='latin1')

    # cloud macroscopic properties
    # cloud base
    CB_mod = new_dict[7]['cloudBase']
    CB_obs = new_dict[8]['cloudBase']
    # cloud top
    CT_mod = new_dict[7]['cloudTop']
    CT_obs = new_dict[8]['cloudTop']

    # relative humidity
    RH_obs = new_dict[6]['relativeHumidity'][:,-1]
    RH_mod = new_dict[5]['relativeHumidity'][:,-1]
    timeRH = new_dict[6]['time']
    height_RH = new_dict[6]['height']
    print(height[0])
    print(height[148])
    print(len(height))
    print('strasuka')


#    CL_mod = new_dict[7]['chordLength']
#    CL_obs = new_dict[8]['chordLength']
    # cloud thickness
    thickness_mod = new_dict[7]['cloudThicknessAll']
    thickness_obs = new_dict[8]['cloudThicknessAll']

    # fluxes
    SHF_mod = new_dict[10]['SHF_iconlem']
    SHF_obs = new_dict[10]['SHF_obs']
    LHF_mod = new_dict[10]['LHF_iconlem']
    LHF_obs = new_dict[10]['LHF_obs']
    timeFluxes = new_dict[10]['datetime_30m']

    CB_mod_df = pd.Series(np.asarray(CB_mod), index=time)
    CB_obs_df = pd.Series(np.asarray(CB_obs), index=time)
    CT_mod_df = pd.Series(np.asarray(CT_mod), index=time)
    CT_obs_df = pd.Series(np.asarray(CT_obs), index=time)
    thickness_mod_df = pd.Series(np.asarray(thickness_mod), index=time)
    thickness_obs_df = pd.Series(np.asarray(thickness_obs), index=time)

    SHF_mod_df = pd.Series(np.asarray(SHF_mod[:-1]), index=timeFluxes[:-1])
    SHF_obs_df = pd.Series(np.asarray(SHF_obs), index=timeFluxes[:-1])
    LHF_mod_df = pd.Series(np.asarray(LHF_mod[:-1]), index=timeFluxes[:-1])
    LHF_obs_df = pd.Series(np.asarray(LHF_obs), index=timeFluxes[:-1])

    RH_mod_df = pd.Series(np.asarray(RH_mod), index=timeRH)
    RH_obs_df = pd.Series(np.asarray(RH_obs), index=timeRH)

    # calculating hourly mean, max, min CB
    timeIncrement = 60 # minutes
    deltaT = datetime.timedelta(minutes=timeIncrement)
    indInt       = 0
    for itime in range(0, NprofilesOut):
        if indInt == 0:
            HourInf = datetime.datetime(time[0].year, time[0].month, time[0].day, 0, 0, 0)
        else:
            HourInf = HourInf + deltaT
        HourSup = HourInf + deltaT
        indInt = indInt + 1
        CB_mod_t = CB_mod_df.loc[(CB_mod_df.index < HourSup) * (CB_mod_df.index >= HourInf)]
        CB_obs_t = CB_obs_df.loc[(CB_obs_df.index < HourSup) * (CB_obs_df.index >= HourInf)]
        CT_mod_t = CT_mod_df.loc[(CT_mod_df.index < HourSup) * (CT_mod_df.index >= HourInf)]
        CT_obs_t = CT_obs_df.loc[(CT_obs_df.index < HourSup) * (CT_obs_df.index >= HourInf)]
        thickness_mod_t = thickness_mod_df.loc[(thickness_mod_df.index < HourSup) * (thickness_mod_df.index >= HourInf)]
        thickness_obs_t = thickness_obs_df.loc[(thickness_obs_df.index < HourSup) * (thickness_obs_df.index >= HourInf)]
        SHF_mod_t = SHF_mod_df.loc[(SHF_mod_df.index < HourSup) * (SHF_mod_df.index >= HourInf)]
        SHF_obs_t = SHF_obs_df.loc[(SHF_obs_df.index < HourSup) * (SHF_obs_df.index >= HourInf)]
        LHF_mod_t = LHF_mod_df.loc[(LHF_mod_df.index < HourSup) * (LHF_mod_df.index >= HourInf)]
        LHF_obs_t = LHF_obs_df.loc[(LHF_obs_df.index < HourSup) * (LHF_obs_df.index >= HourInf)]

        RH_obs_t = RH_obs_df.loc[(RH_obs_df.index < HourSup) * (RH_obs_df.index >= HourInf)]
        RH_mod_t = RH_mod_df.loc[(RH_mod_df.index < HourSup) * (RH_mod_df.index >= HourInf)]

        MeanCB_mod[indFile, itime] = CB_mod_t.mean(skipna=True)
        MeanCT_mod[indFile, itime] = CT_mod_t.mean(skipna=True)
        MeanThickness_mod[indFile, itime] = thickness_mod_t.mean(skipna=True)
        MeanSHF_mod[indFile, itime] = SHF_mod_t.mean(skipna=True)
        MeanLHF_mod[indFile, itime] = LHF_mod_t.mean(skipna=True)
        MeanRH_obs[indFile, itime] = RH_obs_t.mean(skipna=True)
        MeanRH_mod[indFile, itime] = RH_mod_t.mean(skipna=True)

        if len(RH_mod_t) > 0:
            minRH_mod[indFile, itime] = np.nanmin(RH_mod_t)
            maxRH_mod[indFile, itime] = np.nanmax(RH_mod_t)
        else:
            minRH_mod[indFile, itime] = np.nan
            maxRH_mod[indFile, itime] = np.nan

        if len(RH_obs_t) > 0:
            minRH_obs[indFile, itime] = np.nanmin(RH_obs_t)
            maxRH_obs[indFile, itime] = np.nanmax(RH_obs_t)
        else:
            minRH_obs[indFile, itime] = np.nan
            maxRH_obs[indFile, itime] = np.nan

        if len(SHF_mod_t) > 0:
            minSHF_mod[indFile, itime] = np.nanmin(SHF_mod_t)
            maxSHF_mod[indFile, itime] = np.nanmax(SHF_mod_t)
        else:
            minSHF_mod[indFile, itime] = np.nan
            maxSHF_mod[indFile, itime] = np.nan

        if len(LHF_mod_t) > 0:
            minLHF_mod[indFile, itime] = np.nanmin(LHF_mod_t)
            maxLHF_mod[indFile, itime] = np.nanmax(LHF_mod_t)
        else:
            minLHF_mod[indFile, itime] = np.nan
            maxLHF_mod[indFile, itime] = np.nan

        if len(CB_mod_t) > 0:
            minCB_mod[indFile, itime] = np.nanmin(CB_mod_t)
            maxCB_mod[indFile, itime] = np.nanmax(CB_mod_t)
        else:
            minCB_mod[indFile, itime] = np.nan
            maxCB_mod[indFile, itime] = np.nan

        if len(CT_mod_t) > 0:
            minCT_mod[indFile, itime] = np.nanmin(CT_mod_t)
            maxCT_mod[indFile, itime] = np.nanmax(CT_mod_t)
        else:
            minCT_mod[indFile, itime] = np.nan
            maxCT_mod[indFile, itime] = np.nan

        if len(thickness_mod_t) > 0:
            minThickness_mod[indFile, itime] = np.nanmin(thickness_mod_t)
            maxThickness_mod[indFile, itime] = np.nanmax(thickness_mod_t)
        else:
            minThickness_mod[indFile, itime] = np.nan
            maxThickness_mod[indFile, itime] = np.nan

        MeanCB_obs[indFile, itime] = CB_obs_t.mean(skipna=True)
        MeanCT_obs[indFile, itime] = CT_obs_t.mean(skipna=True)
        MeanThickness_obs[indFile, itime] = thickness_obs_t.mean(skipna=True)
        MeanSHF_obs[indFile, itime] = SHF_obs_t.mean(skipna=True)
        MeanLHF_obs[indFile, itime] = LHF_obs_t.mean(skipna=True)


        if len(SHF_obs_t) > 0:
            minSHF_obs[indFile, itime] = np.nanmin(SHF_obs_t)
            maxSHF_obs[indFile, itime] = np.nanmax(SHF_obs_t)
        else:
            minSHF_obs[indFile, itime] = np.nan
            maxSHF_obs[indFile, itime] = np.nan

        if len(LHF_obs_t) > 0:
            minLHF_obs[indFile, itime] = np.nanmin(LHF_obs_t)
            maxLHF_obs[indFile, itime] = np.nanmax(LHF_obs_t)
        else:
            minLHF_obs[indFile, itime] = np.nan
            maxLHF_obs[indFile, itime] = np.nan

        if len(CB_obs_t) > 0:
            minCB_obs[indFile, itime] = np.nanmin(CB_obs_t)
            maxCB_obs[indFile, itime] = np.nanmax(CB_obs_t)
        else:
            minCB_obs[indFile, itime] = np.nan
            maxCB_obs[indFile, itime] = np.nan

        if len(CT_obs_t) > 0:
            minCT_obs[indFile, itime] = np.nanmin(CT_obs_t)
            maxCT_obs[indFile, itime] = np.nanmax(CT_obs_t)
        else:
            minCT_obs[indFile, itime] = np.nan
            maxCT_obs[indFile, itime] = np.nan

        if len(thickness_obs_t) > 0:
            minThickness_obs[indFile, itime] = np.nanmin(thickness_obs_t)
            maxThickness_obs[indFile, itime] = np.nanmax(thickness_obs_t)
        else:
            minThickness_obs[indFile, itime] = np.nan
            maxThickness_obs[indFile, itime] = np.nan

maxArr_mod = [maxCB_mod, maxCT_mod, maxThickness_mod, -maxSHF_mod, -maxLHF_mod]
minArr_mod = [minCB_mod, minCT_mod, minThickness_mod, -minSHF_mod, -minLHF_mod]
meanArr_mod = [MeanCB_mod, MeanCT_mod, MeanThickness_mod, -MeanSHF_mod, -MeanLHF_mod]
maxArr_obs = [maxCB_obs, maxCT_obs, maxThickness_obs, maxSHF_obs, maxLHF_obs]
minArr_obs = [minCB_obs, minCT_obs, minThickness_obs, minSHF_obs, minLHF_obs]
meanArr_obs = [MeanCB_obs, MeanCT_obs, MeanThickness_obs, MeanSHF_obs, MeanLHF_obs]
YlabelArr = ['cloud base [m]', 'cloud top [m]', 'geometrical thickness [m]', 'sensible heat flux', 'latent heat flux']
yminArr = [1000., 1000., 0., -50., 0.]
ymaxArr = [3000., 3500., 2000., 300., 250.]
Nrows = 1
Ncols = 5
Nplots = 5

print(datetime_out)
fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(10,16))
#matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
fontSizeTitle = 16
fontSizeX     = 14
fontSizeY     = 14
#timeTitles = [']

for indPlot in range(1, Nplots+1):
    ax        = plt.subplot(5,1,indPlot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis_date()
    plt.plot(datetime_out, np.nanmean(maxArr_obs[indPlot-1], axis=0), label='max obs', color='black', linestyle=':')
    plt.plot(datetime_out, np.nanmean(minArr_obs[indPlot-1], axis=0), label='min obs', color='black', linestyle=':')
    plt.plot(datetime_out, np.nanmean(meanArr_obs[indPlot-1], axis=0), label='mean obs', color='black')
    plt.plot(datetime_out, np.nanmean(maxArr_mod[indPlot-1], axis=0), label='max mod', color='red', linestyle=':')
    plt.plot(datetime_out, np.nanmean(minArr_mod[indPlot-1], axis=0), label='min mod', color='red', linestyle=':')
    plt.plot(datetime_out, np.nanmean(meanArr_mod[indPlot-1], axis=0), label='mean mod', color='red')
    if indPlot == 5:
        plt.legend()
        plt.hlines(0, datetime_out[0], datetime_out[-1], linestyle=':', color='grey')
    if indPlot == 4:
        plt.hlines(0, datetime_out[0], datetime_out[-1], linestyle=':', color='grey')
    plt.ylim(yminArr[indPlot-1], ymaxArr[indPlot-1])
    plt.xlim(datetime_out[5], datetime_out[17])
    print(YlabelArr[indPlot-1])
    #plt.title('Observations', fontsize=fontSizeTitle)
    plt.xlabel('time [hh:mm]', fontsize=fontSizeX)
    plt.ylabel(YlabelArr[indPlot-1], fontsize=fontSizeY)
    fig.tight_layout()
plt.savefig(pathFig+'figure8_cloudProperties_vs_SHF.png', format='png')