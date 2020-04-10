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
dataset_mean_variance_obs       = []
dataset_mean_variance_mod       = []
dataset_mean_skewness_obs       = []
dataset_mean_skewness_mod       = []

# dynamical properties
varianceW_obs                   = []
varianceW_mod                   = []
stdW_obs                        = []
stdW_mod                        = []
skewnessW_obs                   = []
skewnessW_mod                   = []
date_arr                        = []
MeanCB_mod                      = np.zeros((Nfiles,NprofilesOut))
minCB_mod                       = np.zeros((Nfiles,NprofilesOut))
maxCB_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanCB_obs                      = np.zeros((Nfiles,NprofilesOut))
minCB_obs                       = np.zeros((Nfiles,NprofilesOut))
maxCB_obs                       = np.zeros((Nfiles,NprofilesOut))

datetime_out                    = []
#print(fileListMod)

for indFile in range(Nfiles):
#for indFile in range(1):

    print(indFile)

    filenameMod = fileListMod[indFile]
    filenameObs = fileListObs[indFile]
    #filenameObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/dataset_PBLcloudPaper_ModObs_20130506.p'
    #filenameMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/icon_lem_derivedproperties20130506.nc'
    date = fileListMod[indFile][87:95]
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
#    print(time[0:5])

    height = ncdata.groups['Temp_data'].variables['height2'][:]
    # w_mod         = ncdata.groups['Temp_data'].variables['vertWind']
    # varW_mod      = ncdata.groups['Temp_data'].variables['varianceW']

    # opening the file containing all the data
    infile = open(filenameObs, 'rb')
    new_dict = pickle.load(infile, encoding='latin1')

    # dinamic properties
    varianceW_obs.append(new_dict[3]['varianceW'])
    varianceW_mod.append(new_dict[9]['varianceW'])
    stdW_obs.append(new_dict[3]['stdW'])
    stdW_mod.append(new_dict[9]['stdW'])
    skewnessW_obs.append(new_dict[3]['skewnessW'])
    skewnessW_mod.append(ncdata.groups['Temp_data'].variables['skewnessW'][:])
    CB_mod = new_dict[7]['cloudBase']
    CB_obs = new_dict[8]['cloudBase']

    CB_mod_df = pd.Series(np.asarray(CB_mod), index=time)
    CB_obs_df = pd.Series(np.asarray(CB_obs), index=time)

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
        datetime_out.append(HourInf)
        indInt = indInt + 1
        CB_mod_t = CB_mod_df.loc[(CB_mod_df.index < HourSup) * (CB_mod_df.index >= HourInf)]
        CB_obs_t = CB_obs_df.loc[(CB_obs_df.index < HourSup) * (CB_obs_df.index >= HourInf)]

        MeanCB_mod[indFile, itime] = CB_mod_t.mean(skipna=True)
        if len(CB_mod_t) > 0:
            minCB_mod[indFile, itime] = np.nanmin(CB_mod_t)
            maxCB_mod[indFile, itime] = np.nanmax(CB_mod_t)
        else:
            minCB_mod[indFile, itime] = np.nan
            maxCB_mod[indFile, itime] = np.nan

        MeanCB_obs[indFile, itime] = CB_obs_t.mean(skipna=True)
        if len(CB_obs_t) > 0:
            minCB_obs[indFile, itime] = np.nanmin(CB_obs_t)
            maxCB_obs[indFile, itime] = np.nanmax(CB_obs_t)
        else:
            minCB_obs[indFile, itime] = np.nan
            maxCB_obs[indFile, itime] = np.nan



    # ---- calculating mean variance and standard deviation profiles for each hour of the day for obs and model
    varianceWmean_obs = f_calcMeanStdVarProfiles(new_dict[3]['varianceW'], time[:], height[:], \
                                                 date, yy, mm, dd, NprofilesOut, timeIncrement)
    varianceWmean_mod = f_calcMeanStdVarProfiles(new_dict[9]['varianceW'], time[:], height[:], \
                                                 date, yy, mm, dd, NprofilesOut, timeIncrement)
    skewnessWmean_obs = f_calcMeanStdVarProfiles(new_dict[3]['skewnessW'], time[:], height[:], \
                                             date, yy, mm, dd, NprofilesOut, timeIncrement)
    skewnessWmean_mod = f_calcMeanStdVarProfiles(ncdata.groups['Temp_data'].variables['skewnessW'][:], \
                                                 time[:], height[:],\
                                                date, yy, mm, dd, NprofilesOut, timeIncrement)

    dataset_mean_variance_obs.append(varianceWmean_obs)
    dataset_mean_variance_mod.append(varianceWmean_mod)
    dataset_mean_skewness_obs.append(skewnessWmean_obs)
    dataset_mean_skewness_mod.append(skewnessWmean_mod)





#%%
#=============================================================================
# calculating and plotting global mean profiles and stds of variance of vertical velocity
#=============================================================================

matrixVar_obs    = np.zeros((len(height),len(dataset_mean_variance_obs), NprofilesOut))
matrixVar_mod    = np.zeros((len(height),len(dataset_mean_variance_mod), NprofilesOut))
matrixSkn_obs    = np.zeros((len(height),len(dataset_mean_skewness_obs), NprofilesOut))
matrixSkn_mod    = np.zeros((len(height),len(dataset_mean_skewness_mod), NprofilesOut))
hours_arr        = np.zeros((NprofilesOut))

totMeanCB_mod    = np.nanmean(MeanCB_mod, axis=0)
totMeanCB_obs    = np.nanmean(MeanCB_obs, axis=0)
totMaxCB_mod    = np.nanmean(maxCB_mod, axis=0)
totMaxCB_obs    = np.nanmean(maxCB_obs, axis=0)
totMinCB_mod    = np.nanmean(minCB_mod, axis=0)
totMinCB_obs    = np.nanmean(minCB_obs, axis=0)

print(np.shape(totMaxCB_mod))

meanVariance_obs = np.zeros((len(height),NprofilesOut))
meanSkewness_obs = np.zeros((len(height),NprofilesOut))

for indHour in range(NprofilesOut):
    for indDay in range(len(dataset_mean_variance_obs)):
        matrixVar_obs[:,indDay,indHour] = dataset_mean_variance_obs[indDay]['meanProfiles'][:,indHour]
        matrixVar_mod[:,indDay,indHour] = dataset_mean_variance_mod[indDay]['meanProfiles'][:,indHour]
        matrixSkn_obs[:,indDay,indHour] = dataset_mean_skewness_obs[indDay]['meanProfiles'][:,indHour]
        matrixSkn_mod[:,indDay,indHour] = dataset_mean_skewness_mod[indDay]['meanProfiles'][:,indHour]
        hours_arr[indHour]              = indHour

#print('profiles at the following hours: ', hours_arr)
meanVariance_obs = np.nanmean(matrixVar_obs, axis=1)
stdVariance_obs  = np.nanstd(matrixVar_obs, axis=1)
meanVariance_mod = np.nanmean(matrixVar_mod, axis=1)
stdVariance_mod  = np.nanstd(matrixVar_mod, axis=1)

meanSkewness_obs = np.nanmean(matrixSkn_obs, axis=1)
stdSkewness_obs  = np.nanstd(matrixSkn_obs, axis=1)
meanSkewness_mod = np.nanmean(matrixSkn_mod, axis=1)
stdSkewness_mod  = np.nanstd(matrixSkn_mod, axis=1)

if flagPlot == 1:
    Ncols = 5
    Nrows = 2
    Nplots = 11
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
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
        plt.hlines(totMeanCB_obs[indHourPlotStart], 0.,xmax, color='black', linestyle=':', label='CB-obs')
        plt.hlines(totMeanCB_mod[indHourPlotStart], 0.,xmax, color='red', linestyle=':', label='CB-mod')
        y1        = meanVariance_obs[:,indHourPlotStart]-stdVariance_obs[:,indHourPlotStart]
        y2        = meanVariance_obs[:,indHourPlotStart]+stdVariance_obs[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
        y1        = meanVariance_mod[:,indHourPlotStart]-stdVariance_mod[:,indHourPlotStart]
        y2        = meanVariance_mod[:,indHourPlotStart]+stdVariance_mod[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
        if indPlot == 1:
            plt.legend(loc="upper right", fontsize=12, frameon=False)
        plt.ylim(ymin,ymax)
        plt.xlim(0.,xmax)
        plt.title(str(indHourPlotStart-1)+' UTC', fontsize=fontSizeTitle)
        plt.xlabel('${\sigma^{2}} [m^{2}s^{-2}$]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart+1
    
    plt.savefig(pathFig+'figure1_var_W_profiles.png', format='png')


    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 2500.
    ymin          = 107.
    xmax          = 1.5
    xmin = -1.5
    indHourPlotStart = 8 # corresponding to starting at 8 UTC

    for indPlot in range(1, Nplots):
        ax = plt.subplot(2, 5, indPlot)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        # ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=10)  # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=10)  # sets dimension of ticks in the plots
        plt.plot(meanSkewness_obs[:, indHourPlotStart], height, label='obs', color='black')
        plt.plot(meanSkewness_mod[:, indHourPlotStart], height, label='icon-lem', color='red')
        plt.hlines(totMeanCB_obs[indHourPlotStart], xmin, xmax, color='black', linestyle=':', label='CB-obs')
        plt.hlines(totMeanCB_mod[indHourPlotStart], xmin, xmax, color='red', linestyle=':', label='CB-mod')
        plt.axvline(x=0., color='grey', linestyle='-.')
        y1 = meanSkewness_obs[:, indHourPlotStart] - stdVariance_obs[:, indHourPlotStart]
        y2 = meanSkewness_obs[:, indHourPlotStart] + stdVariance_obs[:, indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2 > y1, facecolor='black', alpha=0.2)
        y1 = meanSkewness_mod[:, indHourPlotStart] - stdVariance_mod[:, indHourPlotStart]
        y2 = meanSkewness_mod[:, indHourPlotStart] + stdVariance_mod[:, indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2 > y1, facecolor='red', alpha=0.2)
        if indPlot == 1:
            plt.legend(loc="upper right", fontsize=12, frameon=False)
        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)
        plt.title(str(indHourPlotStart - 1) + ' UTC', fontsize=fontSizeTitle)
        plt.xlabel('w skewness ', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        indHourPlotStart = indHourPlotStart + 1

    plt.savefig(pathFig + 'figure2_skn_W_profiles.png', format='png')

