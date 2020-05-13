"""
date : wednesday 17 april 2020
author: Claudia Acquistapace
goal: code extracted from PBLpaper_variance_w_profiles.py to produce cloud distributions as a function of height for
 model and observations and will be included in the paper, possibly figure This plot is supposed to justify the choice of a height at which to cut the cloud top/base
    identification:

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
from myFunctions import f_resampleArrays2StandardData
from myFunctions import f_resample2StandardData
from myFunctions import f_plotTest

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']


# setting parameters for calculating averaging and domain size of the model:
NprofilesOut                    = 24  # hourly means
timeIncrement                   = 60  # hourly means
patch                           = 'patch003'
flagPlot                        = 0

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

date_arr                        = []
# dynamical properties
cloudMask_obs                   = np.zeros((Nfiles,9600, 150))
cloudMask_mod                   = np.zeros((Nfiles,9600, 150))
CB_mod                          = np.zeros((Nfiles,9600))
CB_obs                          = np.zeros((Nfiles,9600))
CT_mod                          = np.zeros((Nfiles,9600))
CT_obs                          = np.zeros((Nfiles,9600))
cloudFractionTot_obs            = np.zeros((Nfiles,48, 150))
cloudFractionTot_mod            = np.zeros((Nfiles,48, 150))
MeanCB_mod                      = np.zeros((Nfiles,NprofilesOut))
#minCB_mod                       = np.zeros((Nfiles,NprofilesOut))
#maxCB_mod                       = np.zeros((Nfiles,NprofilesOut))
MeanCB_obs                      = np.zeros((Nfiles,NprofilesOut))
#minCB_obs                       = np.zeros((Nfiles,NprofilesOut))
#maxCB_obs                       = np.zeros((Nfiles,NprofilesOut))


for indFile in range(Nfiles):

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
    height = ncdata.groups['Temp_data'].variables['height2'][:]

    # opening the file containing all the data
    infile = open(filenameObs, 'rb')
    new_dict = pickle.load(infile, encoding='latin1')

    cloudMask_obs[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[8]['cloudMask']), time[:], \
                                                         new_dict[3]['height'], date)
    cloudMask_mod[indFile,:,:]        = f_resample2StandardData(np.asarray(new_dict[7]['cloudMask']), time[:], \
                                                            new_dict[3]['height'], date)
    cloudFractionTot_obs[indFile,:,:] = np.asarray(new_dict[8]['totalCloudFraction'])
    cloudFractionTot_mod[indFile,:,:] = np.asarray(new_dict[7]['totalCloudFraction'])
    height                            = np.asarray(new_dict[7]['heightCloudFraction'])
    CB_mod[indFile,:] = f_resampleArrays2StandardData(np.asarray(new_dict[7]['cloudBase']), time[:], date)
    CT_mod[indFile,:] = f_resampleArrays2StandardData(np.asarray(new_dict[7]['cloudTop']), time[:], date)
    CB_obs[indFile,:] = f_resampleArrays2StandardData(np.asarray(new_dict[8]['cloudBase']), time[:], date)
    CT_obs[indFile,:] = f_resampleArrays2StandardData(np.asarray(new_dict[8]['cloudTop']), time[:], date)
    
    
    
# plotting  cloud vertical distributions as a function of height
# loop on every height. For every height we take all cloud faction calculated values over 30 min and we derive an histogram of them
Nbins              = 25
q_values           = [25, 50, 75, 90]
Npercentiles       = len(q_values)
DistrCF_height_obs = np.zeros((Nbins, len(height)))
x_arr_obs          = np.zeros((Nbins+1, len(height)))
DistrCF_height_mod = np.zeros((Nbins, len(height)))
x_arr_mod          = np.zeros((Nbins+1, len(height)))
percentiles_mod    = np.zeros((Npercentiles, len(height)))
percentiles_obs    = np.zeros((Npercentiles, len(height)))
CFall_height_obs = np.zeros(len(height))
CFall_height_mod = np.zeros(len(height))
CFL_height_obs = np.zeros(len(height))
CFL_height_mod = np.zeros(len(height))
CFI_height_obs = np.zeros(len(height))
CFI_height_mod = np.zeros(len(height))

for indHeight in range(len(height)):
    arrayValues_obs                                         = cloudFractionTot_obs[:,:,indHeight].flatten()
    DistrCF_height_obs[:,indHeight], x_arr_obs[:,indHeight] = np.histogram(arrayValues_obs, bins=Nbins, range=[0.001,0.99], density=True)
    arrayValues_mod                                         = cloudFractionTot_mod[:,:,indHeight].flatten()
    DistrCF_height_mod[:,indHeight], x_arr_mod[:,indHeight] = np.histogram(arrayValues_mod, bins=Nbins, range=[0.001,0.99], density=True)

    arrayCloudMask_obs = cloudMask_obs[:,:,indHeight].flatten()
    arrayCloudMask_mod = cloudMask_mod[:,:,indHeight].flatten()
    # filtering nans
    arrayCloudMask_obs = arrayCloudMask_obs[~np.isnan(arrayCloudMask_obs)]
    arrayCloudMask_mod = arrayCloudMask_mod[~np.isnan(arrayCloudMask_mod)]
    
    
    CFall_height_obs[indHeight] = float(np.count_nonzero(arrayCloudMask_obs)/len(arrayCloudMask_obs))
    CFall_height_mod[indHeight] = float(np.count_nonzero(arrayCloudMask_mod)/len(arrayCloudMask_mod))
    CFL_height_obs[indHeight] = float(np.count_nonzero(arrayCloudMask_obs == 1)/len(arrayCloudMask_obs))
    CFL_height_mod[indHeight] = float(np.count_nonzero(arrayCloudMask_mod == 1)/len(arrayCloudMask_mod))
    CFI_height_obs[indHeight] = float(np.count_nonzero(arrayCloudMask_obs > 1)/len(arrayCloudMask_obs))
    CFI_height_mod[indHeight] = float(np.count_nonzero(arrayCloudMask_mod > 1)/len(arrayCloudMask_mod))    
percentiles_mod     = np.nanpercentile(DistrCF_height_mod, q_values, axis=0)
percentiles_obs     = np.nanpercentile(DistrCF_height_obs, q_values, axis=0)

#%%plotFlag_cloudDistr = 1
if plotFlag_cloudDistr == 1:
    Nrows = 1
    Ncols = 1
    Ymin = 107.
    Ymax = 12000.
    fontSizeTitle = 12
    fontSizeX = 10
    fontSizeY = 10
    fontSizeCbar = 10
    labelsizeaxes = 10
    cbarAspect = 10
    OccMax = 10.

    fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(5, 6))
    ax = plt.subplot(Nrows, Ncols, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    #cax = ax.pcolormesh(xc_arr_obs[:,0], height, DistrCloud_height_obs[:,:].transpose(), vmin=0., vmax=OccMax,
    #                            cmap='Purples')
    ax.plot(CFall_height_obs, height, color='black', label='obs-all')
    ax.plot(CFall_height_mod, height, color='red', label='mod-all')
    ax.plot(CFL_height_obs, height, color='black', label='obs-liquid', linestyle=':')
    ax.plot(CFL_height_mod, height, color='red', label='mod-liquid', linestyle=':')
    ax.plot(CFI_height_obs, height, color='black', label='obs-ice', linestyle='--')
    ax.plot(CFI_height_mod, height, color='red', label='mod-ice', linestyle='--')
    ax.set_ylim(Ymin, Ymax)  # limits of the y-axes
    ax.set_xlim(0.,0.3)                                 # limits of the x-axes
    #ax.set_title("", fontsize=fontSizeTitle)
    ax.set_xlabel(" cloud fraction ", fontsize=fontSizeX)
    ax.set_ylabel("height [m]", fontsize=fontSizeY)
    ax.legend(loc="lower right", fontsize=12, frameon=False)
    #cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
    #cbar.ticks = ([0, 1, 2, 3])
    #cbar.set_label(label="occurrences", size=fontSizeCbar)
    #cbar.ax.tick_params(labelsize=fontSizeCbar)

    #ax = plt.subplot(Nrows, Ncols, 2)
    #ax.spines["top"].set_visible(False)
    #ax.spines["right"].set_visible(False)
    #ax.get_xaxis().tick_bottom()
    #ax.get_yaxis().tick_left()
    #matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    #matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
    #cax = ax.pcolormesh(xc_arr_mod[:,0], height, DistrCloud_height_mod[:,:].transpose(), vmin=0., vmax=OccMax,
    #                            cmap='Purples')
    #ax.set_ylim(Ymin, Ymax)  # limits of the y-axes
    #ax.set_xlim(0.,1.)                                 # limits of the x-axes
    #ax.set_title("model", fontsize=fontSizeTitle)
    #ax.set_xlabel("cloud fraction bin ", fontsize=fontSizeX)
    #ax.set_ylabel("height [m]", fontsize=fontSizeY)
    #cbar = fig.colorbar(cax, orientation='vertical', aspect=cbarAspect)
    #cbar.ticks = ([0, 1, 2, 3])
    #cbar.set_label(label="occurrences", size=fontSizeCbar)
    #cbar.ax.tick_params(labelsize=fontSizeCbar)
    plt.tight_layout()
    plt.savefig(pathFig + 'CF_all_profile.png')

    
#%%




PlotFlag_cloudMask = 0
# plot of cloud fraction time height for all days of simulation/obs
if PlotFlag_cloudMask == 1:
    for indDay in range(len(date_arr)):
        date = date_arr[indDay]
        yy = int(date[0:4])
        mm = int(date[4:6])
        dd = int(date[6:8])
        print('plotting date :' + date)
        timePlotHalfHours = pd.date_range(date, periods=9600, freq='9s')
        timeStart = datetime.datetime(yy,mm,dd,6)
        timeEnd   = datetime.datetime(yy,mm,dd,23)
        Ncols = 1
        Nrows = 2
        Nplots = 2
        fontSizeTitle = 12
        fontSizeX = 10
        fontSizeY = 10
        fontSizeCbar = 10
        labelsizeaxes = 10
        cbarAspect = 10
        Ymax = 11000
        Ymin = 107.
        from matplotlib import rcParams
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Tahoma']
        fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(10,5))
        ax = plt.subplot(Nrows, Ncols, 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        matplotlib.rc('xtick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=labelsizeaxes)  # sets dimension of ticks in the plots
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        cax = ax.pcolormesh(timePlotHalfHours, height, cloudMask_mod[indDay,:,:].transpose(), vmin=0, vmax=3,
                            cmap=plt.cm.get_cmap("GnBu", 4))
        ax.plot(timePlotHalfHours, CB_mod[indDay,:], color='black' )
        ax.plot(timePlotHalfHours, CT_mod[indDay,:], color='red' )

        ax.set_ylim(Ymin, Ymax)  # limits of the y-axes
        ax.set_xlim(timeStart,timeEnd)                                 # limits of the x-axes
        ax.set_title("cloud mask model", fontsize=fontSizeTitle)
        ax.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
        ax.set_ylabel("height [m]", fontsize=fontSizeY)
        cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical', aspect=cbarAspect)
        cbar.ticks = ([0, 1, 2, 3])
        cbar.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
        cbar.set_label(label="cloud type", size=fontSizeCbar)
        cbar.ax.tick_params(labelsize=fontSizeCbar)

        ax1 = plt.subplot(Nrows, Ncols, 2)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.set_xlim(timeStart,timeEnd)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        cax1 = ax1.pcolormesh(timePlotHalfHours, height, cloudMask_obs[indDay,:,:].transpose(), vmin=0, vmax=3,
                            cmap=plt.cm.get_cmap("GnBu", 4))
        ax1.plot(timePlotHalfHours, CB_obs[indDay,:], color='black' )
        ax1.plot(timePlotHalfHours, CT_obs[indDay,:], color='red' )
        ax1.set_ylim(Ymin, Ymax)  # limits of the y-axes
        ax1.set_title("cloud mask obs", fontsize=fontSizeTitle)
        ax1.set_xlabel("time [hh:mm]", fontsize=fontSizeX)
        ax1.set_ylabel("height [m]", fontsize=fontSizeY)
        cbar1 = fig.colorbar(cax1, ticks=[0, 1, 2, 3], orientation='vertical', aspect=cbarAspect)
        cbar1.ticks = ([0, 1, 2, 3])
        cbar1.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
        cbar1.set_label(label="cloud type", size=fontSizeCbar)
        cbar1.ax.tick_params(labelsize=fontSizeCbar)
        plt.tight_layout()
        plt.savefig(pathFig+date+'_cloudMask.png')

#timePlotHalfHours = pd.date_range(date_arr[0], periods=9600, freq='9s')
#result = f_plotTest(cloudMask_obs[0,:,:].transpose(), timePlotHalfHours, height, 'testFunction')

