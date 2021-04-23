


"""
Created on Tue May 26 17:16:28 2020
@author: cacquist
2) plot Thetav profiles as a function of Theta_v(2500m) - Thetav(height_i) and add PBL height to each mean profile


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

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
patch = 'patch003'

flagThermodynVarScatterPlots = 1
# thermodynamic variables
thetaV_radiosObs = []
thetaV_iconlem = []
thetaV_mwrObs = []
rh_radiosObs = []
rh_mod = []
T_mod = []
P_radiosondes = []
T_radiosondes = []
theta_v_radiosondes = []
height_radiosondes = []
time_radiosondes = []
theta_v_mod = []
date_arr = []
time_mod = []
height_mod = []
lcl_radiosondes = []
ccl_radiosondes = []
lts_radiosondes = []
pblHeight_radiosondes = []
lcl_mod = []
ccl_mod = []
lts_mod = []
pblHeightRN_mod = []
pblHeightTW_mod = []
pblHeightRN_windLidarObs = []
pblHeightTW_windLidarObs = []

# directories where data are stored
# path = '/Volumes/NO NAME/PBlcloudPaper/statistics/dataset_obs_model/'
pathObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_' + patch + '/'
# '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_' + patch + '/'
pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/' + patch + '/figures_JAMES/'
# path = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
fileListObs = sorted(glob.glob(pathObs + '*.p'))
fileListMod = sorted(glob.glob(pathMod + '*icon_lem*.nc'))
Nfiles = len(fileListObs)

for indFile in range(Nfiles):
    print(indFile)
    date = fileListObs[indFile][81:89]
    yy = int(date[0:4])
    mm = int(date[4:6])
    dd = int(date[6:8])
    date_arr.append(date)

    print('processing date ' + date)
    ncdata = Dataset(pathMod + 'icon_lem_derivedproperties' + date + '.nc', mode='r')
    time = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                        ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
    height = ncdata.groups['Temp_data'].variables['height2'][:]

    # opening the file containing all the data
    infile = open(pathObs + 'dictionaries_ModObs_' + date + '.p', 'rb')
    new_dict = pickle.load(infile, encoding='latin1')

    # reading data radiosoundings
    from myFunctions import f_reshapeRadiosondes

    RadiosondeFormatted = f_reshapeRadiosondes(new_dict)

    P_radiosondes.append(RadiosondeFormatted['P'])
    T_radiosondes.append(RadiosondeFormatted['T'])
    theta_v_radiosondes.append(RadiosondeFormatted['theta_v'])
    height_radiosondes.append(RadiosondeFormatted['height'])
    time_radiosondes.append(RadiosondeFormatted['time'])
    lcl_radiosondes.append(RadiosondeFormatted['lcl'])
    ccl_radiosondes.append(RadiosondeFormatted['ccl'])
    lts_radiosondes.append(RadiosondeFormatted['lts'])
    pblHeight_radiosondes.append(RadiosondeFormatted['pblHeight'])
    rh_radiosObs.append(RadiosondeFormatted['RH'])

    # reading model data from lem for fluxes at surface and atm indeces
    theta_v_mod.append(new_dict[5]['virtualPotentialTemperature'])
    time_mod.append(new_dict[5]['time'])
    height_mod.append(new_dict[5]['height'])
    lcl_mod.append(new_dict[5]['lclHeight'])
    ccl_mod.append(new_dict[5]['cclHeight'])
    lts_mod.append(new_dict[5]['LTS'])
    pblHeightTW_mod.append(new_dict[9]['PBLHeightTW'])
    pblHeightRN_mod.append(new_dict[9]['PBLHeightRN'])
    rh_mod.append(new_dict[5]['relativeHumidity'])
    T_mod.append(new_dict[9]['T_iconlem'])

# =============================================================================
# calculating and plotting potential temperature , temperature and relative humidity profiles
# =============================================================================
from myFunctions import f_calculateMeanThetaVModelProfiles

theta_v_dict_obs_mod_arr = f_calculateMeanThetaVModelProfiles(time_radiosondes, \
                                                              theta_v_radiosondes, \
                                                              T_radiosondes, \
                                                              rh_radiosObs, \
                                                              height_radiosondes, \
                                                              lcl_radiosondes, \
                                                              ccl_radiosondes, \
                                                              lts_radiosondes, \
                                                              pblHeight_radiosondes, \
                                                              theta_v_mod, \
                                                              T_mod, \
                                                              rh_mod, \
                                                              time_mod, \
                                                              height_mod, \
                                                              lcl_mod, \
                                                              ccl_mod, \
                                                              lts_mod, \
                                                              pblHeightRN_mod)

from myFunctions import f_calculateMeanProfilesPlotThetaVRadiosondes
result = f_calculateMeanProfilesPlotThetaVRadiosondes(theta_v_dict_obs_mod_arr, height_mod)

MatrixHourMeanProfileThetaRad = result[0]
MatrixHourStdProfileThetaRad  = result[1]
listHourDict                  = result[2]
MatrixHourMeanProfileTRad     = result[3]
MatrixHourStdProfileTRad      = result[4]
MatrixHourMeanProfileRHRad    = result[5]
MatrixHourStdProfileRHRad     = result[6]
gridHeight                    = height_mod[0]
Height_2500 = np.min(gridHeight[gridHeight > 2500])
indHeight_2500 = np.argmin(gridHeight[gridHeight > 2500])

print(np.shape(MatrixHourMeanProfileThetaRad))
# hours with radiosondes: 3,5,7,8,9,11,13,14,15,16,17,19,20,21,23
# indeces               : 0,1,2,3,4, 5, 6, 7, 8, 9,10,11,12,13,14
#%%
#indexPlot = [1,3,4,5,6,7,11]
indexPlot = [2,4,5,6,8,10,12,14]
HourPlot  = ['7:00','9:00','11:00','13:00','15:00','17:00','20:00','23:00']
lclx_obs_Arr = [-7.6,  -6.8, -5.6, -3.8, -5.15, -7., -7.5, -7.55]
lclx_mod_Arr = [-10., -7.55, -6.5,  -5.,   -5., -5., -7.15, -7.6]
flagPlotThetaVglobalProfiles = 1
matplotlib.rc('xtick', labelsize=10)  # sets dimension of ticks in the plots
matplotlib.rc('ytick', labelsize=10)  # sets dimension of ticks in the plots

if flagPlotThetaVglobalProfiles == 1:
    Ncols = 4
    Nrows = 2
    Nplots = 11

    fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14, 10))
    # matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax = 2500.
    ymin = height_mod[0][-1]
    xmin = -15.#283.
    xmax = 5.#298.
    fontSizeTitle = 16
    fontSizeX = 12
    fontSizeY = 12
    # timeTitles = [']
    # axes coordinates for plotting
    indxArr = [0, 0, 0, 0, 1, 1, 1, 1]
    indyArr = [0, 1, 2, 3, 0, 1, 2, 3]

    print('numer of plots :', len(indxArr))
    for indPlot in range(len(indxArr)):

        # reading number of profiles, pbl heights for the selected hour
        Nprofiles = listHourDict[indexPlot[indPlot]]['Nprofiles']

        pbl_rad_plot = np.nanmean(listHourDict[indexPlot[indPlot]]['pblHeight_rad_hour'])
        pbl_rad_plot_std = np.nanstd(listHourDict[indexPlot[indPlot]]['pblHeight_rad_hour'])
        pbl_mod_plot = np.nanmean(listHourDict[indexPlot[indPlot]]['pblHeight_mod_hour'])
        pbl_mod_plot_std = np.nanstd(listHourDict[indexPlot[indPlot]]['pblHeight_mod_hour'])
        print(pbl_rad_plot, pbl_mod_plot)
        stringHourProcess = HourPlot[indPlot]
        print('processing '+stringHourProcess)

        # assigning indeces for subplot positions
        indx = indxArr[indPlot]
        indy = indyArr[indPlot]

        # removing subplot box top and right lines
        ax[indx, indy].spines["top"].set_visible(False)
        ax[indx, indy].spines["right"].set_visible(False)
        ax[indx, indy].get_xaxis().tick_bottom()
        ax[indx, indy].get_yaxis().tick_left()
        ax[indx, indy].text(284, 2000., 'N = ' + str(Nprofiles), fontsize=10)


        prof_mod = listHourDict[indexPlot[indPlot]]['meanProfile_mod']
        prof_obs = MatrixHourMeanProfileThetaRad[:, indexPlot[indPlot]]
        thetaV_ref_obs = np.repeat(prof_obs[indHeight_2500], len(prof_obs))
        thetaV_ref_mod = np.repeat(prof_mod[indHeight_2500], len(prof_mod))
        std_mod = listHourDict[indexPlot[indPlot]]['stdProfileMod']
        labelHour = listHourDict[indexPlot[indPlot]]['hour']
        std_obs = MatrixHourStdProfileThetaRad[:, indexPlot[indPlot]]
        #print(np.nanmin(thetaV_ref_obs-prof_obs), np.nanmax(thetaV_ref_obs-prof_obs))
        #print(thetaV_ref_obs-prof_obs)
        ax[indx, indy].plot(-thetaV_ref_obs+prof_obs, gridHeight, label='obs ' + str(labelHour) + ' UTC', color='black')
        ax[indx, indy].plot(-thetaV_ref_mod+prof_mod, height_mod[0], label='icon-lem', color='red')
        y1 = -thetaV_ref_obs+prof_obs - std_obs
        y2 = -thetaV_ref_obs+prof_obs + std_obs
        ax[indx, indy].fill_betweenx(gridHeight, y1, y2, where=y2 > y1, facecolor='black', alpha=0.2)
        y1 = -thetaV_ref_mod+prof_mod - std_mod
        y2 = -thetaV_ref_mod+prof_mod + std_mod
        ax[indx, indy].fill_betweenx(height_mod[0], y1, y2, where=y2 > y1, facecolor='red', alpha=0.2)
        ax[indx, indy].legend(loc='upper left', fontsize=12, frameon=False)
        ax[indx, indy].set_ylim(ymin, ymax)
        ax[indx, indy].set_xlim(xmin, xmax)
        # overplotting PBL height for model and obs
        ax[indx,indy].errorbar(lclx_obs_Arr[indPlot], pbl_rad_plot, fmt='ko', xerr=0, yerr=pbl_rad_plot_std, ecolor='black')
        ax[indx,indy].errorbar(lclx_mod_Arr[indPlot], pbl_mod_plot, fmt='ro', xerr=0, yerr=pbl_mod_plot_std, ecolor='red')
        # plt.title('8:00 UTC', fontsize=fontSizeTitle)
        ax[indx, indy].set_xlabel('${\Theta_v}-{\Theta_v[2500m]}$ [K]', fontsize=fontSizeX)
        ax[indx, indy].set_ylabel('height [m]', fontsize=fontSizeY)
    fig.subplots_adjust(hspace=0.15, bottom=0.1, left=0.05)
    #ax[1, 3].set_visible(False)  # to remove last plot
    fig.savefig(pathFig + 'theta_v_AllDays.png', format='png')



