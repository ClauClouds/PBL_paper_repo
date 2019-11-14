#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:50:14 2019

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
from myFunctions import f_closest
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
from cloudnetFunctions import f_calculateCloudMaskCloudnet
#from myFunctions import f_calcCloudBaseTopPBLclouds
from myFunctions import f_resamplingMatrixCloudnet
from myFunctions2 import hourDecimal_to_datetime
from myFunctions import f_selectingPBLcloudWindow
from cloudnetFunctions import f_calculateCloudMaskCloudnet


def f_calcCloudBaseTopPBLcloudsV2(cloudMask, dimTime, dimHeight, height, cloudTimeArray, time):
    # cloud mask for identifying cloud base and cloud top of PBL clouds
    meanCloudThickness = 600.
    minCBheight = 2500.
    # filtering clouds above 5000mt
    ind_above = np.where(height > 5000.)
    cloudMask[:, ind_above] = 0.
    
    
    # converting cloud mask to 1 / 0 matrices
    BinaryMatrix = np.zeros((dimTime, dimHeight))
    for itime in range(dimTime):
        for iH in range(dimHeight):
            if cloudMask[itime,iH] != 0.:
                BinaryMatrix[itime,iH] = 1
            
    # calculating gradient of binary cloud mask
    gradBinary = np.diff(BinaryMatrix,axis=1)

    # counting max number of cloud base/cloud top found 
    numberCB = []
    numberCT = []
    for itime in range(dimTime):
        column = gradBinary[itime,:]
        numberCB.append(len(np.where(column == 1.)[0][:]))
        numberCT.append(len(np.where(column ==-1.)[0][:]))  

    NCB=max(numberCB)   
    NCT=max(numberCT)
    
    # generating cloud base and cloud top arrays 
    CBarray      = np.zeros((dimTime,NCB))
    CBarray.fill(np.nan)
    CTarray      = np.zeros((dimTime,NCT))
    CTarray.fill(np.nan)
    NlayersArray = np.zeros((dimTime))
    NlayersArray.fill(np.nan)
    
    # if no cloud bases or no cloud tops are found, then CB and CT are assigned to nan
    if (NCB == 0) or (NCT == 0):
        CB_collective = np.zeros((dimTime))
        CB_collective.fill(np.nan)
        CT_collective = np.zeros((dimTime))
        CT_collective.fill(np.nan)
        CB_PBL_out = np.zeros((dimTime))
        CB_PBL_out.fill(np.nan)
        CT_PBL_out = np.zeros((dimTime))
        CT_PBL_out.fill(np.nan)
    else:
    # if some cloud base / cloud tops are found, all the found values are stored
        # storing cloud base and cloud top arrays
        for iTime in range(dimTime):
            column                    = gradBinary[iTime,:]
            indCB                     = np.where(column == -1.)[0][:]
            NfoundCB                  = len(indCB)
            indCT                     = np.where(column == 1.)[0][:] 
            NfoundCT                  = len(indCT)
            CBarray[iTime,0:NfoundCB] = height[indCB]
            CTarray[iTime,0:NfoundCT] = height[indCT]
            NlayersArray[iTime]       = numberCB[iTime]
        
        
        # we define a collective cloud base/top to consider multilayer PBL clouds as one
        # we assign min CB and max CT for each PBl cloud found.
        CB_collective = np.asarray(CBarray[:,0])
        CT_collective = np.asarray(CTarray[:,0])
        CB_PBL_out    = np.repeat(np.nan, len(time))
        CT_PBL_out    = np.repeat(np.nan, len(time))
        
        for ind in range(dimTime):
        #    if (np.isnan(CB[ind,0]) == True):
            CB_collective[ind] = np.nanmin(CBarray[ind,:])
            CT_collective[ind] = np.nanmax(CTarray[ind,:])
            print(time[ind])
            #selecting temporal window in which cloud top and base for PBL clouds have to be calculated
            if (time[ind] > cloudTimeArray[0]) * (time[ind] < cloudTimeArray[1]):
                print('time matched')
                if (CB_collective[ind] < minCBheight):
                    # for boundary layer clouds, we can assume the lowest cloud base is correct
                    # we can also assume that from the lowest cloud base, the cloud does not extend 
                    # in the vertical for more than 1500 mt. If the max cloud top is above 1500 mt then 
                    # we select among cloud tops, those that are located withing such distance from cloud base 
                    maxCTheightPBL = np.nanmin(CBarray[ind,:]) + meanCloudThickness
                    #print('max cloud top', maxCTheightPBL)
                    if (np.nanmax(CTarray[ind,:]) > maxCTheightPBL):
                        findLowerCT = np.where(CTarray[ind,:] < maxCTheightPBL)
                        if (len(findLowerCT[0]) == 0): # no elements are found below the maximum allowed height for cloud top
                            print(CBarray[ind,:])
                            print(time[ind])
                            print(len(findLowerCT[0]))
                            print(CTarray[ind,:])
                            CT_PBL_out[ind] = np.nan
                            CB_PBL_out[ind] = np.nan
                        else:
                            print('sono qui')
                            CT_PBL_out[ind] = np.nanmin(CTarray[ind,findLowerCT]) # assigning minmum cloud top
                            CB_PBL_out[ind] = CB_collective[ind] # assigning cloud base if it is below 2500 mt

                            
                        print(np.nanmax(CTarray[ind,:]))
                        print(CT_PBL_out[ind])
                        print('*******************')
        # filtering clouds in PBL using human filtering for hours
        #if  np.count_nonzero(~np.isnan(CB_collective)) != 0:
# =============================================================================
#         timeStart = cloudTimeArray[0]
#         timeEnd   = cloudTimeArray[1]
#         CB_PBL = pd.Series(np.repeat(np.nan, len(time)), index=time)
#         maskt = (CB_PBL.index > timeStart) * (CB_PBL.index < timeEnd)
#         CB_PBL.loc[maskt] = CB_collective[maskt]
#         CT_PBL = pd.Series(np.repeat(np.nan, len(time)), index=time)
#         maskt = (CT_PBL.index > timeStart) * (CT_PBL.index < timeEnd)
#         CT_PBL.loc[maskt] = CT_collective[maskt]
# =============================================================================

    
    return (CBarray, CTarray, NlayersArray, CB_PBL_out, CT_PBL_out, CB_collective, CT_collective)

date='20130427'
# set, based on human inspection, the time interval in which to select clouds
cloudTimeArray = f_selectingPBLcloudWindow(date)
print(cloudTimeArray[0])
print(cloudTimeArray[1])


iconLemData           = Dataset('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/old/icon_lem_derivedproperties'+date+'.nc', mode='r')
datetime_ICON  = nc4.num2date(iconLemData.groups['Temp_data'].variables['datetime_ICON'][:], \
                                  iconLemData.groups['Temp_data'].variables['datetime_ICON'].units) 
height_ICON    = iconLemData.groups['Temp_data'].variables['height2'][:].copy()
cloudMask_ICON = iconLemData.groups['Temp_data'].variables['cloudMask'][:].copy()
Hsurf          = 97.4#iconLemData.groups['Temp_data'].variables['HeightSurface'][:].copy()

#cloudMask, dimTime, dimHeight, height\


result = f_calcCloudBaseTopPBLclouds(cloudMask_ICON, len(datetime_ICON), \
            len(height_ICON), height_ICON, cloudTimeArray, datetime_ICON)
#(cloudMask, dimTime, dimHeight, height, cloudTimeArray, time)

#%%
timeStart = cloudTimeArray[0]
timeEnd   = cloudTimeArray[1]
CB = result[5]
CT = result[6]
N  = result[2]
CB_PBL = result[3]
CT_PBL = result[4]

#plt.plot(datetime_ICON,CB, label='all cb')
plt.plot(datetime_ICON, CB_PBL, label='pbl cb')
#plt.plot(datetime_ICON,CT, label='all ct')
plt.plot(datetime_ICON, CT_PBL, label='pbl ct')
plt.legend()
#%%
# plot of cloud mask, cloud base and cloud top normal and filtered for PBL
fig, ax = plt.subplots(figsize=(10,4))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
ax.xaxis_date()
cax = ax.pcolormesh(datetime_ICON, height_ICON, cloudMask_ICON.transpose(), vmin=0, vmax=3, cmap=plt.cm.get_cmap("RdPu", 4))
ax.set_ylim(0,12000.)                                               # limits of the y-axes
#ax.set_xlim(0,24)                                 # limits of the x-axes
ax.set_title("cloud mask", fontsize=14)
ax.set_xlabel("time ", fontsize=12)
ax.set_ylabel("height [m]", fontsize=12)
plt.plot(datetime_ICON, CB_PBL, color='black', label='CB PBL')
plt.plot(datetime_ICON, CT_PBL, color='black', linestyle=':', label='CT PBL')
plt.legend()
cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical')
cbar.ticks=([0,1,2,3])
cbar.ax.set_yticklabels(['no cloud','liquid','ice', 'mixed phase'])
cbar.set_label(label="cloud type",size=12)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=80
#plt.savefig(pathDebugFig+'cloudMask_'+stringData+'_'+date+'.png', format='png')
    
#%%
path_cloudnet_cat = '/data/hatpro/jue/cloudnet/juelich/processed/categorize/2013/'
CLOUDNET_data = Dataset(path_cloudnet_cat+date+'_juelich_categorize.nc', mode='r')

# ----- reading CLOUDNET data variables
time_CLOUDNET         = CLOUDNET_data.variables['time'][:].copy()
height_CLOUDNET       = CLOUDNET_data.variables['height'][:].copy()
datetime_CLOUDNET     = hourDecimal_to_datetime(int(2013), int(4), int(27), time_CLOUDNET)
cloudnet              = CLOUDNET_data.variables['category_bits'][:].copy()
# test the same on cloudnet
cloudnet_res = f_resamplingMatrixCloudnet(datetime_CLOUDNET, height_CLOUDNET, \
                                              cloudnet, datetime_ICON, height_ICON, cloudMask_ICON) 
CategoryCN_res    = cloudnet_res
cloudMask_cloudnet         = f_calculateCloudMaskCloudnet(datetime_ICON, height_ICON, \
                                                          CategoryCN_res.transpose().astype(int))
result_obs         = f_calcCloudBaseTopPBLclouds(cloudMask, len(datetime_ICON), len(height_ICON), \
                                                 height_ICON,cloudTimeArray, datetime_ICON)

#%%
CB = result_obs[3]
CT = result_obs[4]
N  = result_obs[2]
#plt.plot(datetime_ICON,CB, label='all cb')
plt.plot(datetime_ICON, CB, label='pbl cb')
#plt.plot(datetime_ICON,CT, label='all ct')
plt.plot(datetime_ICON, CT, label='pbl ct')
plt.legend()
#plt.plot(datetime_ICON, CB_Final, color='yellow')
#%%

fig, ax = plt.subplots(figsize=(10,4))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left() 
ax.xaxis_date()
cax = ax.pcolormesh(datetime_ICON, height_ICON, cloudMask_cloudnet.transpose(), vmin=0., vmax=3., cmap=plt.cm.get_cmap("jet", 4))
#ax.set_ylim(Hsurf,4000.)                                               # limits of the y-axes
#ax.set_xlim(0,24)                                 # limits of the x-axes
plt.plot(datetime_ICON, CB, color='white')
plt.plot(datetime_ICON, CT, color='red')

ax.set_title("PBL classification", fontsize=14)
ax.set_xlabel("time [UTC]", fontsize=12)
ax.set_ylabel("height [m]", fontsize=12)
cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3, 4], orientation='vertical')
cbar.ticks=([0,1,2,3,4])
cbar.ax.set_yticklabels(['no cloud','liquid','ice', 'ice and liquid'])
cbar.set_label(label="cloud mask",size=12)
cbar.ax.tick_params(labelsize=12)
cbar.aspect=20
#plt.savefig(pathDebugFig+'PBLclassification_iconlem_'+date+'.png', format='png')