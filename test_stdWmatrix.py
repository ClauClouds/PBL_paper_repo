#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 13:36:23 2019

@author: cacquist
"""



# ---- importing libraries
import numpy as np
import matplotlib
import scipy
import numpy.ma as ma
import pandas as pd
import netCDF4 as nc4
import glob
import datetime
from netCDF4 import Dataset
import matplotlib.dates as mdates
 
from myFunctions import f_closest
import matplotlib.pyplot as plt
from myFunctions import f_calcPblHeightRN
from myFunctions import f_calcWvariance
from myFunctions import f_runningMeanSkewnessVarianceStd_W
from myFunctions import f_PBLClass
from myFunctions import f_calcCloudBaseTopPBLcloudsV2
from myFunctions import f_calcCloudBaseTopPBLclouds
from myFunctions import f_calcPblHeightTW
from myFunctions import f_cloudmask


path = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_patch003/'
filename = 'icon_lem_derivedproperties20130502.nc'
ncdata = Dataset(path+filename, mode='r')
time          = nc4.num2date(ncdata.groups['Temp_data'].variables['datetime_ICON'][:], \
                     ncdata.groups['Temp_data'].variables['datetime_ICON'].units)
height        = ncdata.groups['Temp_data'].variables['height2'][:]
std_matrix    = ncdata.groups['Temp_data'].variables['stdWmatrix'][:]
PBL_h         = ncdata.groups['Temp_data'].variables['PBLHeightArrRN'][:]
#PBL_h2        = ncdata.groups['Temp_data'].variables['PBLHeightArrTW'][:]

# =============================================================================
PBLHeightArrTW = np.zeros((len(time)))
PBLHeightArrTW.fill(np.nan)
# 
device         = 'mod'
sigmaThreshold = 0.4 #modelInputParameters['SigmaWThresStd'] #  m/s, threshold for std of w from Schween et al, 2014.AMT
dimTime     = len(time)
# =============================================================================
# 
# 
# #std_matrix[:,height < height[142]] = 0.
# for ind in range(len(time)):
#     if device == 'mod':
#         column   = std_matrix[ind,:]
#         aboveThr = column > sigmaThreshold
#             
#         #selecting heights below 2000
#         Hsel = height[aboveThr]
#         Hbelow = Hsel[Hsel < 2000.]
#         print(Hbelow)
#         if np.count_nonzero((Hbelow)) != 0:
#             PBLHeightArrTW[ind] = np.nanmax(Hbelow)
# 
# 
# =============================================================================
PBLHeightArrTW = f_calcPblHeightTW(std_matrix,sigmaThreshold,height,time, device)
#print(PBLHeightArrTW)
# =============================================================================
pathFig = '/work/cacquist/HDCP2_S2/statistics/debug/'
#%%

# plotting Qc and Qi
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,6))
#Qrplot = np.ma.array(std_matrix[:,:], mask=std_matrix[:,:]<0.4)

Qrplot = np.ma.array(std_matrix[:,:], mask=np.isnan(std_matrix[:,:]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis_date()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel("height [m] ", fontsize=14)
ax.set_title('height [m]', fontsize=14)
ax.set_xlabel("time [hh:mm] ", fontsize=14)
cax1 = ax.pcolormesh(time, height[:], Qrplot.transpose(), vmin=np.nanmin(Qrplot), \
                     vmax=np.nanmax(Qrplot), cmap='viridis')
plt.plot(time, PBL_h, color='white')
plt.plot(time, PBLHeightArrTW, color='yellow')
plt.xlim(datetime.datetime(2013,5,2,0,15,0), datetime.datetime(2013,5,2,23,45,0))
plt.ylim(107.,4000.)
plt.legend(loc='upper left', fontsize=12, frameon=False)
#plt.xlim(timeStart,timeEnd)
cbar = fig.colorbar(cax1, orientation='vertical')
cbar.set_label(label=" $\sigma $ [$ms^{-1}$]",size=14)
cbar.ax.tick_params(labelsize=14)
cbar.aspect=80
plt.savefig(pathFig+'sigmaW_ZMLH.png', format='png')


    #column[np.where(column )]