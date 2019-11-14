#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:58:40 2019
code to read all surface fluxes data for the days of the dataset.te
@author: cacquist
"""

import numpy as np
import matplotlib
import scipy
import pylab
import netCDF4 as nc4
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2num
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
from myFunctions import f_resamplingfield
import xlrd
from matplotlib import pyplot as pl

dateArr = ['20130424', '20130425', '20130427', '20130429','20130501','20130502' ]
dateTime_fluxes  = [datetime.datetime(2013,1,1,0,0,0) + datetime.timedelta(minutes=30*x) for x in range(0, 17519)]
filenameList     = ['RU_EC_001_fluxes_2013.xlsx','ME_EC_001_fluxes_2013.xlsx','SE_EC_001_fluxes_2013.xlsx']
StringStatArr    = ['RU', 'ME', 'SE']
path             = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
pathFig          = '/work/cacquist/HDCP2_S2/statistics/figs/'


for indDate in range(len(dateArr)):
    date             = dateArr[indDate]
    print('processing day ' +date)
    
    yyyy             = date[0:4]
    mm               = date[4:6]
    dd               = date[6:8]
    if mm[0] == '0':
        mm = mm[-1]
    if dd[0] == '0':
        dd = dd[-1]
    yyyy             = int(yyyy)
    mm               = int(mm)
    dd               = int(dd)
    
    
    datetimeDay      = [datetime.datetime(yyyy,mm,dd,0,0,0) + datetime.timedelta(minutes=30*x) for x in range(0, 48)]
            
    
    for indFile in range(len(filenameList)):
        fileName         = filenameList[indFile]
        stringStat       = StringStatArr[indFile]
        print('processing filename ' +fileName)
        book             = xlrd.open_workbook(path+fileName)
        data_sheet       = book.sheet_by_index(2)
        
        SHF              = data_sheet.col_values(37)
        LHF              = data_sheet.col_values(39)
        SHF_flag         = data_sheet.col_values(43)
        LHF_flag         = data_sheet.col_values(45)
        SHF_relErr       = data_sheet.col_values(54)
        LHF_relErr       = data_sheet.col_values(55)
        
        # filtering nan values and values flagged by the quality check and calculating absolute error
        Ntime = len(dateTime_fluxes)
        for ind in range(Ntime):
            if SHF[ind] == 42:
                SHF[ind] = np.nan
            if LHF[ind] == 42:
                LHF[ind] = np.nan
            if SHF_flag[ind] == 2:
                SHF[ind] = np.nan
            if LHF_flag[ind] == 2:
                LHF[ind] = np.nan

        Err_LHF          = np.asarray(LHF[1:]) * np.asarray(LHF_relErr[1:])
        Err_SHF          = np.asarray(SHF[1:]) * np.asarray(SHF_relErr[1:]) 

        # selecting only time stamps corresponding to the selected day
        LHF_dataframe        = pd.Series(LHF[1:], index=dateTime_fluxes)
        SHF_dataframe        = pd.Series(SHF[1:], index=dateTime_fluxes)
        LHF_Err_dataframe    = pd.Series(Err_LHF, index=dateTime_fluxes)
        SHF_Err_dataframe    = pd.Series(Err_SHF, index=dateTime_fluxes)        
        
        mask_t = (LHF_dataframe.index >= datetimeDay[0]) * (LHF_dataframe.index <= datetimeDay[-1])
        LHF_day = LHF_dataframe[mask_t].values
        LHF_Err_day = LHF_Err_dataframe[mask_t].values
        
        mask_t = (SHF_dataframe.index >= datetimeDay[0]) * (SHF_dataframe.index <= datetimeDay[-1])
        SHF_day = SHF_dataframe[mask_t].values
        SHF_Err_day = SHF_Err_dataframe[mask_t].values

        
        if stringStat == 'RU':
            SHF_RU              = SHF_day
            Err_SHF_RU          = SHF_Err_day
            LHF_RU              = LHF_day
            Err_LHF_RU          = LHF_Err_day
            string_RU           = stringStat
    
        if stringStat == 'ME':
            SHF_ME              = SHF_day
            Err_SHF_ME          = SHF_Err_day
            LHF_ME              = LHF_day
            Err_LHF_ME          = LHF_Err_day
            string_ME           = stringStat
    
        if stringStat == 'SE':
            SHF_SE              = SHF_day
            Err_SHF_SE          = SHF_Err_day
            LHF_SE              = LHF_day
            Err_LHF_SE          = LHF_Err_day
            string_SE           = stringStat 
            
    
    # calculating mean fluxes over the three stations for the selected day 
    SHF_mean = []
    LHF_mean = [] 
    SHF_std  = []
    LHF_std  = []  
         
    for indTime in range(len(datetimeDay)):
        SHF_mean.append(np.nanmean([SHF_SE[indTime], SHF_ME[indTime], SHF_RU[indTime]]))
        LHF_mean.append(np.nanmean([LHF_SE[indTime], LHF_ME[indTime], LHF_RU[indTime]]))
        SHF_std.append(np.nanstd([SHF_SE[indTime], SHF_ME[indTime], SHF_RU[indTime]]))
        LHF_std.append(np.nanstd([LHF_SE[indTime], LHF_ME[indTime], LHF_RU[indTime]]))            
    
    
#%%       
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(9,8))
    label_size = 16
    #mpl.rcParams['xtick.labelsize'] = label_size 
    #mpl.rcParams['ytick.labelsize'] = label_size
    matplotlib.rcParams['savefig.dpi'] = 100
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ax = plt.subplot(211)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()  
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.ylim(-100., 150.)
    plt.title(' sensible heat fluxes from 3 different stations - '+date, fontsize=18)
    plt.ylabel(' fluxes [W/m^2]', fontsize=16)
    plt.xlabel('time [hh:mm]', fontsize=16)
    plt.axhline(y=0., color='black', linestyle='--')
    plt.xlim(datetime.datetime(yyyy,mm,dd,1,0,0),datetime.datetime(yyyy,mm,dd,23,0,0))
    
    #pl.fill_between(datetimeDay, SHF_SE[:]-Err_SHF_SE[:], SHF_SE[:]+Err_SHF_SE[:], alpha=0.1, color = 'red')
    plt.plot(datetimeDay, SHF_SE, label='sensible heat flux - '+string_SE, color='red')
    
    #pl.fill_between(datetimeDay, SHF_ME[:]-Err_SHF_ME[:], SHF_ME[:]+Err_SHF_ME[:], alpha=0.1, color = 'blue')
    plt.plot(datetimeDay, SHF_ME, label='sensible heat flux - '+string_ME, color='blue')
    
    #pl.fill_between(datetimeDay, SHF_RU[:]-Err_SHF_RU[:], SHF_RU[:]+Err_SHF_RU[:], alpha=0.1, color = 'green')
    plt.plot(datetimeDay, SHF_RU, label='sensible heat flux - '+string_RU, color='green')
    plt.plot(datetimeDay, SHF_mean,label='sensible heat flux - mean ', color = 'black' )
    pl.fill_between(datetimeDay, np.asarray(SHF_mean)[:]-np.asarray(SHF_std)[:],\
                    np.asarray(SHF_mean)[:]+np.asarray(SHF_std)[:], alpha=0.1, color = 'black')
    plt.legend(frameon=False)
    plt.tight_layout()
    ax1 = plt.subplot(212)  
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)  
    ax1.get_xaxis().tick_bottom()  
    ax1.get_yaxis().tick_left()  
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.ylim(-100., 350.)
    plt.title(' Latent heat fluxes from 3 different stations - '+date, fontsize=18)
    plt.ylabel(' fluxes [W/m^2]', fontsize=16)
    plt.xlabel('time [hh:mm]', fontsize=16)
    plt.axhline(y=0., color='black', linestyle='--')
    plt.xlim(datetime.datetime(yyyy,mm,dd,1,0,0),datetime.datetime(yyyy,mm,dd,23,0,0))
    plt.plot(datetimeDay, LHF_SE, label='latent heat flux - '+string_SE, color='red')
    plt.plot(datetimeDay, LHF_ME, label='latent heat flux - '+string_ME, color='blue')
    plt.plot(datetimeDay, LHF_RU, label='latent heat flux - '+string_RU, color='green')
    plt.plot(datetimeDay, LHF_mean,label='latent heat flux - mean ', color = 'black' )
    pl.fill_between(datetimeDay, np.asarray(LHF_mean)[:]-np.asarray(LHF_std)[:],\
                    np.asarray(LHF_mean)[:]+np.asarray(LHF_std)[:], alpha=0.1, color = 'black')
    plt.legend(frameon=False)
    plt.tight_layout()
    
    
    
    
    plt.savefig(pathFig+'sensibleLatentHeatFluxes_'+date+'_allStations.png', format='png')
        
    
#%% 
    print('saving file for day '+date)
    units = "seconds since "+str(yyyy)+'-'+str(mm)+'-'+str(dd)+' '+"00:00:00"
    calendar = 'standard'
    # ------- saving mean outputs for the selected day as ncdf for Joyce
    f                 = nc4.Dataset(path+'meanSurface_LHSH_Fluxes_ME_RU_SE_'+date+'.nc',\
                                    'w', format='NETCDF4')      # creates a netCDF file for writing 
        
    # specify dimensions of the data ( each dimension of multidimensiona array needs to be given a name and a length)
    f.createDimension('dimT', len(datetimeDay))             # dimension for time
    
    # preallocating netCDF variables for data storage
    LHF_RU_ncdf       = f.createVariable('LHF_RU_ncdf', 'f4','dimT')
    LHF_ME_ncdf       = f.createVariable('LHF_ME_ncdf', 'f4','dimT')
    LHF_SE_ncdf       = f.createVariable('LHF_SE_ncdf', 'f4','dimT')
    SHF_RU_ncdf       = f.createVariable('SHF_RU_ncdf', 'f4','dimT')
    SHF_ME_ncdf       = f.createVariable('SHF_ME_ncdf', 'f4','dimT')
    SHF_SE_ncdf       = f.createVariable('SHF_SE_ncdf', 'f4','dimT')
    ERR_LHF_RU_ncdf   = f.createVariable('ERR_LHF_RU_ncdf', 'f4','dimT')
    ERR_LHF_ME_ncdf   = f.createVariable('ERR_LHF_ME_ncdf', 'f4','dimT')
    ERR_LHF_SE_ncdf   = f.createVariable('ERR_LHF_SE_ncdf', 'f4','dimT')
    ERR_SHF_RU_ncdf   = f.createVariable('ERR_SHF_RU_ncdf', 'f4','dimT')
    ERR_SHF_ME_ncdf   = f.createVariable('ERR_SHF_ME_ncdf', 'f4','dimT')
    ERR_SHF_SE_ncdf   = f.createVariable('ERR_SHF_SE_ncdf', 'f4','dimT')
    SHF_MEAN_ncdf     = f.createVariable('SHF_MEAN_ncdf', 'f4','dimT')
    LHF_MEAN_ncdf     = f.createVariable('LHF_MEAN_ncdf', 'f4','dimT')
    ERR_SHF_MEAN_ncdf = f.createVariable('ERR_SHF_MEAN_ncdf', 'f4','dimT')
    ERR_LHF_MEAN_ncdf = f.createVariable('ERR_LHF_MEAN_ncdf', 'f4','dimT')
    time_ncdf         = f.createVariable(varname='time', dimensions=('dimT',),datatype='float64')    
   
    # passing data into the variables
    SHF_MEAN_ncdf[:]     = SHF_mean
    LHF_MEAN_ncdf[:]     = LHF_mean
    ERR_SHF_MEAN_ncdf[:] = SHF_std
    ERR_LHF_MEAN_ncdf[:] = LHF_std
    
    LHF_RU_ncdf[:]       = LHF_RU
    LHF_ME_ncdf[:]       = LHF_ME
    LHF_SE_ncdf[:]       = LHF_SE
    SHF_RU_ncdf[:]       = SHF_RU
    SHF_ME_ncdf[:]       = SHF_ME
    SHF_SE_ncdf[:]       = SHF_SE
    ERR_LHF_RU_ncdf[:]   = Err_LHF_RU
    ERR_LHF_ME_ncdf[:]   = Err_LHF_ME
    ERR_LHF_SE_ncdf[:]   = Err_LHF_SE
    ERR_SHF_RU_ncdf[:]   = Err_SHF_RU
    ERR_SHF_ME_ncdf[:]   = Err_SHF_ME
    ERR_SHF_SE_ncdf[:]   = Err_SHF_SE

    time_ncdf[:]         = nc4.date2num(datetimeDay, units=units, calendar=calendar)
    time_ncdf.units      = units
    
    #Add global attributes
    f.description        = "surface fluxes extracted from EC stations around joyce"
    f.history            = "Created by Claudia Acquistapace cacquist@meteo.uni-koeln.de - University of Cologne"
    
    #Add local attributes to variable instances
    
    
        
    # closing ncdf file
    f.close()
    print('file saved for day '+date)
