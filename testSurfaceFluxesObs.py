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
#'20130503','20130504', '20130505','20130506','20130509','20130510',
#'20130424', '20130425', '20130427', '20130429','20130501','20130502', '20130503' ]

# new days : 20130414, 20130420, 20130426 20130428 20130430 20130524 20130525 20130527 20130528
#dateArr          = ['20130414', '20130420', '20130426', '20130428', '20130430', '20130524', '20130525', '20130527', '20130528']
dateArr          = ['20130414','20130420', \
 '20130424', '20130425', '20130426','20130427', '20130428', '20130429', '20130430', \
 '20130501', '20130502', '20130503','20130504','20130505','20130506','20130509', \
 '20130510', '20130518','20130524', '20130525','20130527','20130528']
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
        data_sheet
        SHF              = data_sheet.col_values(37)
        LHF              = data_sheet.col_values(39)
        SHF_flag         = data_sheet.col_values(43)
        LHF_flag         = data_sheet.col_values(45)
        SHF_relErr       = data_sheet.col_values(54)
        LHF_relErr       = data_sheet.col_values(55)
        WVdensity        = data_sheet.col_values(10)
        P                = data_sheet.col_values(11)
        T                = data_sheet.col_values(9)


        # removing first element of the line with the name variable
        SHF.remove(SHF[0])
        LHF.remove(LHF[0])
        SHF_flag.remove(SHF_flag[0])
        LHF_flag.remove(LHF_flag[0])
        SHF_relErr.remove(SHF_relErr[0])
        LHF_relErr.remove(LHF_relErr[0])
        WVdensity.remove(WVdensity[0])
        P.remove(P[0])
        T.remove(T[0])


        # converting lists to arrays for algebric operations (RH calculation)
#        SHF = np.asarray(SHF, dtype=np.float64)
#        LHF = np.asarray(LHF, dtype=np.float64)
#        SHF_flag = np.asarray(SHF_flag, dtype=np.float64)
#        LHF_flag = np.asarray(LHF_flag, dtype=np.float64)
#        SHF_relErr = np.asarray(SHF_relErr, dtype=np.float64)
#        LHF_relErr = np.asarray(LHF_relErr, dtype=np.float64)
        WVdensity = np.asarray(WVdensity, dtype=np.float64)  # [g/m^3]
        P = np.asarray(P, dtype=np.float64)  # [hpa]
        T = np.asarray(T, dtype=np.float64)  # [C]
        T = T + 273.15  # [K]


        # calculating relative humidity based on the formula  RH = 100% P a/(ro my Esat(Ta))
        # with P air pressure, ro air density, my = 622 the mol mass ratio between vapor and Air, Esat(Ta) the vapor
        # saturation pressure at air temperature Ta
        my = 0.622
        T0 = 273.  # [K]
        Psat0 = 6.11  # [hPa] saturation pressure at T=T0
        Rd = 287.0  # [J/Kg K]
        L = 2.50 * 10 ** 6  # [J/Kg]
        Rv = 461.5  # [J/Kg K]
        RH = np.zeros((len(T)))
        RH.fill(np.nan)
        for indTime in range(len(T)):
            # print(T[indTime])
            # print(P[indTime])
            rho_air = P[indTime] * 100. / (Rd * T[indTime])  # calculation of density of dry air
            # print(rho_air)
            Psat = Psat0 * np.exp((L / Rv) * ((1 / T0) - (1 / T[indTime])))
            # print(Psat)
            RH[indTime] = (P[indTime] * WVdensity[indTime] * 0.001) / (rho_air * my * Psat)
            #print(RH)



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


        Err_LHF          = np.asarray(LHF[:]) * np.asarray(LHF_relErr[:])
        Err_SHF          = np.asarray(SHF[:]) * np.asarray(SHF_relErr[:])

        # selecting only time stamps corresponding to the selected day
        LHF_dataframe        = pd.Series(LHF[:], index=dateTime_fluxes)
        SHF_dataframe        = pd.Series(SHF[:], index=dateTime_fluxes)
        LHF_Err_dataframe    = pd.Series(Err_LHF, index=dateTime_fluxes)
        SHF_Err_dataframe    = pd.Series(Err_SHF, index=dateTime_fluxes)        
        P_dataframe          = pd.Series(P, index=dateTime_fluxes)
        T_dataframe          = pd.Series(T, index=dateTime_fluxes)
        RH_dataframe         = pd.Series(RH, index=dateTime_fluxes)


        mask_t = (LHF_dataframe.index >= datetimeDay[0]) * (LHF_dataframe.index <= datetimeDay[-1])
        LHF_day = LHF_dataframe[mask_t].values
        LHF_Err_day = LHF_Err_dataframe[mask_t].values
        
        mask_t = (SHF_dataframe.index >= datetimeDay[0]) * (SHF_dataframe.index <= datetimeDay[-1])
        SHF_day = SHF_dataframe[mask_t].values
        SHF_Err_day = SHF_Err_dataframe[mask_t].values

        mask_t = (P_dataframe.index >= datetimeDay[0]) * (P_dataframe.index <= datetimeDay[-1])
        P_day = P_dataframe[mask_t].values

        mask_t = (T_dataframe.index >= datetimeDay[0]) * (T_dataframe.index <= datetimeDay[-1])
        T_day = T_dataframe[mask_t].values
        mask_t = (P_dataframe.index >= datetimeDay[0]) * (P_dataframe.index <= datetimeDay[-1])
        P_day = P_dataframe[mask_t].values
        mask_t = (RH_dataframe.index >= datetimeDay[0]) * (RH_dataframe.index <= datetimeDay[-1])
        RH_day = RH_dataframe[mask_t].values


        if stringStat == 'RU':
            SHF_RU              = SHF_day
            Err_SHF_RU          = SHF_Err_day
            LHF_RU              = LHF_day
            Err_LHF_RU          = LHF_Err_day
            P_RU                = P_day
            T_RU                = T_day
            RH_RU               = RH_day
            string_RU           = stringStat
    
        if stringStat == 'ME':
            SHF_ME              = SHF_day
            Err_SHF_ME          = SHF_Err_day
            LHF_ME              = LHF_day
            Err_LHF_ME          = LHF_Err_day
            P_ME                = P_day
            T_ME                = T_day
            RH_ME               = RH_day
            string_ME           = stringStat
    
        if stringStat == 'SE':
            SHF_SE              = SHF_day
            Err_SHF_SE          = SHF_Err_day
            LHF_SE              = LHF_day
            Err_LHF_SE          = LHF_Err_day
            P_SE                = P_day
            T_SE                = T_day
            RH_SE               = RH_day
            string_SE           = stringStat 
            
    
    # calculating mean fluxes over the three stations for the selected day 
    SHF_mean = []
    LHF_mean = [] 
    SHF_std  = []
    LHF_std  = []  
    P_mean   = []
    T_mean   = []
    RH_mean  = []
    for indTime in range(len(datetimeDay)):
        SHF_mean.append(np.nanmean([SHF_SE[indTime], SHF_ME[indTime], SHF_RU[indTime]]))
        LHF_mean.append(np.nanmean([LHF_SE[indTime], LHF_ME[indTime], LHF_RU[indTime]]))
        SHF_std.append(np.nanstd([SHF_SE[indTime], SHF_ME[indTime], SHF_RU[indTime]]))
        LHF_std.append(np.nanstd([LHF_SE[indTime], LHF_ME[indTime], LHF_RU[indTime]]))
        P_mean.append(np.nanmean([P_SE[indTime], P_ME[indTime], P_RU[indTime]]))
        T_mean.append(np.nanmean([T_SE[indTime], T_ME[indTime], T_RU[indTime]]))
        RH_mean.append(np.nanmean([RH_SE[indTime], RH_ME[indTime], RH_RU[indTime]]))

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
    P_MEAN_ncdf       = f.createVariable('P_MEAN_ncdf', 'f4','dimT')
    T_MEAN_ncdf       = f.createVariable('T_MEAN_ncdf', 'f4','dimT')
    RH_MEAN_ncdf       = f.createVariable('RH_MEAN_ncdf', 'f4','dimT')

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
    P_MEAN_ncdf[:]       = P_mean
    T_MEAN_ncdf[:]       = T_mean
    RH_MEAN_ncdf[:]       = RH_mean

    time_ncdf[:]         = nc4.date2num(datetimeDay, units=units, calendar=calendar)
    time_ncdf.units      = units
    
    #Add global attributes
    f.description        = "surface fluxes and P, T, RH time series extracted from EC stations around joyce"
    f.history            = "Created by Claudia Acquistapace cacquist@meteo.uni-koeln.de - University of Cologne"
    
    #Add local attributes to variable instances
    
    
        
    # closing ncdf file
    f.close()
    print('file saved for day '+date)
