#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 16:18:01 2019

@author: cacquist
"""


# coding: utf-8

# In[1]:


# ------------------------------------------------------------------------
# date      : 12.04.2018
# author    : Claudia Acquistapace
# goal      : routine to read 1D meteogram for a given date and site ( Joyce ) and extract data for the site and also level2 variables for the site Store them in a ncdf file to be copied on ostro for comparison 1to1 with observations from the ground
# DAYS WITH BOUNDARY LAYER CLOUDS OF INTEREST:
# - 20130502 (folder 20130502-default )
# - 20130505 (folder 20130505-default-redone_v1) 
# - 20130511 (folder 20130511-default )
# - 20160603 (folder 20160603-default-redone_v2 )

# ------------------------------------------------------------------------


# In[1]:


# ---- importing libraries
import numpy as np
import matplotlib
import scipy
import numpy.ma as ma
import pandas as pd
import netCDF4 as nc4
import glob
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
from myFunctions import f_calcWindSpeed_Dir


def f_processModelOutput(path_icon, \
                         iconFilename, \
                         modelInputParameters, \
                         date, \
                         cloudTimeArray, \
                         debuggingFlag, \
                         verboseFlag, \
                         pathDebugFig, \
                         pathOut, \
                         domSel):
   
    print('processing meteograms for the '+date)
    
    
    
    # ---- reading datafile selected
    data          = Dataset(path_icon+iconFilename, mode='r')
    time          = data.variables['time'][:].copy()
    datetime_ICON = nc4.num2date(data.variables['time'][:],data.variables['time'].units)
    Qi            = data.variables['QI'][:].copy()
    Qc            = data.variables['QC'][:].copy()
    T             = data.variables['T'][:].copy() # in [K]
    zonalWind     = data.variables['U'][:].copy()
    merWind       = data.variables['V'][:].copy()
    vertWind      = data.variables['W'][:].copy()
    LWP           = data.variables['TQC'][:].copy()
    IWV           = data.variables['TQV'][:].copy()
    thetaV        = data.variables['THETAV'][:].copy()
    height        = data.variables['height'][:].copy()
    P             = data.variables['P'][:].copy()        # [Pa]
    RH            = data.variables['REL_HUM'][:].copy()
    q             = data.variables['QV_DIA'][:].copy()   # [kg/kg]
    Hsurf         = float(data.station.split('_hsurf=')[-1].split('\n')[0])
    height2       = data.variables['height_2'][:].copy()
    rho           = data.variables['RHO'][:].copy()  
    SWSurfFlux    = data.variables['SOBS'][:].copy() # shortwave net flux at surface
    LWSurfFlux    = data.variables['THBS'][:].copy() # longwave net flux at surface
    LHFL          = data.variables['LHFL'][:].copy() # latent heat flux (surface)
    SHFL          = data.variables['SHFL'][:].copy() # sensible heat flux (surface)
    TempSurf      = data.variables['T_S'][:]

    print(Hsurf)
    print(height2[-1])
    print(height2[0])
    print(len(height2))
    # subtracting from model height arrays the height of the ground level at JOYCE
    # and make it comparable with the observations
    height2 = height2 - np.repeat(Hsurf, len(height2))
    height  = height -np.repeat(Hsurf, len(height))
    
    
    # --- reading dimension of height and time arrays
    dimTime       = len(datetime_ICON)
    dimHeight     = len(height2)
    
    if verboseFlag == 1: 
        print('variable extracted from the data')
        print('data loaded for '+date)
        print('dimension for height_2 :', dimHeight)
        print('dimension for time :', dimTime)
    
    # ------------------------------------------------------------------                
    # plot meridional and zonal wind for checking fields
    # ------------------------------------------------------------------
    if debuggingFlag == 1:
        if verboseFlag == 1: 
            print('no plots')


    
# =============================================================================
#         fig, ax = plt.subplots(figsize=(14,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("Fluxes at the surface [W/m2]", fontsize=16)
#         plt.plot(datetime_ICON, SWSurfFlux, label='Shortwave net flux')
#         plt.plot(datetime_ICON, LWSurfFlux, label='Longwave net flux')
#         plt.legend()
#         plt.savefig(pathDebugFig+'surface_LWSW_surfFlux_iconlem_'+date+'.png', format='png')
# 
# =============================================================================

# =============================================================================
#         fig, ax = plt.subplots(figsize=(14,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("Latent/sensible heat fluxes at the surface [W/m2]", fontsize=16)
#         plt.plot(datetime_ICON, LHFL, label='Latent heat flux')
#         plt.plot(datetime_ICON, SHFL, label='Sensible heat flux')
#         plt.legend()
#         plt.savefig(pathDebugFig+'surface_LatentSensible_heatFlux_iconlem_'+date+'.png', format='png')
# =============================================================================
        if verboseFlag == 1: 
            print('end of plotting graphs for debugging in debugging mode')         
            
            
    # ------------------------------------------------------------------                
    # defining constants needed for calculations
    # ------------------------------------------------------------------ 
    Rw    = 462.
    Rl    = 287.
    g     = 9.81
    P_0   = 100*1000.
    const = 0.286 # R/Cp
    P_0   = 100*1000.
    const = 0.286 # R/Cp
    Lv    = 2260 # J / kg
    Cp    = 1005.7 # /K Kg
    
    # ------------------------------------------------------------------                
    # derivation of water vapor mixing ratio
    # ------------------------------------------------------------------                
    r     = np.zeros((dimTime, dimHeight))
    for itempo in range(dimTime):
        for ih in range(dimHeight):
            r[itempo,ih] = q[itempo,ih]/(1. - q[itempo,ih] )
    if verboseFlag == 1: 
        print('water vapor mixing ratio calculated')         
    
    
    # ------------------------------------------------------------------
    # --- calculating cloud mask for ice and liquid clouds using thresholds on Qi, Qc
    # ------------------------------------------------------------------
    QcThreshold            = modelInputParameters['QcThresholdVar']
    QiThreshold            = modelInputParameters['QiThresholdVar']
    cloudMask              = f_cloudmask(time,height2,Qc,Qi,QiThreshold,QcThreshold)
    if (cloudTimeArray[2] == 'minmax'):
        result                 = f_calcCloudBaseTopPBLclouds(cloudMask, len(datetime_ICON), \
                                                    len(height2), height2, cloudTimeArray, datetime_ICON)
    else:
        result                 = f_calcCloudBaseTopPBLcloudsV2(cloudMask, len(datetime_ICON), \
                                                    len(height2), height2, cloudTimeArray, datetime_ICON)
    # return (CBarray, CTarray, NlayersArray, CB_collective, CT_collective)
    #CBMatrix_ICON          = result[0]
    #CTMatrix_ICON          = result[1]
    CT_array_ICON          = np.empty(len(datetime_ICON))
    CB_array_ICON          = np.empty(len(datetime_ICON))
    CT_array_ICON[:]       = np.nan
    CB_array_ICON[:]       = np.nan
    NcloudLayers           = result[2]
    CB_array_ICON[:]       = result[5]
    CT_array_ICON[:]       = result[6]
    CT_PBL          = np.empty(len(datetime_ICON))
    CB_PBL          = np.empty(len(datetime_ICON))
    CT_PBL[:]       = np.nan
    CB_PBL[:]       = np.nan
    CB_PBL[:]       = result[3]
    CT_PBL[:]       = result[4]
    print(CB_PBL)
    print(CT_PBL)
    print('*************************************************')
# =============================================================================
#     
#     for indT in range(len(datetime_ICON)):#
#         if (~np.isnan(CBMatrix_ICON[indT,0]) == True) and (~np.isnan(CTMatrix_ICON[indT,0])== True):
#             
#             indCB          = f_closest(height, CBMatrix_ICON[indT,0])
#             indCT          = f_closest(height, CTMatrix_ICON[indT,0])
#             
#             if (indCB == 0) or (indCT == 0):
#                 CT_array_ICON[indT] = np.nan
#                 CB_array_ICON[indT] = np.nan
#             else:
#                 CT_array_ICON[indT] = height[indCT]                                 # saving cloud top height
#                 CB_array_ICON[indT] = height[indCB]                                 # saving cloud base height
#                 
# =============================================================================
    print('cloud base and cloud top for ICON-LEM calculated ')
        
    
    # ------------------------------------------------------------------                
    # ---- calculating potential temperature and equivalent potential temperature
    # ------------------------------------------------------------------
    theta   = np.zeros((dimTime, dimHeight))
    theta_e = np.zeros((dimTime, dimHeight))
    theta.fill(np.nan)
    theta_e.fill(np.nan)
    
    for iTime in range(dimTime):
        for iHeight in range(dimHeight):
            if height[iHeight] < Hsurf:
                theta[iTime, iHeight] = 0.
            else:
                theta[iTime, iHeight] = T[iTime, iHeight] * (float(P_0)/float(P[iTime, iHeight]))**(const)  
    if verboseFlag == 1: 
        print('potential temperature calculated')   
    
    for iTime in range(dimTime):
        for iHeight in range(dimHeight):
            lv = (2500.-2.42*(T[iTime, iHeight]-273.15))*1000. # latent heat of vaporization in J/kg
            theta_e[iTime, iHeight] = theta[iTime, iHeight]+(lv*r[iTime, iHeight]/Cp)*         (np.power(100000./P[iTime, iHeight], Rl/Cp)) # equivalent potential temperature in K
    if verboseFlag == 1: 
        print('equivalent potential temperature calculated')   
    
    
    # ------------------------------------------------------------------
    # --- Calculating Boundary layer height using the richardson number derivation according to Seidel Et al, 2010
    # ------------------------------------------------------------------
    device = 'mod'
    PBLHeightArrRN    = f_calcPblHeightRN(thetaV,zonalWind,merWind,height2,time, device)
    if verboseFlag == 1: 
        print('height of the PBL (RN) calculated')

    # ------------------------------------------------------------------
    # --- calculation of the variance, std, skewness of the vertical velocity using a running mean window 
    # ------------------------------------------------------------------
    timeWindowSk    = modelInputParameters['timeWindowSkVar']
    runningWindow   = modelInputParameters['runningWindowVar']
    resultDyn       = f_runningMeanSkewnessVarianceStd_W(time, timeWindowSk, runningWindow, height2, vertWind)
    # output of the function : varianceW, stdWmatrix, SKmatrix
    varianceWmatrix = resultDyn[0]
    stdWmatrix      = resultDyn[1]
    SKmatrix        = resultDyn[2]
        
    if verboseFlag == 1: 
        print('variance, std and skewness of w calculated')
    
    print('std max = '+str(np.nanmax(stdWmatrix)))
    # ------------------------------------------------------------------
    # --- Calculating Boundary layer height using the threshold on variance of w ()
    # ------------------------------------------------------------------
    device         = 'mod'
    sigmaW         = stdWmatrix
    sigmaThreshold = modelInputParameters['SigmaWThresStd'] #  m/s, threshold for std of w from Schween et al, 2014.AMT
    PBLHeightArrTW = f_calcPblHeightTW(sigmaW,sigmaThreshold,height2,time, device)
    if verboseFlag == 1: 
        print('height of the PBL (TW) calculated')
        

        
    # ------------------------------------------------------------------
    # --- Calculating variance over the timewindow using running mean
    # ------------------------------------------------------------------
    #timewindow      = modelInputParameters['timewindowVar']
    #varianceWmatrix = f_calcWvariance(vertWind,time,height2,timewindow)
    #if verboseFlag == 1: 
    #    print('variance of vertical velocity calculated')  

    # ------------------------------------------------------------------
    # --- calculation of the connection of the turbulence to the surface. 
    # ------------------------------------------------------------------
    #Turbulence is connected to the surface if checks if variance at 200 m of height is greater than 0.03 for turbulence
    # calculating the time serie of difference of the sigmaW and the threshold value at 200 m height
    deltaSigma         = np.subtract(varianceWmatrix, 0.03)[:,f_closest(height,200.)]
    connection2Surface = []                    # array indicating connection of the turbulence to the surface
    # calculating connection to the surface. =0 ( not connected, if sigmaW(200)-sigmaGround)<0,
    # =1 (connected thus turbulent, if sigmaW(200)-sigmaGround)>0)
    for itime in range(dimTime):
        if deltaSigma[itime] < 0.:
            connection2Surface.append(0)
        else:
            connection2Surface.append(1)
    if verboseFlag == 1: 
        print('connection of turbulence with the surface calculated')        
            
    # ------------------------------------------------------------------        
    #  ---- calculation of the stability array
    # ------------------------------------------------------------------
    stabilityArr = []
    # difference of temperature between 150m and closest level to surface
    deltaT       = np.subtract(T, T[f_closest(height,Hsurf),:])[:,f_closest(height,150.)]
    for itime in range(dimTime):
        
        #print(Tarray[indRef]-Tarray[indGround])
        if deltaT[itime] < 0.3:
            stabilityArr.append(1)
        else:
            stabilityArr.append(0)
    if verboseFlag == 1: 
        print('stability at the surface calculated')        



    # ------------------------------------------------------------------
    # --- Calculation of wind shear as done for PBL  ( running mean over 30 min of sqrt(Delta U^2 + delta V^2))/delta H 
    # where variations are calculated over 5 range gates 
    # ------------------------------------------------------------------
    windData      = f_calcWindSpeed_Dir(datetime_ICON, height2, zonalWind, merWind)
    windSpeed     = windData['windSpeed']
    windDirection = windData['windDirection']
    

# =============================================================================
# --- calculating shear of horizontal wind 
    u_rm = np.zeros((len(datetime_ICON), len(height2)))
    v_rm = np.zeros((len(datetime_ICON), len(height2)))
# --- defining running mean values of zonal and meridional wind
    for indH in range(0,len(height2)):
         zonal = pd.Series(zonalWind[:,indH])
         mer   = pd.Series(merWind[:,indH])
         #u_rm[:,indH] = pd.rolling_mean(zonalWind[:,indH], window=200) 
         #v_rm[:,indH] = pd.rolling_mean(merWind[:,indH], window=200) 
         u_rm[:,indH] = zonal.rolling(200).mean()
         v_rm[:,indH] = mer.rolling(200).mean()
#         
# calculating wind shear and horizontal wind
    shear_ICON = np.zeros((len(datetime_ICON), len(height2)))
    for indT in range(0,len(datetime_ICON)):
        for indH in range(0,len(height2)):

             if (indH < 2.) or (indH > len(height2)-3):
                 shear_ICON[indT, indH] = 0.
             else:
                 deltaV = (np.absolute(v_rm[indT, indH+2] - v_rm[indT, indH-2]))**2
                 deltaU = (np.absolute(u_rm[indT, indH+2] - u_rm[indT, indH-2]))**2
                 deltaH = np.absolute(height[indH+2] - height[indH-2])
                 shear_ICON[indT, indH] = (np.sqrt(deltaU + deltaV))/deltaH
# =============================================================================
    if verboseFlag == 1:
        print('horizontal wind speed, direction and shear calculated')
    # ------------------------------------------------------------------        
    # ----calculating boundary layer classification (version from submitted paper)
    # ------------------------------------------------------------------
    ylim        = np.repeat(3000, dimTime)     # defining array of heights up to which PBL classification is calculated
    gradWindThr = 0.01
    SigmaWThres = 0.2
    outputClass = f_PBLClass(datetime_ICON, \
                             height2, \
                             gradWindThr, \
                             SigmaWThres, \
                             ylim, \
                             cloudMask, \
                             varianceWmatrix, \
                             SKmatrix, \
                             stabilityArr, \
                             connection2Surface, \
                             shear_ICON, \
                             CB_array_ICON)
    PBLclass     = outputClass[0]
    
    if verboseFlag == 1: 
        print('PBL classification calculated')

        
    if debuggingFlag == 1:
        print('dimensions of PBL class')
        print(np.shape(PBLclass))
        
# =============================================================================
#         # plotting classification 
#         fig, ax = plt.subplots(figsize=(10,4))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         cax = ax.pcolormesh(datetime_ICON, height2, PBLclass.transpose(), vmin=0., vmax=6., cmap=plt.cm.get_cmap("jet", 7))
#         ax.set_ylim(Hsurf,3000.)                                               # limits of the y-axes
#         #ax.set_xlim(0,24)                                 # limits of the x-axes
#         ax.set_title("PBL classification", fontsize=14)
#         ax.set_xlabel("time [UTC]", fontsize=12)
#         ax.set_ylabel("height [m]", fontsize=12)
#         cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3, 4, 5, 6], orientation='vertical')
#         cbar.ticks=([0,1,2,3,4,5,6])
#         cbar.ax.set_yticklabels(['no class','in cloud','non turb','cloud driven','convective', 'intermittent','wind shear'])
#         cbar.set_label(label="PBL type",size=12)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.aspect=20
#         plt.savefig(pathDebugFig+'PBLclassification_iconlem_'+date+'.png', format='png')
# =============================================================================
    
    

    
    # ------------------------------------------------------------------
    # --- calculation of the LCL 
    # ------------------------------------------------------------------
    # determining P, T and RH at the surface 
    Psurf    =  data.variables['P_SFC'][:].copy()
    Tsurf    = data.variables['T2M'][:].copy()
    RHsurf   = RH[:,149]
    LCLarray = []
    from myFunctions import lcl
    for iTime in range(dimTime):
        LCLarray.append(lcl(Psurf[iTime],Tsurf[iTime],RHsurf[iTime]/100.))
    if verboseFlag == 1: 
        print('LCL calculated')   
    
    
    
    
    # ------------------------------------------------------------------
    # calculate LTS index for lower tropospheric stability (Wood and Bretherton, 2006)
    # ------------------------------------------------------------------
    LTS  = np.zeros(dimTime)
    H700 = np.zeros(dimTime)
    Pthr = 700 * 100. # Pressure level of 700 Hpa used as a reference
    # calculating height of the surface
     
    indSurf = 146# f_closest(height,Hsurf)
    for iTime in range(dimTime):
        indP700 = f_closest(P[iTime,:],Pthr)
        LTS[iTime] = theta[iTime, indP700] - theta[iTime, indSurf]
        H700[iTime] = height[indP700]
    if verboseFlag == 1: 
        print('LTS calculated')
    #print(theta[4500, indP700])
    #print(theta[4500, indSurf])
    #print(theta[4500, :])
    
        
    # ------------------------------------------------------------------                
    # ---- calculating liquid potential temperature
    # ------------------------------------------------------------------
    theta_liquid = np.zeros((dimTime, dimHeight))
    theta_liquid.fill(np.nan)
    
    for iTime in range(dimTime):
        for iHeight in range(dimHeight):
            if height[iHeight] < Hsurf:
                theta_liquid[iTime, iHeight] = 0.
            else:
                theta_liquid[iTime, iHeight] = theta[iTime, iHeight] - (Lv/Cp)* Qc[iTime, iHeight] 
                
    if verboseFlag == 1: 
        print('liquid potential temperature calculated')   
    
    


    
    
    # ------------------------------------------------------------------
    # ------- saving mean outputs as ncdf file 
    # ------------------------------------------------------------------
    
    f = nc4.Dataset(pathOut+'icon_lem_derivedproperties'+date+'.nc','w', format='NETCDF4')      # creates a netCDF file for writing 
    tempgrp = f.createGroup('Temp_data')                                                    # creates a group: A netCDF group is basically a directory or folder within the netCDF dataset
    
    # specify dimensions of the data ( each dimension of multidimensiona array needs to be given a name and a length)
    tempgrp.createDimension('dimH', len(height2))                                 # dimension for height
    tempgrp.createDimension('dimHlong', len(height))                                 # dimension for height
    
    tempgrp.createDimension('dimHsurf', 1)                                         # dimension for scalar values
    tempgrp.createDimension('dimT', len(datetime_ICON))                                   # dimension for time
    tempgrp.createDimension('NclassesPBL', 8)                # dimension for the number of cloud layers found
    tempgrp.createDimension('dimHlarger', len(height))                                 # dimension for height
    tempgrp.createDimension('nchar', 5)
    
    # preallocating netCDF variables for data storage
    varHeight2             = tempgrp.createVariable('height2', 'f4', 'dimH')
    varHeight             = tempgrp.createVariable('height', 'f4', 'dimHlong')
    vardomain             = tempgrp.createVariable('domain', 'S1', 'nchar')
    vartime               = tempgrp.createVariable('datetime_ICON', 'f4', 'dimT')
    varLTS                = tempgrp.createVariable('LTS', 'f4', 'dimT')
    varPBLheight          = tempgrp.createVariable('PBLHeightArrRN', 'f4', 'dimT')
    varPBLheight2         = tempgrp.createVariable('PBLHeightArrTW', 'f4', 'dimT')
    varCBheight           = tempgrp.createVariable('cloudBase', 'f4', 'dimT')
    varCTheight           = tempgrp.createVariable('cloudTop', 'f4', 'dimT')
    varCloudLayers        = tempgrp.createVariable('NcloudLayers', 'f4', 'dimT')
    varCB_PBL             = tempgrp.createVariable('cloudBasePBL', 'f4', 'dimT')
    varCT_PBL             = tempgrp.createVariable('cloudTopPBL', 'f4', 'dimT')
    varHsurf              = tempgrp.createVariable('HeightSurface', 'f4', 'dimHsurf')
    varLCL                = tempgrp.createVariable('LCLarray', 'f4', 'dimT')
    varLWP                = tempgrp.createVariable('LWP', 'f4', 'dimT')
    varIWV                = tempgrp.createVariable('IWV', 'f4', 'dimT')
    varLHFL               = tempgrp.createVariable('LHFL', 'f4', 'dimT')
    varSHFL               = tempgrp.createVariable('SHFL', 'f4', 'dimT')
    varLWSF               = tempgrp.createVariable('LWSurfFlux', 'f4', 'dimT')
    varSWSF               = tempgrp.createVariable('SWSurfFlux', 'f4', 'dimT')    
    # PBL class and connected flags, LTS clouds, SW clouds, PBL height, CB height
    varPBL_class          = tempgrp.createVariable('PBLclass', 'f4', ('dimT','dimH'))
    varflagCloud          = tempgrp.createVariable('flagCloud', 'f4', ('dimT','dimH'))
    varQc                 = tempgrp.createVariable('Qc', 'f4', ('dimT','dimH'))
    varQi                 = tempgrp.createVariable('Qi', 'f4', ('dimT','dimH'))
    varflagTurb           = tempgrp.createVariable('flagTurb', 'f4', ('dimT','dimH'))
    varflagcloudDriven    = tempgrp.createVariable('flagcloudDriven', 'f4', ('dimT','dimH'))
    varflagInstability    = tempgrp.createVariable('flagInstability', 'f4',('dimT','dimH'))
    varflagWindShear      = tempgrp.createVariable('flagWindShear', 'f4', ('dimT','dimH'))
    varflagSurfDriven     = tempgrp.createVariable('flagSurfaceDriven', 'f4', ('dimT','dimH'))
    varvarianceW          = tempgrp.createVariable('varianceW', 'f4', ('dimT','dimH'))
    varHwind              = tempgrp.createVariable('windSpeed', 'f4', ('dimT','dimH'))
    varWindDirection      = tempgrp.createVariable('windDirection', 'f4', ('dimT','dimH'))
    varShearHwind         = tempgrp.createVariable('shearHwind', 'f4', ('dimT','dimH'))
    varcloudMask          = tempgrp.createVariable('cloudMask', 'f4', ('dimT','dimH'))
    varthetaPot           = tempgrp.createVariable('theta', 'f4', ('dimT','dimH'))
    varskewnessW          = tempgrp.createVariable('skewnessW', 'f4', ('dimT','dimH'))
    varstdWmatrix         = tempgrp.createVariable('stdWmatrix', 'f4', ('dimT','dimH'))
    varMixingRatio        = tempgrp.createVariable('r', 'f4', ('dimT','dimH'))
    varthetaL             = tempgrp.createVariable('theta_liquid', 'f4', ('dimT','dimH'))
    varthetaPot_e         = tempgrp.createVariable('theta_e', 'f4', ('dimT','dimH'))
    varw                  = tempgrp.createVariable('vertWind', 'f4', ('dimT','dimHlarger'))
    varP                  = tempgrp.createVariable('P', 'f4', ('dimT','dimH'))
    varRH                 = tempgrp.createVariable('RH', 'f4', ('dimT','dimH'))
    varQ                  = tempgrp.createVariable('q', 'f4', ('dimT','dimH'))
    varT                  = tempgrp.createVariable('T', 'f4', ('dimT','dimH'))
    varMerWind            = tempgrp.createVariable('merWind', 'f4', ('dimT','dimH'))
    varZonWind            = tempgrp.createVariable('zonalWind', 'f4', ('dimT','dimH'))
    varRho                = tempgrp.createVariable('rho', 'f4', ('dimT','dimH'))
    varT_surf             = tempgrp.createVariable('TempSurf', 'f4', 'dimT') 
    
    # passing data into the variables
    varHeight2[:]         = height2
    varHeight[:]          = height
    vardomain             = domSel
    vartime[:]            = time
    varLTS[:]             = LTS
    varPBLheight[:]       = PBLHeightArrRN
    varPBLheight2[:]      = PBLHeightArrTW
    varCBheight[:]        = CB_array_ICON
    varCTheight[:]        = CT_array_ICON
    varCloudLayers[:]     = NcloudLayers
    varCB_PBL[:]          = CB_PBL
    varCT_PBL[:]          = CT_PBL
    varHsurf              = Hsurf
    varLCL[:]             = LCLarray
    varLWP[:]             = LWP
    varIWV[:]             = IWV
    varLHFL[:]            = LHFL
    varSHFL[:]            = SHFL
    varLWSF[:]            = LWSurfFlux
    varSWSF[:]            = SWSurfFlux
    varPBL_class[:,:]     = PBLclass
    varflagCloud[:]       = outputClass[1]
    varflagTurb[:]        = outputClass[2]
    varflagcloudDriven[:] = outputClass[3]
    varflagInstability[:] = outputClass[4]
    varflagWindShear[:]   = outputClass[5]
    varflagSurfDriven[:]  = outputClass[6]
    varvarianceW[:,:]     = varianceWmatrix
    varHwind[:,:]         = windSpeed
    varWindDirection[:,:] = windDirection
    varShearHwind[:,:]    = shear_ICON
    varcloudMask[:,:]     = cloudMask
    varthetaPot[:,:]      = theta
    varskewnessW[:,:]     = SKmatrix
    varstdWmatrix[:,:]    = stdWmatrix
    varMixingRatio[:,:]   = r
    varthetaL[:,:]        = theta_liquid
    varthetaPot_e[:,:]    = theta_e
    varw[:,:]             = vertWind
    varP[:,:]             = P
    varRH[:,:]            = RH
    varQ[:,:]             = q    
    varT[:,:]             = T
    varMerWind[:,:]       = merWind
    varZonWind[:,:]       = zonalWind
    varRho[:,:]           = rho
    varQc[:,:]            = Qc
    varQi[:,:]            = Qi
    varT_surf[:]          = TempSurf
    
    
    #Add global attributes
    f.description = "icon lem model derived physical quantities and PBL classification"
    f.history = "Created by Claudia Acquistapace cacquist@meteo.uni-koeln.de - University of Cologne"
    
    #Add local attributes to variable instances
    varPBL_class.units = '1=in cloud, 2=non turb, 3=cloud driven, 4=convective, 5=intermittent, 6=wind shear'
    vartime.units = 'seconds since '+date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' 00:00:00'
    
    # closing ncdf file
    f.close()
                        
    print('File Saved ')