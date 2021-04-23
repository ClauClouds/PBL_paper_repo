#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:49:53 2019

@ author: cacquist
@  date : 16 juli 2019
@  goal : built statistics of PBL observations and ICON LEM model outputs. The 
code reads in the data from observations and model from the three different sites,
check if data are not there and sets a flag for that, For each positive data 
flag, it resamples with respect to time and height grid each variable. It then,
in the end, stores all resampled data in daily ncdf files with specific 
variable types, Obs data are from:
    - tower obs
    - mwr radiometer data la, lb, lc
    - PBL classification data
    - cloudnet classification
    - radiosondes
    - gps
    - sun photo cloud mode
    
"""
# ----- importing libraries needed
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
import math

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
    
from f_processModelOutput import f_processModelOutput
from myFunctions2 import f_readingTowerData
from myFunctions import f_resamplingfield
from myFunctions2 import hourDecimal_to_datetime
from myFunctions import f_resampling_twoD_Field
from myFunctions import f_downscaleScalarfield
from myFunctions import f_downscalevectorfield
#from myFunctions import f_processRadiosondesDay
from myFunctions import f_calcThermodynamics
from cloudnetFunctions import f_calculateCloudMaskCloudnet
from myFunctions import f_resamplingMatrixCloudnet
from myFunctions import f_calculateAllCloudQuantities
from myFunctions import f_selectingPBLcloudWindow

# ----- user defined parameters for model data reading and processing
station       = 'joyce'
debuggingFlag = 0
verboseFlag   = 1
reprocICON    = 0        # flag for reprocessing the icon output data
timeSpinOver  = 0.0      # ending time of the spin up of the model (corresponds to the time at which we start to calculate cloud fraction
intTime       = 400.     # time interval over which calculate cloud fraction and coupling parameters [seconds] corresponding to minutes with a resolution of model output of 9 sec (5 min = 50), (15 min = 150)
QcThreshold   = 10**(-7) # threshold in Qc to identify liquid cloud presence
QiThreshold   = 10**(-7) # threshold in Qc to identify ice cloud presence   
SigmaWThres   = 0.4      # threshold in the standard deviation of the vertical velocity, used as a threshold for identifying turbulent regions
                         # the value is taken from the Schween et al, 2014 amt paper, 
nTimeMean     = 200      # number of time stamps to be collected to average over 30 min=1800 s
timeStep      = 33       # time step for running mean
timewindow    = 200      # time window for calculating running mean of variance corresponding to 30 min with 9 sec resolution (200*9=1800 sec = 30 min)
gradWindThr   = 0.01     # Threshold for wind gradient to determine wind shear presence in the PBLclas
timeWindowSk  = 33       # number of time stamps corresponding to 5 minutes in the ICON file
runningWindow = 200*4    # number of time stamps corresponding to 30 minutes running window for the average


# ----- creating dictionary of input parameters to process icon lem model output
modelInputParameters = {'timeSpinOverVar':timeSpinOver, 'intTimeVar':intTime, 'QcThresholdVar':QcThreshold, \
                  'QiThresholdVar':QiThreshold, 'SigmaWThresStd':SigmaWThres, 'nTimeMeanVar':nTimeMean, \
                  'timeStepVar':timeStep, 'timewindowVar':timewindow, 'gradWindThrVar':gradWindThr, \
                  'timeWindowSkVar':timeWindowSk, 'runningWindowVar':runningWindow}


# ----- define list of days to be processed 
#dayList           = ['20130420', '20130426','20130428', '20130430','20130524','20130525','20130527', '20130528']
#dayList           = ['20130424','20130425', '20130427','20130429','20130501',\
#                     '20130502','20130503','20130504', '20130505','20130506',\
#                     ,'20130518']
# days processed on 17 april: '20130426''20130414','20130420','20130424', '20130425',
# '20130427','20130428','20130429','20130430','20130501','20130502','20130503','20130504', '20130505','20130506','20130509','20130510'
#
#
# ,'20130424', '20130425', '20130427','20130429','20130501','20130502', '20130518'
dayList = ['20130414','20130420', \
 '20130424', '20130425', '20130426','20130427', '20130428', '20130429', '20130430', \
 '20130501', '20130502', '20130503','20130504','20130505','20130506','20130509', \
 '20130510', '20130518','20130524', '20130525','20130527','20130528']
# new days : = ['20130414','20130420', '20130426','20130428', '20130430','20130524','20130525','20130527', '20130528']
Ndays             = len(dayList)
print('total number of days to process : ', Ndays)
# dayListAll = ['20130413','20130414','20130417','20130418','20130419','20130420', \
#,'20130421','20130422','20130423','20130424', '20130425', '20130427', 20130429', \
# '20130501', '20130502', '20130503','20130504','20130505','20130506','20130509', \
# '20130510','20130518','20130519','20130524','20130525','20130527','20130528', \
# '20130530',]


# ----- defining input directories for data
path_tower          = '/data/hatpro/jue/hdcp2/meteo_data/'#'/data/TR32/D2/data/juelich/met_mast/'
path_cloudnet_cat   = '/data/hatpro/jue/cloudnet/juelich/processed/categorize/2013/'
path_icon           = '/data/inscape/icon/experiments/juelich/meteo-4nest/'
path_bl_class       = '/data/hatpro/jue/cloudnet/juelich/products/bl-classification/2013/'
path_windLidarCeilo = '/data/TR32/D2/data/wind_lidar/data/nc/'
path_mwr_joyce      = '/data/hatpro/hps/data/level2/'
path_mwr_kit        = '/data/hatpro/hpk/data/level2/'
path_mwr_lacross    = '/data/hatpro/hpl/data/level2/'
path_gps            = '/data/hatpro/jue/data/gps/'
path_radiation      = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/'
path_radiosondes    = '/home/juelich/rs_hope/KIT/'
path_LH_SH_fluxes   = '/work/cacquist/HDCP2_S2/statistics/dataset_obs_model/'
path_LWC            = '/data/hatpro/jue/cloudnet/juelich/products/lwc-scaled-adiabatic/'
patch               = 'patch003' # patch002, patch003, patch004
domSel              = 'DOM03'
pathOut             = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'
# ----- defining data flags
dataFlagArr         = np.repeat(1, 14)
dataFlagLabel       = ['tower', 'cloudnetClass', 'PBLclass', 'MWR_joyce','radiosoundings', \
                     'LWC_Cloudnet_prod', 'windLidar_ceilo', 'windLidarScans']


# ----- loop on the number of days
for iDay in range(Ndays):


    # set day, month and year string and output path for plots for the selected day
    date              = dayList[iDay]
    yy                = dayList[iDay][0:4]
    mm                = dayList[iDay][4:6]
    dd                = dayList[iDay][6:8]
    pathDebugFig      = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/'+date+'/'

    # set, based on human inspection, the time interval in which to select clouds
    humanInfo = f_selectingPBLcloudWindow(date)

    print('PROCESSING DAY:'+date)
    print('file '+str(iDay)+'of '+str(Ndays))
    if (reprocICON == 1) or (os.path.isfile('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/icon_lem_derivedproperties'+date+'.nc') == True):
        if verboseFlag == 1:
            print('reading icon lem data for this day already processed')
        iconLemData = Dataset('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/icon_lem_derivedproperties'+date+'.nc', mode='r')
    else:
        if verboseFlag == 1:
            print('reading meteograms')
        # ---- reading icon lem model output for the selected day
        iconFilename          = 'METEOGRAM_'+patch+'_'+date+'_'+station+'.nc'
        readingFile           = f_processModelOutput(path_icon, iconFilename, \
                                                 modelInputParameters, \
                                                 date, \
                                                 humanInfo, \
                                                 debuggingFlag, \
                                                 verboseFlag, \
                                                 pathDebugFig, \
                                                 pathOut, \
                                                 domSel)
        # assigning a string name to the file just created for the next steps of the code
        iconLemData           = Dataset('/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_'+patch+'/icon_lem_derivedproperties'+date+'.nc', mode='r')
        print('reprocessing or processing the day for the first time')


    
    # -----------------------------------------------------------------------------------    
    # ---- reading icon lem data and building dataframes for regridding the observational
    # datasets on the same time/height grid of the model.
    # -----------------------------------------------------------------------------------
    time_ICON      = iconLemData.groups['Temp_data'].variables['datetime_ICON'][:].copy()
    datetime_ICON  = nc4.num2date(iconLemData.groups['Temp_data'].variables['datetime_ICON'][:], \
                                  iconLemData.groups['Temp_data'].variables['datetime_ICON'].units) 
    height_ICON    = iconLemData.groups['Temp_data'].variables['height2'][:].copy()
    cloudMask_ICON = iconLemData.groups['Temp_data'].variables['cloudMask'][:].copy()
    Hsurf          = 97.4#iconLemData.groups['Temp_data'].variables['HeightSurface'][:].copy()

    # ---- defining ICON data as dataframe reference
    ICON_DF        = pd.DataFrame(cloudMask_ICON, index=datetime_ICON, columns=height_ICON)     
    ICON_DF_T      = pd.DataFrame(cloudMask_ICON.transpose(), index=height_ICON, columns=datetime_ICON)
    print('icon file read')

    print('dimensione time icon, datetime', len(time_ICON))


    # -----------------------------------------------------------------------------------
    # ---- LWC data from scaled adiabatic cloudnet profiles
    # -----------------------------------------------------------------------------------       
    fileLWC      = date+'_juelich_lwc-scaled-adiabatic.nc'
    if (os.path.isfile(path_LWC+str(yy)+'/'+fileLWC)):
        dataFlagArr[5] = 1
        print('LWC cloudnet scaled adiabatic data found: reading and resampling data')
        dataLWC      = Dataset(path_LWC+str(yy)+'/'+fileLWC, mode='r')
        LWC          = dataLWC.variables['lwc'][:]
        datetime_LWC = nc4.num2date(dataLWC.variables['time'][:], \
                        dataLWC.variables['time'].units)
        height_LWC   = dataLWC.variables['height'][:]
        LWC_obs_res  = f_resampling_twoD_Field(LWC, datetime_LWC, \
                                              height_LWC, ICON_DF, ICON_DF_T)
    else:
        dataFlagArr[5]    = 0
        
    # test plots

    #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
    #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    
    # -----------------------------------------------------------------------------------
    # ---- radiosoundings
    # -----------------------------------------------------------------------------------   
    pathIn   = path_radiosondes+yy+mm+dd+'/'
    fileList = glob.glob(pathIn+'*.txt')
    
    for indRemove in range(len(fileList)):
        if fileList[indRemove] == 'KIT_HOPE_2013042013.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042411.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042509.txt':
            del fileList[indRemove]            
        if fileList[indRemove] == 'KIT_HOPE_2013042511.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013042616.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013043011.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013050311.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013050515.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013052513.txt':
            del fileList[indRemove]
        if fileList[indRemove] == 'KIT_HOPE_2013052715.txt':
            del fileList[indRemove]
    
    if os.listdir(pathIn):    
        # folder full, process radiosondes that are in
        from myFunctions import f_processRadiosondesDay
        radiosondeList = f_processRadiosondesDay(fileList, yy, mm, dd)
        dataFlagArr[4]    = 1        
    else: 
        # folder empty
        dataFlagArr[4]    = 0
    
    # -----------------------------------------------------------------------------------
    # ---- tower observations 
    # -----------------------------------------------------------------------------------
    pathIn   = path_tower+'/'+yy+'/'
    #filename = 'meteoFZJ'+dayList[iDay]+'.dat'
    filename = 'sups_joy_mett00_l1_any_v00_'+date+'000000.nc'
    if (os.path.isfile(pathIn+filename)):
        print('tower data found: reading and resampling data')
        dataFlagArr[0] = 1
        tower_dict     = f_readingTowerData(date, pathIn)
        
        # ---- process tower observations
        P_surf = tower_dict['P'][:]
        T_surf = tower_dict['Tsurf'][:]
        datetime_tower = tower_dict['time'][:]

        # ---- resampling tower observations ( used only for surface values, so no resampling needed in height)
        print('resampling tower observations on ICON time resolution (surf meas, no need for height resampling')
        Psurf_res = f_resamplingfield(P_surf, datetime_tower, ICON_DF)
        Tsurf_res = f_resamplingfield(T_surf, datetime_tower, ICON_DF)
        print('resampled tower observations: Psurf, Tsurf,')

    else:
        print('tower data not found')
        dataFlagArr[0] = 0
        
        
    # -----------------------------------------------------------------------------------
    # ---- cloudnet classification
    # -----------------------------------------------------------------------------------
    pathIn   = path_cloudnet_cat
    filename = date+'_juelich_categorize.nc'
    if (os.path.isfile(pathIn+filename)):

        print('cloudnet class found: reading and resampling data')
        dataFlagArr[1] = 1
        CLOUDNET_data         = Dataset(pathIn+filename, mode='r')
    
        # ----- reading CLOUDNET data variables
        time_CLOUDNET         = CLOUDNET_data.variables['time'][:].copy()
        height_CLOUDNET       = CLOUDNET_data.variables['height'][:].copy()
        datetime_CLOUDNET     = hourDecimal_to_datetime(int(yy), int(mm), int(dd), time_CLOUDNET)
        cloudnet              = CLOUDNET_data.variables['category_bits'][:].copy()
        #Ze_CLOUDNET           = CLOUDNET_data.variables['Z'][:].copy()
        #Vd_CLOUDNET           = CLOUDNET_data.variables['v'][:].copy()
        #Sw_CLOUDNET           = CLOUDNET_data.variables['width'][:].copy()
        P_CLOUDNET            = CLOUDNET_data.variables['pressure'][:].copy()          # [Pa]
        T_CLOUDNET            = CLOUDNET_data.variables['temperature'][:].copy()       # [K]
        Q_CLOUDNET            = CLOUDNET_data.variables['specific_humidity'][:].copy() # [Kg/Kg]
        model_Height_CLOUDNET = CLOUDNET_data.variables['model_height'][:].copy()
        Z_CLOUDNET            = CLOUDNET_data.variables['Z'][:].copy()
        CB_CLOUDNET           = CLOUDNET_data.variables['Z'][:].copy()

    # ---- resampling cloudnet observations ( used only for surface values, so no resampling needed in height)
    print('resampling Pressure, temperature and humidity from COSMO on ICON time and height resolution')
    P_cosmo_res = f_resampling_twoD_Field(P_CLOUDNET, datetime_CLOUDNET, \
                                              model_Height_CLOUDNET, ICON_DF, ICON_DF_T)    
    T_cosmo_res = f_resampling_twoD_Field(T_CLOUDNET, datetime_CLOUDNET, \
                                              model_Height_CLOUDNET, ICON_DF, ICON_DF_T) 
    Q_cosmo_res = f_resampling_twoD_Field(Q_CLOUDNET, datetime_CLOUDNET, \
                                              model_Height_CLOUDNET, ICON_DF, ICON_DF_T) 
    cloudnet_res = f_resamplingMatrixCloudnet(datetime_CLOUDNET, height_CLOUDNET, \
                                              cloudnet, datetime_ICON, height_ICON, cloudMask_ICON) 
    Z_cloudnet_res = f_resamplingMatrixCloudnet(datetime_CLOUDNET, height_CLOUDNET, \
                                              Z_CLOUDNET, datetime_ICON, height_ICON, cloudMask_ICON)
     
    dictCosmo = {'pressure':P_cosmo_res.values.transpose(),
                 'temperature':T_cosmo_res.values.transpose(),
                 'absoluteHumidity':Q_cosmo_res.values.transpose(),
                 'cloudnetCategorization':cloudnet_res.data.transpose(),
                 'Reflectivity':Z_cloudnet_res.data.transpose()
                 }

    # -----------------------------------------------------------------------------------
    # ---- pbl classification
    # -----------------------------------------------------------------------------------
    pathIn   = path_bl_class
    filename = date+'_bl_classification_juelich_t_3min.nc'
    if (os.path.isfile(pathIn+filename)):    
    
         
        print('boundary layer class found: reading and resampling data')
        dataFlagArr[2] = 1
        OBS_PBL_data = Dataset(pathIn+filename, mode='r')        
        
        # ----- reading PBL class data variables
        time_PBL_class        = OBS_PBL_data.variables['time'][:].copy()
        time_dwl              = OBS_PBL_data.variables['time_dwl'][:].copy()
        datetime_dwl          = hourDecimal_to_datetime(int(yy), int(mm), int(dd), time_dwl)
        height_PBL_class      = OBS_PBL_data.variables['height'][:].copy()
        PBLclass              = OBS_PBL_data.variables['bl_classification'][:].copy()
        #PBLclass = PBLclass[np.where(time_PBL_class < 2.)[0], :]
        datetime_PBL          = hourDecimal_to_datetime(int(yy), int(mm), int(dd), time_PBL_class)
        beta                  = OBS_PBL_data.variables['beta'][:].copy()
        eps                   = OBS_PBL_data.variables['eps'][:].copy()
        shear                 = OBS_PBL_data.variables['shear'][:].copy()
        Hwind                 = OBS_PBL_data.variables['speed'][:].copy()
        wDir                  = OBS_PBL_data.variables['dir'][:].copy()


    # -----------------------------------------------------------------------------------
    # ---- wind lidar/ceilometer data
    # -----------------------------------------------------------------------------------
    pathIn   = path_windLidarCeilo+'/'+yy+'/'+mm+'/'+dd+'/'
    # filename for LH estimate based on standard deviation over 30min of vertical velocity, mean w, skewness etc.
    filename = 'sups_joy_dlidST00_l2_zmlaw_v00_'+date+'000000.nc'  
    if (os.path.isfile(pathIn+filename)):
        
         
        print('wind lidar/ceilometer data found: reading and resampling data')
        dataFlagArr[6] = 1
        WIND_CEILO_data = Dataset(pathIn+filename, mode='r')        
    
        beta                  = WIND_CEILO_data.variables['beta'][:].copy()
        w                     = WIND_CEILO_data.variables['w'][:].copy()
        sigma_w               = WIND_CEILO_data.variables['w_sdev'][:].copy()
        w_skew                = WIND_CEILO_data.variables['w_skew'][:].copy()
        zmlaw                 = WIND_CEILO_data.variables['zmlaw'][:].copy()
        zmlawErr              = WIND_CEILO_data.variables['zmlaw_uncty'][:].copy()
        datetime_Lidar        = nc4.num2date(WIND_CEILO_data.variables['time'][:], \
                                             WIND_CEILO_data.variables['time'].units)
        heightLidar           = WIND_CEILO_data.variables['height'][:].copy()
        
    # ---- resampling cloudnet observations ( used only for surface values, so no resampling needed in height)
    print('resampling PBL variables and lidar on ICON time and height resolution')
    w_obs_res        = f_resampling_twoD_Field(w, datetime_Lidar, \
                                              heightLidar, ICON_DF, ICON_DF_T)
    Hwind_obs_res    = f_resampling_twoD_Field(Hwind, datetime_dwl, \
                                             height_PBL_class, ICON_DF, ICON_DF_T)
    skew_obs_res     = f_resampling_twoD_Field(w_skew, datetime_Lidar, \
                                             heightLidar, ICON_DF, ICON_DF_T)    
    PBLclass_obs_res = f_resampling_twoD_Field(PBLclass, datetime_PBL, \
                                              height_PBL_class, ICON_DF, ICON_DF_T)     
    shear_obs_res    = f_resampling_twoD_Field(shear, datetime_PBL, \
                                              height_PBL_class, ICON_DF, ICON_DF_T)     
    wDir_obs_res     = f_resampling_twoD_Field(wDir, datetime_dwl, \
                                              height_PBL_class, ICON_DF, ICON_DF_T)  
    zmlaw_obs_res    = f_resamplingfield(zmlaw, datetime_Lidar, ICON_DF)
    zmlawErr_obs_res = f_resamplingfield(zmlawErr, datetime_Lidar, ICON_DF)
    sigma_w_obs_res  = f_resampling_twoD_Field(sigma_w, datetime_Lidar, \
                                             heightLidar, ICON_DF, ICON_DF_T)
    print('PBL variables resampled: w, sigmaw, Hwind, PBLclass, skewnessW, shear, wdir, mixing layer height')

    # filename for LH estimate based on standard deviation over 30min of vertical velocity, mean w, skewness etc.
    filename = 'wind_dbs-3_'+date+'.nc'
    if (os.path.isfile(pathIn+filename)):
        if verboseFlag == '1':
            print('wind lidar/ceilometer data found: reading and resampling data')
        dataFlagArr[7] = 1
        windLidarScan_data = Dataset(pathIn+filename, mode='r')
        time_windScan      = windLidarScan_data.variables['time'][:].copy()
        speed_windScan     = windLidarScan_data.variables['speed'][:].copy()
        dir_windScan       = windLidarScan_data.variables['dir'][:].copy()
        windVec_windScan   = windLidarScan_data.variables['wind_vec'][:].copy()
        height_windScan    = windLidarScan_data.variables['height'][:].copy()
        # converting time in datetime array
        T_unix             = (time_windScan-2440587.5)*86400
        T_units            = 'seconds since 1970-01-01 00:00:00'
        datetime_windScan  = nc4.num2date(T_unix[:],T_units)
        
        # reading u/v components from the wind vector array
        u_windScan         = windVec_windScan[0,:,:].T
        v_windScan         = windVec_windScan[1,:,:].T
        
        # resampling u, v components on icon time/height grid
        u_obs_res          = f_resampling_twoD_Field(u_windScan, datetime_windScan, \
                                              height_windScan, ICON_DF, ICON_DF_T)
        v_obs_res          = f_resampling_twoD_Field(v_windScan, datetime_windScan, \
                                              height_windScan, ICON_DF, ICON_DF_T)
    # -----------------------------------------------------------------------------------
    # ---- microwave radiometer joyce
    # -----------------------------------------------------------------------------------
    pathIn             = path_mwr_joyce+yy[2:4]+mm+'/'
    file_lwp_iwv_joyce = date[2:8]+'_hps_l2a.nc'
    file_qtprof_joyce  = date[2:8]+'_hps_l2b.nc'    
    if (os.path.isfile(pathIn+file_lwp_iwv_joyce)) and (os.path.isfile(pathIn+file_qtprof_joyce)):
        
    
        print('MWR data LWP/IWV sunhat found: reading and resampling data')
        dataFlagArr[3]    = 1
        LWP_IWV_obs_joyce = Dataset(pathIn+file_lwp_iwv_joyce, mode='r')
        QTprof_obs_joyce  = Dataset(pathIn+file_qtprof_joyce, mode='r')
    
    
        # ----- reading MWR radiometer data variables for la and lb formats
        datetime_lwp_iwv_joyce    = nc4.num2date(LWP_IWV_obs_joyce.variables['time'][:], \
                                         LWP_IWV_obs_joyce.variables['time'].units) 
        IWV_obs_joyce             = LWP_IWV_obs_joyce.variables['atmosphere_water_vapor_content'][:].copy()      # integrated water vapor [kg m^-2]
        LWP_obs_joyce             = LWP_IWV_obs_joyce.variables['atmosphere_liquid_water_content'][:].copy()    # liquid water path [kg m^-2]
        RHsurf_obs_joyce          = LWP_IWV_obs_joyce.variables['relative_humidity'][:].copy()
        Tprof_obs_joyce           = QTprof_obs_joyce.variables['tprof'][:].copy()         # air temperature [K]
        Qprof_obs_joyce           = QTprof_obs_joyce.variables['qprof'][:].copy()  
        datetime_qtprof_joyce     = nc4.num2date(QTprof_obs_joyce.variables['time'][:], \
                                                 QTprof_obs_joyce.variables['time'].units)
        height_qtprof_joyce       = QTprof_obs_joyce.variables['z'][:].copy() 


        # resampling LWP and IWV on ICON time resolution and on height resolution for matrices of interest (T, Q)
        LWP_obs_res    = f_downscaleScalarfield(LWP_obs_joyce, datetime_lwp_iwv_joyce, ICON_DF)
        IWV_obs_res    = f_downscaleScalarfield(IWV_obs_joyce, datetime_lwp_iwv_joyce, ICON_DF)
        tProf_obs_time = f_downscalevectorfield(Tprof_obs_joyce, datetime_qtprof_joyce, height_qtprof_joyce, ICON_DF)        
        qProf_obs_time = f_downscalevectorfield(Qprof_obs_joyce, datetime_qtprof_joyce, height_qtprof_joyce, ICON_DF)    
        #tflip          = np.flip(tProf_obs_time.values,1) # reversing order of heigth for being consisten witn ICON LEM
        #qflip          = np.flip(qProf_obs_time.values,1) # reversing order of heigth for being consisten witn ICON LEM
        qProf_obs_res  = f_resampling_twoD_Field(qProf_obs_time, datetime_ICON, height_qtprof_joyce, ICON_DF, ICON_DF_T)
        tProf_obs_res  = f_resampling_twoD_Field(tProf_obs_time, datetime_ICON, height_qtprof_joyce, ICON_DF, ICON_DF_T)
    
    dictObsWindLidarMwr = {
            'datetime_obs':datetime_lwp_iwv_joyce,
            'verticalWind':w_obs_res.values.transpose(),
            'horizontalWind':Hwind_obs_res.values.transpose(),
            'zonalWind':u_obs_res.values.transpose(),
            'meridionalWind':v_obs_res.values.transpose(), 
            'stdW':sigma_w_obs_res.values.transpose(),
            'varianceW':(sigma_w_obs_res.values.transpose())**2,
            'skewnessW':skew_obs_res.values.transpose(),
            'PBLclassObs':PBLclass_obs_res.values.transpose(),
            'shear':shear_obs_res.values.transpose(),
            'windDirection':wDir_obs_res.values.transpose(),
            'mixingLayerHeight_w_obs': zmlaw_obs_res.values.transpose(),
            'mixingLayerHeightErr_w_obs': zmlawErr_obs_res.values.transpose(),
            'absoluteHumidity':qProf_obs_res.values.transpose(),
            'temperature':tProf_obs_res.values.transpose(),
            'IWV_mwr':IWV_obs_res, 
            'LWP_mwr':LWP_obs_res,
            'LWC_cloudnet':LWC_obs_res,
            'height':height_ICON
            }
    

    
#%%
    
    # ----------------------------------------------------------------------------------
    # ---- calculating cloud properties for iconlem and obs
    # ---------------------------------------------------------------------------------
    # obs
    device            = 'obs'
    # conversion of Ze into linear scale
    #Ze = 10.*np.log10(Ze_lin) formula reminder
    Ze_lin = np.zeros((len(datetime_ICON), len(height_ICON)))
    Ze_lin.fill(np.nan)
    Ze_lin[:,:] = ((10.)**(Z_cloudnet_res.data/10.)).T
    
    LWC_obs           = Ze_lin
#    LWC_obs           = LWC_obs_res.values.transpose()
    cloudDict_obs, clouds_obs, PBLclouds_obs     = f_calculateAllCloudQuantities(cloudnet_res, \
                                                      datetime_ICON, \
                                                      height_ICON, \
                                                      LWP_obs_res, \
                                                      LWC_obs, \
                                                      humanInfo, \
                                                      Hwind_obs_res.values.transpose(), \
                                                      w_obs_res.values.transpose(),\
                                                      int(yy), \
                                                      int(dd), \
                                                      int(mm), \
                                                      QiThreshold, \
                                                      QcThreshold, \
                                                      iconLemData, \
                                                      device, \
                                                      verboseFlag, \
                                                      debuggingFlag, \
                                                      pathDebugFig)
    if device == 'obs':
        print('shape of meanHeight out of the function')
        print(np.shape(cloudDict_obs['meanheightFromCB']))

    print('dimensione time in cloudDict', len(cloudDict_obs['timeSerie']))
    print('dimensione CB in clouddict', len(cloudDict_obs['timeSerie']))

    #CloudInfo, time, height, LWP, LWC, Hwind, Wwind, yy, dd, mm, QiThreshold, QcThreshold, iconLemData, device, verboseFlag):
    # iconlem
    device            = 'iconlem' 
    w_iconlem         = iconLemData.groups['Temp_data'].variables['vertWind'][:].copy()
    LWP_iconlem       = iconLemData.groups['Temp_data'].variables['LWP'][:].copy()   
    LWC_iconlem       = iconLemData.groups['Temp_data'].variables['Qc'][:].copy()  
    LWC_iconlem[LWC_iconlem < QcThreshold] = np.nan
    Hwind_iconlem     = iconLemData.groups['Temp_data'].variables['windSpeed'][:].copy()
    
    cloudDict_iconlem, clouds_iconlem, PBLclouds_iconlem = f_calculateAllCloudQuantities(cloudMask_ICON, \
                                                      datetime_ICON, \
                                                      height_ICON, \
                                                      LWP_iconlem, \
                                                      LWC_iconlem, \
                                                      humanInfo, \
                                                      Hwind_iconlem, \
                                                      w_iconlem,\
                                                      int(yy), \
                                                      int(dd), \
                                                      int(mm), \
                                                      QiThreshold, \
                                                      QcThreshold, \
                                                      iconLemData, \
                                                      device, \
                                                      verboseFlag, \
                                                      debuggingFlag, \
                                                      pathDebugFig)

    

    # -----------------------------------------------------------------------------------
    # ---- sensible and latent heat fluxes from pyranometers at joyce
    # -----------------------------------------------------------------------------------

    # reading latent and sensible heat fluxes for observations and model output
    SHFL_iconlem      = iconLemData.groups['Temp_data'].variables['SHFL'][:].copy()
    LHFL_iconlem      = iconLemData.groups['Temp_data'].variables['LHFL'][:].copy()
    
    # resampling ICON lem variables on the same time resolution of the observations (mean half hourly)
    SHFL_DF           = pd.Series(SHFL_iconlem, index=datetime_ICON)
    SHFL_30min        = SHFL_DF.resample('30min').mean()    
    LHFL_DF           = pd.Series(LHFL_iconlem, index=datetime_ICON)
    LHFL_30min        = LHFL_DF.resample('30min').mean()        
    datetime_30m      = [datetime.datetime(int(yy),int(mm),int(dd),0,0,0) + \
                    datetime.timedelta(minutes=30*x) for x in range(0, 49)]

    # reading surface latent and sensible heat fluxes from the observations
    FluxesData        = Dataset(path_LH_SH_fluxes+'meanSurface_LHSH_Fluxes_ME_RU_SE_'+date+'.nc')
    SHF_obs           = FluxesData.variables['SHF_MEAN_ncdf'][:]
    SHF_Err_obs       = FluxesData.variables['ERR_SHF_MEAN_ncdf'][:]
    LHF_obs           = FluxesData.variables['LHF_MEAN_ncdf'][:]
    LHF_Err_obs       = FluxesData.variables['ERR_LHF_MEAN_ncdf'][:]


    # -----------------------------------------------------------------------------------
    # ---- radiation fluxes (longwave/shortwave)
    # -----------------------------------------------------------------------------------
    TempSurficonLem  = iconLemData.groups['Temp_data'].variables['TempSurf'][:] 
    sigma            = 5.67*10**(-8) # W/m2K4
    UpLWfluxiconLem  = sigma*TempSurficonLem**4 
    LWiconLem        = iconLemData.groups['Temp_data'].variables['LWSurfFlux'][:]
    SWiconLem        = iconLemData.groups['Temp_data'].variables['LWSurfFlux'][:]
    LW_mod           = LWiconLem+UpLWfluxiconLem
    SW_DF_mod        = pd.Series(SWiconLem, index=datetime_ICON)
    SW_mod_30min     = SW_DF_mod.resample('30min').mean() 
    LW_DF_mod        = pd.Series(LW_mod, index=datetime_ICON)
    LW_mod_30min     = LW_DF_mod.resample('30min').mean() 
    datetime_30m     = [datetime.datetime(int(yy),int(mm),int(dd),0,0,0) + \
                        datetime.timedelta(minutes=30*x) for x in range(0, 49)]    

    # correction to fluxes in case there are less elements: 
    # we add nans so to make them comparable to observed time series
    if (len(SHFL_30min) < 49):
        NumberNans = 49 - len(SHFL_30min)
        outSerieSHFL = np.append(np.asarray(SHFL_30min.values), np.repeat(np.nan, float(NumberNans)))
        SHFL_30min = pd.Series(outSerieSHFL, index = datetime_30m[:])
        
    
    if (len(LHFL_30min) < 49):
        NumberNans = 49 - len(LHFL_30min)
        outSerieLHFL = np.append(np.asarray(LHFL_30min.values), np.repeat(np.nan, float(NumberNans)))
        LHFL_30min = pd.Series(outSerieLHFL, index = datetime_30m[:])
    
   
    pathLWflux = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/2013/'
    fileLWObs  = 'sups_joy_pyrg00_l1_rlds_v00_'+date+'000000.nc'
    LWobsData  = Dataset(pathLWflux+fileLWObs, mode='r')
    LWobs      = LWobsData.variables['rlds'][:]
    LWobsErr   = LWobsData.variables['rlds_error'][:]
    timeLWObs  = nc4.num2date(LWobsData.variables['time'][:], LWobsData.variables['time'].units) 
    
    pathSWflux = '/data/data_hatpro/jue/hdcp2/radiation_hdcp2/2013/'
    fileSWObs  = 'sups_joy_pyr00_l1_rsd_v00_'+date+'000000.nc'
    SWobsData  = Dataset(pathSWflux+fileSWObs, mode='r')
    SWobs      = SWobsData.variables['rsd'][:]
    SWobsErr   = SWobsData.variables['rsd_error'][:]
    timeSWObs  = nc4.num2date(SWobsData.variables['time'][:], SWobsData.variables['time'].units) 


    # resampling obvs to icon resolution
    LWobs_res = f_resamplingfield(LWobs, timeLWObs, ICON_DF)
    SWobs_res = f_resamplingfield(SWobs, timeSWObs, ICON_DF)
    LWobsErr_res = f_resamplingfield(LWobsErr, timeLWObs, ICON_DF)
    SWobsErr_res = f_resamplingfield(SWobsErr, timeSWObs, ICON_DF)
    
    
    # calculating half hourly mean values of the longwave and shortwave surface fluxes for obs and model 
    SW_DF_obs        = pd.Series(SWobs_res.values[:,0], index=datetime_ICON)
    SW_obs_30min     = SW_DF_obs.resample('30min').mean() 
    LW_DF_obs        = pd.Series(LWobs_res.values[:,0], index=datetime_ICON)
    LW_obs_30min     = LW_DF_obs.resample('30min').mean() 
    
    LW_DF_obs_Err    = pd.Series(LWobsErr_res.values[:,0], index=datetime_ICON)
    LW_obs_Err_30min = LW_DF_obs_Err.resample('30min').mean() 
    SW_DF_obs_Err    = pd.Series(SWobsErr_res.values[:,0], index=datetime_ICON)
    SW_obs_Err_30min = SW_DF_obs_Err.resample('30min').mean() 
    


    dict_surface_fluxes = {
            'SHF_iconlem':SHFL_30min.values,
            'LHF_iconlem':LHFL_30min.values,
            'datetime_30m':datetime_30m, 
            'SHF_obs':SHF_obs, 
            'LHF_obs':LHF_obs, 
            'SHF_Err_obs':SHF_Err_obs, 
            'LHF_Err_obs':LHF_Err_obs,
            'LW_iconlem':LW_mod_30min.values,
            'SW_iconlem':SW_mod_30min.values,
            'LW_obs':LW_obs_30min.values,
            'SW_obs':SW_obs_30min.values,
            'SW_Err_obs':SW_obs_Err_30min.values,
            'LW_Err_obs':LW_obs_Err_30min.values,
            }
#%%    
    
    # -----------------------------------------------------------------------------------
    # ---- reading and storing additional icon lem variables
    # -----------------------------------------------------------------------------------    
    print('reading dynamics properties for iconlem')
    IWV_iconlem             = iconLemData.groups['Temp_data'].variables['IWV'][:].copy()  # in [K]
    LTS_iconlem             = iconLemData.groups['Temp_data'].variables['LTS'][:].copy()
    PBLheightRN_iconlem     = iconLemData.groups['Temp_data'].variables['PBLHeightArrRN'][:].copy()
    PBLheightTW_iconlem     = iconLemData.groups['Temp_data'].variables['PBLHeightArrRN'][:].copy()
    Tmatrix_iconlem         = iconLemData.groups['Temp_data'].variables['T'][:].copy()  # in [K]
    w_iconlem               = iconLemData.groups['Temp_data'].variables['varianceW'][:].copy()
    u_iconlem               = iconLemData.groups['Temp_data'].variables['zonalWind'][:].copy()
    v_iconlem               = iconLemData.groups['Temp_data'].variables['merWind'][:].copy()
    varianceW_icon_lem      = iconLemData.groups['Temp_data'].variables['varianceW'][:].copy()
    stdW_icon_lem           = iconLemData.groups['Temp_data'].variables['stdWmatrix'][:].copy()
    windSpeed_icon_lem      = iconLemData.groups['Temp_data'].variables['windSpeed'][:].copy()
    windDirection_icon_lem  = iconLemData.groups['Temp_data'].variables['windDirection'][:].copy()
    skenwessW_icon_lem      = iconLemData.groups['Temp_data'].variables['skewnessW'][:].copy()



    dict_iconlem_variables = {
            'IWV_iconlem':IWV_iconlem, 
            'LTS_iconlem':LTS_iconlem,
            'PBLheightRN_iconlem':PBLheightRN_iconlem,
            'PBLheightTW_iconlem':PBLheightTW_iconlem,
            'datetime_iconlem':datetime_ICON,
            'T_iconlem':Tmatrix_iconlem, 
            'u_iconlem':u_iconlem,
            'v_iconlem':v_iconlem,
            'w_iconlem':w_iconlem,
            'varianceW':varianceW_icon_lem, 
            'stdW':stdW_icon_lem,
            'skewnessW':skenwessW_icon_lem,
            'PBLHeightRN':PBLheightRN_iconlem,
            'PBLHeightTW':PBLheightTW_iconlem,
            'windSpeed':windSpeed_icon_lem,#windData_ICON['windSpeed'], 
            'windDirection':windDirection_icon_lem, #windData_ICON['windDirection'], 
            }
    
    # -----------------------------------------------------------------------------------
    # ---- calculating thermodynamics: mixing ratio, RH, virtual temp, virtual pot temp, ccl, lcl 
    # -----------------------------------------------------------------------------------           
    print('calculating thermodynamic properties for COSMO')
    # reading P, Q, T for cosmo
    model             = 'cosmo'
    Qmatrix_cosmo     = Q_cosmo_res.values.transpose()
    Pmatrix_cosmo     = P_cosmo_res.values.transpose()
    Tmatrix_cosmo     = T_cosmo_res.values.transpose()
    #print('P cosmo ', Pmatrix_cosmo[0,0], Pmatrix_cosmo[0,-1])
    Thermodyn_cosmo   = f_calcThermodynamics(Pmatrix_cosmo, Qmatrix_cosmo, Tmatrix_cosmo, \
                                             LTS_iconlem, datetime_ICON, height_ICON, Hsurf, date)
    #(P,Q,T, LTS, time, height, Hsurf
    #print('ccl cosmo')
    #print(Thermodyn_cosmo['cclHeight'])

    #print('calculating thermodynamic properties for observations')
    ## reading P, from iconlem, q, t from microwave radiometer
    print('calculating thermodynamic properties for MWR radiometer obs')
    model             = 'mwr'
    Pmatrix_obs       = iconLemData.groups['Temp_data'].variables['P'][:].copy()
    Qmatrix_obs       = qProf_obs_res.values.transpose()
    Tmatrix_obs       = tProf_obs_res.values.transpose()
    Thermodyn_obs     = f_calcThermodynamics(Pmatrix_obs, Qmatrix_obs, Tmatrix_obs, \
                                            LTS_iconlem, datetime_ICON, height_ICON, Hsurf, date)
    #print('P obs ', Pmatrix_obs[0,0], Pmatrix_obs[0,-1])
    # to debug: lcl array from microwave comes out of 3459 elements instead of 9600 and I don't know why
    #print('ccl obs')
    #print(Thermodyn_obs['cclHeight'])

    print('calculating thermodynamic properties for iconlem')
    # reading P, Q, T for iconlem
    model             = 'iconlem'
    Pmatrix_iconlem   = iconLemData.groups['Temp_data'].variables['P'][:].copy()  # in [Pa]
    Qmatrix_iconlem   = iconLemData.groups['Temp_data'].variables['q'][:].copy()  # in [kg/kg]
    LTS_iconlem             = iconLemData.groups['Temp_data'].variables['LTS'][:].copy()

    #print('P icon ', Pmatrix_iconlem[0,0], Pmatrix_iconlem[0,-1])
    Thermodyn_iconlem = f_calcThermodynamics(Pmatrix_iconlem, Qmatrix_iconlem, Tmatrix_iconlem, \
                                             LTS_iconlem, datetime_ICON, height_ICON, Hsurf, date)
    #print('ccl iconlem')
    #print(Thermodyn_iconlem['cclHeight'])

    # ----------------------------------------------------------------------------------
    # ---- writing data produced in a pickle file and compressing it
    # ---------------------------------------------------------------------------------    
    outputArray   = [radiosondeList, \
                     tower_dict, \
                     dictCosmo, \
                     dictObsWindLidarMwr, \
                     Thermodyn_cosmo, \
                     Thermodyn_iconlem, \
                     Thermodyn_obs, \
                     cloudDict_iconlem, \
                     cloudDict_obs, \
                     dict_iconlem_variables, \
                     dict_surface_fluxes]
    fileOutPickle1 = pathOut+'dictionaries_ModObs_'+date+'.p'

    outfile1       = open(fileOutPickle1,'wb')
    pickle.dump(outputArray, outfile1)
    outfile1.close()

    clouds_obs.to_netcdf(pathOut+'Clouds_Obs_'+date+'.nc')
    clouds_iconlem.to_netcdf(pathOut+'Clouds_iconlem_'+date+'.nc')
    PBLclouds_obs.to_netcdf(pathOut+'PBLClouds_Obs_'+date+'.nc')
    PBLclouds_iconlem.to_netcdf(pathOut+'PBLClouds_iconlem_'+date+'.nc')

    #with open(pathOut+'dataset_PBLcloudPaper_ModObs_'+date+'.p', 'wb') as fp:
    #    pickle.dump(radiosondeList, tower_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
    #import bz2
    #fileOutPickleSmaller = pathOut+'dataset_PBLcloudPaper_ModObs_'+date+'smaller.p'
    #sfile                = bz2.BZ2File(fileOutPickleSmaller, 'w')
    #pickle.dump(outputArray, sfile)
    

    
    # ----------------------------------------------------------------------------------
    # ---- plotting quantities for testing data 
    # ---------------------------------------------------------------------------------      
# =============================================================================
#     if debuggingFlag == 1:
#         
#         # ---------------cloud fraction --------------------------
#         datetime_out_CF  = cloudDict_iconlem['datetimeCloudFraction']
#         height           = cloudDict_iconlem['heightCloudFraction']
#         mean_CF_tot_ICON = cloudDict_iconlem['totalCloudFraction']
#         fig, ax          = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         label_size = 16
#         mpl.rcParams['xtick.labelsize'] = label_size 
#         mpl.rcParams['ytick.labelsize'] = label_size
#         cax1             = ax.pcolormesh(datetime_out_CF, height, mean_CF_tot_ICON.transpose(), vmin=0., vmax=1., cmap='BuPu')
#         ax.set_ylim(400.,4000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("cloud fraction iconlem - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         #cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="cloud fraction ",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'cloudFraction_tot_iconlem_'+date+'.png', format='png')    
# 
#         # ------------------cloud properties -------------------------
#         duration_obs = cloudDict_obs['duration']
#         chordLength_obs = cloudDict_obs['chordLength']
#         massFlux_obs = cloudDict_obs['massFlux']
#         cloudLWP_obs = cloudDict_obs['cloudLWP']
#         duration_iconlem = cloudDict_iconlem['duration']
#         chordLength_iconlem = cloudDict_iconlem['chordLength']
#         massFlux_iconlem = cloudDict_iconlem['massFlux']
#         cloudLWP_iconlem = cloudDict_iconlem['cloudLWP']        
#         nbins      = 10
#         ymax=0.05
#         fig = plt.figure(figsize=(15,7))
#         plt.gcf().subplots_adjust(bottom=0.15)
#         #fig.tight_layout()
#         stringplot = ['a)','b)','c)','d)']
#         stringmidHours = ['cloud duration [s]', 'chord length [m]', 'Mass flux [Kg/ms]', 'mean LWP [g/m^2]']
#         HistQuantities_obs = [duration_obs, chordLength_obs, \
#                               massFlux_obs, np.asarray(cloudLWP_obs)*1000.]
#         HistQuantities_iconlem = [duration_iconlem, chordLength_iconlem, \
#                                   massFlux_iconlem, np.asarray(cloudLWP_iconlem)*1000.]
#         
#         Rangearr = [[0., 1500.], [0., 6000.], [-1500., 1500.], [0., 300.]]
#         ymaxArr = [0.0025, 0.0009, 0.0015, 0.025]
#         xlabelArr = [1400., 5000., 1000., 250.]
#         ylabelArr = [0.007, 0.0022, 0.0037, 0.056]
#         colorArr = ['blue', 'red', 'green', 'purple']
#         for indPlot in range(1,5):
#             ax = plt.subplot(2,2,indPlot)  
#             ax.spines["top"].set_visible(False)  
#             ax.spines["right"].set_visible(False)  
#             ax.get_xaxis().tick_bottom()    
#             ax.get_yaxis().tick_left() 
#             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
#         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
#         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
#             ax.set_ylabel('occurrences')
#             ax.set_xlabel(stringmidHours[indPlot-1])   
#             plt.hist(HistQuantities_obs[indPlot-1], \
#                      bins=nbins, \
#                      normed=True, \
#                      color='black', \
#                      range=Rangearr[indPlot-1], \
#                      cumulative=False, \
#                      alpha=0.5)    
#             plt.hist(HistQuantities_iconlem[indPlot-1], \
#                      bins=nbins, \
#                      normed=True, \
#                      color=colorArr[indPlot-1], \
#                      range=Rangearr[indPlot-1], \
#                      cumulative=False, \
#                      alpha=0.5)  
#         plt.savefig(pathDebugFig+'histograms_cloudProperties_'+date+'.png', format='png')
#             
# 
#         datetime_out_CF  = cloudDict_obs['datetimeCloudFraction']
#         height           = cloudDict_obs['heightCloudFraction']
#         mean_CF_tot_ICON = cloudDict_obs['totalCloudFraction']
#         fig, ax          = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         label_size = 16
#         mpl.rcParams['xtick.labelsize'] = label_size 
#         mpl.rcParams['ytick.labelsize'] = label_size
#         cax1             = ax.pcolormesh(datetime_out_CF, height, mean_CF_tot_ICON.transpose(), vmin=0., vmax=1., cmap='BuPu')
#         ax.set_ylim(400.,4000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("cloud fraction obs - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         #cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="cloud fraction ",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'cloudFraction_tot_obs_'+date+'.png', format='png')                
#     
# #    if debuggingFlag == 1: 
# =============================================================================
        
# =============================================================================
#         # ---------------ccl / lcl / pbl height time series  --------------------------
#         T_ccl_iconlem  = Thermodyn_iconlem['cclTemperature']
#         T_surf_iconlem = Thermodyn_iconlem['surfaceTemperature']
#         lcl_iconlem    = Thermodyn_iconlem['lclHeight']
#         T_ccl_cosmo    = Thermodyn_cosmo['cclTemperature']
#         T_surf_cosmo   = Thermodyn_cosmo['surfaceTemperature']   
#         lcl_cosmo      = Thermodyn_cosmo['lclHeight']
#         #T_ccl_obs      = Thermodyn_obs['cclTemperature']
#         #T_surf_obs     = Thermodyn_obs['surfaceTemperature']   
#         #lcl_obs        = Thermodyn_obs['lclHeight']        
#         
#         T_ccl_radios     = []
#         T_surf_radios    = []
#         datetimeRadios   = []
#         lcl_radios       = []
#         PBLHeight_radios = []
#         for ind in range(len(radiosondeList)):
#             T_ccl_radios.append(radiosondeList[ind]['T_ccl'])
#             T_surf_radios.append(radiosondeList[ind]['surfaceTemperature'])
#             datetimeRadios.append(radiosondeList[ind]['time'])
#             lcl_radios.append(radiosondeList[ind]['z_lcl'])
#             PBLHeight_radios.append(radiosondeList[ind]['PBLheight'])
# 
#         PBLheight_iconlem = iconLemData.groups['Temp_data'].variables['PBLHeightArr'][:].copy()
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         label_size = 16
#         mpl.rcParams['xtick.labelsize'] = label_size 
#         mpl.rcParams['ytick.labelsize'] = label_size
#         plt.title('PBLheight '+str(date), fontsize=18)
#         plt.ylabel('PBLheight [m]', fontsize=16)
#         plt.ylim(0., 12000.)
#         plt.xlabel('time [hh:mm]', fontsize=16)
#         #plt.xlim(datetime.datetime(2013,5,2,6,0,0),datetime.datetime(2013,5,2,20,0,0))
#         plt.plot(datetime_ICON, PBLheight_iconlem, label='ICON-LEM', color='red')
#         #plt.plot(datetime_ICON, T_ccl_obs-T_surf_obs, label='mwr-obs', color='black')
#         plt.plot(np.asarray(datetimeRadios), np.asarray(PBLHeight_radios), 'o', markersize=12, \
#                  label='radiosondes', color='black')
#         plt.axhline(y=0., color='black', linestyle='--')
#         plt.plot()
#         plt.legend(loc='lower left', fontsize=12)
#         plt.tight_layout()
#         plt.savefig(pathDebugFig+'PBLheight_comparison_'+date+'.png', format='png')            
#         
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         label_size = 16
#         mpl.rcParams['xtick.labelsize'] = label_size 
#         mpl.rcParams['ytick.labelsize'] = label_size
#         plt.title(' Tccl - Tsurf  '+str(date), fontsize=18)
#         plt.ylabel('Tccl - Tsurf [K]', fontsize=16)
#         plt.ylim(-5, 25)
#         plt.xlabel('time [hh:mm]', fontsize=16)
#         #plt.xlim(datetime.datetime(2013,5,2,6,0,0),datetime.datetime(2013,5,2,20,0,0))
#         plt.plot(datetime_ICON, T_ccl_iconlem-T_surf_iconlem, label='ICON-LEM', color='red')
#         #plt.plot(datetime_ICON, T_ccl_cosmo-T_surf_cosmo, label='COSMO-EU', color='green')
#         #plt.plot(datetime_ICON, T_ccl_obs-T_surf_obs, label='mwr-obs', color='black')
#         plt.plot(np.asarray(datetimeRadios), np.asarray(T_ccl_radios)-np.asarray(T_surf_radios), 'o', markersize=12, label='radiosondes', color='black')
#         plt.axhline(y=0., color='black', linestyle='--')
#         plt.plot()
#         plt.legend(loc='lower left', fontsize=12)
#         plt.tight_layout()
#         plt.savefig(pathDebugFig+'TCCL_Tsurf_comparison_'+date+'.png', format='png')
#         
# 
# 
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.spines["top"].set_visible(False)  
#         ax.spines["right"].set_visible(False)  
#         ax.get_xaxis().tick_bottom()  
#         ax.get_yaxis().tick_left() 
#         ax.xaxis_date()
#         label_size = 16
#         mpl.rcParams['xtick.labelsize'] = label_size 
#         mpl.rcParams['ytick.labelsize'] = label_size
#         plt.title(' LCL height '+str(date), fontsize=18)
#         plt.ylabel('LCL height [m]', fontsize=16)
#         plt.ylim(0.,12000.)
#         plt.xlabel('time [hh:mm]', fontsize=16)
#         #plt.xlim(datetime.datetime(2013,5,2,6,0,0),datetime.datetime(2013,5,2,20,0,0))
#         plt.plot(datetime_ICON, lcl_iconlem, label='ICON-LEM', color='red')
#         #plt.plot(datetime_ICON, lcl_cosmo, label='COSMO-EU', color='green')
#         #plt.plot(datetime_ICON, lcl_obs, label='mwr-obs', color='black')
#         plt.plot(np.asarray(datetimeRadios), np.asarray(lcl_radios), 'o', markersize=12, label='radiosondes', color='black')
#         plt.axhline(y=0., color='black', linestyle='--')
#         plt.plot()
#         plt.legend(loc='lower left', fontsize=12)
#         plt.tight_layout()
#         plt.savefig(pathDebugFig+'LCLheight_comparison_'+date+'.png', format='png')        
# 
# =============================================================================


# =============================================================================
#     if debuggingFlag == 1: 
# =============================================================================
    
# =============================================================================
#         # ----------------- horizontal wind and absolute humidity obs/iconlem ------------------------------
#         print('plotting time/height plot for Ze')
#         w_plot = np.ma.array(Hwind_obs_res.values.transpose(), mask=np.isnan(Hwind_obs_res.values.transpose()))
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis_date()
#         cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
#         ax.set_ylim(400.,3000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("Horizontal wind observations (wind lidar) - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="horizontal wind [m/s]",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'Hwind_obs_joyce_'+date+'.png', format='png')
#         
#         w_plot = np.ma.array(Qmatrix_obs, mask=np.isnan(Qmatrix_obs))
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis_date()
#         cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
#         ax.set_ylim(400.,3000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("q  mwr - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="q [K]",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'q_obs_joyce_'+date+'.png', format='png')
# 
#     
#         w_plot = np.ma.array(Hwind_iconlem, mask=np.isnan(Hwind_iconlem))
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis_date()
#         cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
#         ax.set_ylim(400.,3000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("Horizontal wind iconlem - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="horizontal wind [m/s]",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'Hwind_iconlem_joyce_'+date+'.png', format='png')
#         
#         w_plot = np.ma.array(Qmatrix_iconlem, mask=np.isnan(Qmatrix_iconlem))
#         fig, ax = plt.subplots(figsize=(12,6))
#         ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
#         ax.xaxis_date()
#         cax1 = ax.pcolormesh(datetime_ICON, height_ICON, w_plot.transpose(), vmin=np.nanmin(w_plot), vmax=np.nanmax(w_plot), cmap='plasma')
#         ax.set_ylim(400.,3000.)                                               # limits of the y-axe
#         ax.set_xlim()                                                        # limits of the x-axes
#         ax.set_title("q iconlem - JOYCE", fontsize=16)
#         ax.set_xlabel("time [hh:mm]", fontsize=16)
#         ax.set_ylabel("height [m]", fontsize=16)
#         #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
#         #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
#         plt.legend(loc='upper left')
#         cbar = fig.colorbar(cax1, orientation='vertical')
#         cbar.ticks=([0,1,2,3])
#         #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
#         cbar.set_label(label="q [K]",size=14)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.aspect=80
#         fig.tight_layout()
#         plt.savefig(pathDebugFig+'q_iconlem_joyce_'+date+'.png', format='png')
# =============================================================================

    print('FINE PROCESSING '+date)
    