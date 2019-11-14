import numpy as np
import matplotlib
import scipy
import netCDF4 as nc4
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import struct
import glob
import pandas as pd
import datetime as dt
#import xarray as xr
from numpy import convolve


# closest function
#---------------------------------------------------------------------------------
# date :  16.10.2017
# author: Claudia Acquistapace
# goal: return the index of the element of the input array that in closest to the value provided to the function
def f_closest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx  



# function to calculate the convective condensation level height and temperature
#---------------------------------------------------------------------------------
# date :  17.05.2018
# author: Claudia Acquistapace
# goal: function that calculates the convective condensation level (CCL) height and temperature. 
# for the definition of the CCL check this: https://en.wikipedia.org/wiki/Convective_condensation_level
# input: 
# - T field (time, height)
# - RH field (time, height)
# - P field (time, height)
# - height [in m]
# - datetime [in datetime format]
# output:
# - Z_ccl (time) time serie of the height of the CCL 
# - T_CCl (time) time serie of the temperature a parcel should have to reach the height for condensation
# - T_dew point (time, height ) dew point temperature field for every time, height
# method: the function first calculates the dew point temperature field and the dew point at the surface for every time. Then, we derive the saturation mixing ratio at the surface for T=Td. Then, we calculate the mixing ratio field
#--------------------------------------------------------------------------------
def f_CCL(T_ICON, P_ICON, RH_ICON, height_ICON, datetime_ICON, Hsurf):

    # defining constants
    cost_rvl = np.power(5423,-1.) #K
    E0 = 0.611 # Kpa
    T0 = 273. # K 
    Rv = 461 # J K^-1 Kg^-1
    epsilon = 0.622 
    Ad_rate = -9.8 # K/Km

    # ---- substituting RH = 0. to RH = nan to avoid log(0) cases
    RH_ICON [ RH_ICON == 0.] = np.nan
    T_ICON [ T_ICON == 0.] = np.nan

    # ---- calculating due point temperature profile for each time (The dew point is \
    # the temperature to which air must be cooled to become saturated with water vapor. )
    Td = np.power(np.power(T_ICON,-1.)-cost_rvl*np.log(RH_ICON/100.),-1.)

    # ---- calculating mixing ratio at the surface for T dew point 
    Td_surf = Td[:, 149]
    P_surf = P_ICON[:,149]
    RH_surf = RH_ICON[:,149]


    # ---- calculating the saturation mixing ratio for td at the surface (assuming RH =100%)
    M0 = epsilon*E0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))) / (P_surf - E0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))))

    # ---- calculating mixing ratio profile for each P,T, RH using profile RH
    m = (RH_ICON/100.)*epsilon*E0*np.exp((1./Rv)*(T0**(-1.)-T_ICON**(-1.))) / (P_ICON - (RH_ICON/100.)*E0*np.exp((1./Rv)*(T0**(-1.)-T_ICON**(-1.))))

    
    #fig, ax = plt.subplots(figsize=(12,5))

    #plt.plot(m[5000,:],height_ICON)
    #plt.plot(np.repeat(M0[5000],len(height_ICON)), height_ICON)
    #plt.ylim(0,6000)
    
    
    # ---- calculating indeces of height of Z_ccl and Z_CCL height
    ind_CCL = []
    dimHeight=len(height_ICON)
    print(dimHeight)
    for indTime in range(len(datetime_ICON)):
        for indHeight in range(dimHeight-1,1,-1):
            #print(height_ICON[indHeight])
            if (m[indTime,indHeight] < M0[indTime] and m[indTime,indHeight-1] > M0[indTime] ):
                ind_CCL.append(indHeight)
                break
            if indHeight == 1:
                ind_CCL.append(150)
        
    z_ccl = height_ICON[ind_CCL]

    
    # ---- finding z(CCL) using the dry adiabatic lapse rate
    T_ground_CCL = []
    for indTime in range(len(datetime_ICON)):
        T_top = T_ICON[indTime, ind_CCL[indTime]]
        T_ground_CCL.append(T_top - Ad_rate* z_ccl[indTime]*10.**(-3))
    
    dict_out={'z_ccl':z_ccl, 
              'T_ccl':T_ground_CCL, 
              'Td':Td
    }
    return(dict_out)


# variance of vertical velocity calculation
#---------------------------------------------------------------------------------
# date :  17.01.2018
# author: Claudia Acquistapace
# goal: function that calculates the variance of the vertical velocity matrix
# input: 
# - matrix of vertical velocity
# - time array
# - height array
# - time window for the running mean (30 min for comparing to obs)
#--------------------------------------------------------------------------------
def f_calcWvariance(Wwind,time,height,window):

    dimTime = len(time)
    dimHeight = len(height)
    variance = np.zeros((dimTime,dimHeight))
    
    for iheight in range(dimHeight):
        
        # reading array of w values at a given height
        Warray = Wwind[:,iheight]
        
        # calculating running mean variance over the selected array
        # variance[:,iheight] = runningMeanVariance(Warray,window)\
        variance[:,iheight] = np.power(pd.rolling_std(Warray, window),2) 
    
    return variance




def f_calcTheta_MWR(QA_resampled, TA_resampled, PbaroSurf, datetime_ICON, height_TA):

    # date: 11/05/2018
    # author: Claudia Acquistapace
    # goal: calculate P, theta, theta equivalent, mixing ratio using microwave radiometer (absolute humidity and air temperature)
    # data and tower surface pressure obs
    # input: 
    # - absolute humidity QA_resampled on ICON res
    # - air temperature TA_resampled on ICON res
    # - surface pressure from tower measurements resampled on ICON res
    # - datetime array for ICON res
    # - height array of mucrowave obs 
    # output: a dictionary containing:
    #         'height':height_TA,
    #         "P_baro":P_baro,
    #         'theta':theta,
    #         'theta_e':theta_e,
    #         'mixingRatio':mr
    #         'Td':Td,
    #         'z_ccl':z_ccl,
    #         'T_ccl':T_ground_CCL,       
    # --- Defining constants needed for calc
    Cp = 1004.
    Rw = 462.
    Rl = 287.
    g = 9.81
    cost_rlv= 5423.**(-1) # K^-1 (Lv/Rv)
    e0 = 0.611 # Kpa
    T0 = 273.15 # K   

    # --- calculating virtual temperature in K
    T = TA_resampled.values
    q = QA_resampled.values
    Tv = np.zeros((len(datetime_ICON),len(height_TA)))
    for indTime in range(len(datetime_ICON)):
        for indHeight in range(len(height_TA)):
            Tv[indTime,indHeight] = T[indTime,indHeight] + 0.608 * q[indTime,indHeight]


    
    # --- calculating pressure at each level using barometric height formula in Pascal
    Pbaro = np.zeros((len(datetime_ICON),len(height_TA)))

    for indTime in range(len(datetime_ICON)):
        for indHeight in range(len(height_TA)):
            if indHeight == 0:
                Pbaro[indTime, indHeight] = PbaroSurf[indTime]
            else:
                dz = height_TA[indHeight] - height_TA[indHeight-1]
                deltaT = Tv[indTime,indHeight]+Tv[indTime,indHeight-1]
                #print(dz)
                #print(deltaT)
                Pbaro[indTime, indHeight] = Pbaro[indTime, indHeight-1] *np.exp(-g*(dz)/(Rl*(deltaT)/2))
    
    # --- calculating theta and theta_eq
    Theta = np.zeros((len(datetime_ICON),len(height_TA)))
    Theta_e = np.zeros((len(datetime_ICON),len(height_TA)))
    mr = np.zeros((len(datetime_ICON),len(height_TA)))
    Td = np.zeros((len(datetime_ICON),len(height_TA)))

    for indTime in range(len(datetime_ICON)):
        for indHeight in range(len(height_TA)):
            Theta[indTime, indHeight] = T[indTime, indHeight]*(np.power(100000./Pbaro[indTime, indHeight], Rl/Cp)) # potential temperature in K
            e = 462.*T[indTime, indHeight]*q[indTime, indHeight] # partial water vapor pressure in Pa
            mr[indTime, indHeight] = 0.622*e/(Pbaro[indTime, indHeight]-e) # water vapor mixing ratio in kg/kg
            lv = (2500.-2.42*(T[indTime, indHeight]-273.15))*1000. # latent heat of vaporization in J/kg
            Theta_e[indTime, indHeight] = Theta[indTime, indHeight]+(lv*mr[indTime, indHeight]/Cp) \
            *(np.power(100000./Pbaro[indTime, indHeight], Rl/Cp)) # equivalent potential temperature in K
            Td[indTime, indHeight] = (T0**(-1) - cost_rlv * (np.log((e*(10.**(-3)))/e0)))**(-1)

    # ---- calculating the saturation mixing ratio for td at the surface (assuming RH =100%)
    # defining constants
    cost_rvl = np.power(5423,-1.) #K
    Rv = 461 # J K^-1 Kg^-1
    epsilon = 0.622  # g (vapor)/ g (dry air)

    M0 =[]
    for IndTime in range(len(datetime_ICON)):
        M0.append(epsilon*e0*np.exp((1./Rv)*(T0**(-1.)-Td[IndTime,0]**(-1.))) / \
        (PbaroSurf[IndTime]*(10**(-3)) - e0*np.exp((1./Rv)*(T0**(-1.)-Td[IndTime,0]**(-1.)))) )
        
    # --- calculating z_ccl and T_ccl
    # calculating at which height the mixing ratio crosses the surface value
    ind_CCL = []
    for indTime in range(len(datetime_ICON)):
        for indHeight in range(1, len(height_TA)-1):
            #print(height_ICON[indHeight])
            if ((mr[indTime,indHeight] < M0[indTime]) and (mr[indTime,indHeight-1] > M0[indTime])):
                ind_CCL.append(indHeight)
                break

    z_ccl = height_TA[ind_CCL]
    Ad_rate = -9.8 # K/Km
    T_ground_CCL = []
    # ---- finding z(CCL) using the dry adiabatic lapse rate
    for indTime in range(len(datetime_ICON)):
        T_top = T[indTime, ind_CCL[indTime]]
        T_ground_CCL.append(T_top - Ad_rate* z_ccl[indTime]*10.**(-3))        

    dictOut={
        'height':height_TA,
        "Pbaro":Pbaro,
        'theta':Theta,
        'theta_e':Theta_e,
        'mixingRatio':mr,
        'Td':Td,
        'z_ccl':z_ccl,
        'T_ccl':T_ground_CCL,       
    }

    return(dictOut)

def f_readingTowerData(date, PathIn):

    # reading file dat to have rel humidity data
    #filename=PathIn+'meteoFZJ'+date+'.dat'
    #array_smart = np.loadtxt(filename) # numpy é intelligente e riconosce da solo che le linee che iniziano con # sono header

    # read corresponding ncdf file to have all the other variables
    ncfile = 'sups_joy_mett00_l1_any_v00_'+date+'000000.nc'
    ncData = Dataset(PathIn+ncfile, mode='r')
    
    # reading variables from ncdf
    datetime_tower = nc4.num2date(ncData.variables['time'][:],ncData.variables['time'].units) 
    height = ncData.variables['zag'][:].copy()
    T = ncData.variables['ta'][:].copy()
    P = ncData.variables['pa'][:].copy()
    windSpeed = ncData.variables['wspeed'][:].copy()
    wDir = ncData.variables['wdir'][:].copy()
    cols = ['epoch time (Sekunden seit 1.1.1970)','Jahr','Monat [1..12]','Tag [1..31]','Stunde MEZ [0..23]',\
        'Minute [0..59]','Windgeschwindigkeit(2m)','Windgeschwindigkeit(10m)','Windgeschwindigkeit(20m)',\
        'Windgeschwindigkeit(30m)','Windgeschwindigkeit(50m)','Windgeschwindigkeit(80m)','Windgeschwindigkeit(100m)',\
        'Windgeschwindigkeit(120m)','Temperatur(PT100,2m)','Temperatur(PT100,10m)','Temperatur(PT100,20m)',\
        'Temperatur(PT100,30m)','Temperatur(PT100,50m)','Temperatur(PT100,80m)','Temperatur(PT100,100m)',\
        'Temperatur(PT100,120m)','rel.Feuchte(2m)','rel.Feuchte(10m)','rel.Feuchte(20m)','rel.Feuchte(30m)',\
        'rel.Feuchte(50m)','rel.Feuchte(80m)','rel.Feuchte(100m)','rel.Feuchte(120m)','Windrichtung(30m)',\
        'Windrichtung(50m)','Windrichtung(120m)','Einstrahlung(30m)','Ausstrahlung(30m)',\
        'Temperatur_Strahlungsbilanz(30m)','Einstrahlung(120m)','Ausstrahlung(120m)',\
        'Temperatur_Strahlungsbilanz(120m)','Sonnenscheindauer(20m)','Niederschlag(1m)','Luftdruck(2m)',\
        'FFmaxboe(30m)']
    #RHArray = array_smart[0:len(datetime_tower),23:30]
    # reading surface values
    #RHsurf = RHArray[:,0]
    Tsurf = T[:,0]
    
    
    # saving output in a dictionary
    dictOut={
        'time':datetime_tower, 
        'T':T, 
        'P':P, 
        'windSpeed':windSpeed,
        'wDir':wDir,
    #    'RH':RHArray,
        'height':height,
        'Tsurf':Tsurf,
    #    'RHsurf':RHsurf
    }
    return(dictOut)
    
    
def hourDecimal_to_datetime(year, month, day, time):
    a = float('nan')
    timeOut=[]
    for indTime in range(len(time)):
        hours = (int(time[indTime]))
        if hours > 23:
            #print(hours)
            hours = 0
            timeOut.append(dt.datetime(year=year, month=month, day=day, hour=hours, minute=minutes, second=seconds))
        else:
            minutes = int((time[indTime]*60) % 60)
            seconds = int((time[indTime]*3600) % 60)
            timeOut.append(dt.datetime(year=year, month=month, day=day, hour=hours, minute=minutes, second=seconds))

    return(timeOut)



def getNearestIndex(timeRef, timeStamp):
# this function finds the nearest element of timeRef array to the value timeStamp within the given tolerance 
# and returns the index of the element found. If non is found within the given tolerance, it returns nan.
    try:
        index = timeRef.index.get_loc(timeStamp, method='nearest')
    
    except:

        index = np.nan
    
    return index

def getIndexList(dataTable, reference):
#   this function reads the less resolved time array (dataTable) and the time array to be used as reference (reference) 
#   and the tolerance. then for every value in the reference array, it finds the index of the nearest element of 
#   dataTable for a fixed tolerance. It provides as output the list of indeces of dataTable corresponding 
#   to the closest elements of the reference array. 
    #print(len(reference))
    indexList = []
    for value in reference:
        #print(value)
        index = getNearestIndex(dataTable, value)
        indexList.append(index)

    return indexList

def getIndexListsuka(dataTable, reference):
    indexList = []
    for value in reference:
        #print(value)
        index = dataTable.index.get_loc(value, method='nearest')
        indexList.append(index)

    return indexList    
    

def getResampledDataPd(emptyDataFrame, LessResolvedDataFrame, indexList):
#  it reads the dataframe to be filled with the resampled data (emptyDataFrame), then the originally less resolved
#  data (dataDataFrame) and the list of indeces of the less resolved time array upsampled 
#  to the highly resolved resolutions. Then, with a loop on the indeces of the indexList, 
#  It assigns to the emptydataframe the values of the less resolved dataframe called by the   
#  corresponding index of the indexlist. The output is the filled emptydataFrame
 

    for i, index in enumerate(indexList):

        try:
            emptyDataFrame.iloc[i]=LessResolvedDataFrame.iloc[index]
        
        except:
            pass
        
    return emptyDataFrame







 


# function to define cb and ct of the first cloud layer in a column appearing when reading from the top
#---------------------------------------------------------------------------------
# date :  31.01.2019
# author: Claudia Acquistapace
# goal: return the index of the element of the input array that in closest to the value provided to the function
def f_calcCloudBaseTop(cloudMask, dimTime, dimHeight, height):
    
    
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
        numberCB.append(len(np.where(column == -1.)[0][:]))
        numberCT.append(len(np.where(column == 1.)[0][:]))  

    NCB=max(numberCB)   
    NCT=max(numberCT)

    # generating cloud base and cloud top arrays 
    CBarray = np.zeros((dimTime,NCB))
    CBarray.fill(np.nan)
    CTarray = np.zeros((dimTime,NCT))
    CTarray.fill(np.nan)

    # storing cloud base and cloud top arrays
    for iTime in range(dimTime):
        column = gradBinary[iTime,:]
        indCB=np.where(column == 1.)[0][:]
        NfoundCB=len(indCB)
        indCT=np.where(column == -1.)[0][:] 
        NfoundCT=len(indCT)
        CBarray[iTime,0:NfoundCB]=height[indCB]
        CTarray[iTime,0:NfoundCT]=height[indCT]
    
    
    
    #fig, ax = plt.subplots(figsize=(10,4))
    #cax = ax.pcolormesh(np.arange(dimTime), height, BinaryMatrix.transpose(), vmin=0, vmax=1, cmap='cool')
    #ax.set_ylim(0,5000.)                                               # limits of the y-axes
    #ax.set_xlim(2000,2100)                                                 # limits of the x-axes
    #plt.plot(np.arange(dimTime), CBarray, color='black', linestyle='-', label='cloud base')
    #plt.plot(np.arange(dimTime), CTarray, color='red',  linestyle=':',label='cloud top')
    #plt.legend(loc='upper left')
    #ax.set_title("cloud mask", fontsize=14)
    #ax.set_xlabel("time ", fontsize=12)
    #ax.set_ylabel("height [m]", fontsize=12)
    #cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical')
    #cbar.ticks=([0,1,2,3])
    #cbar.ax.set_yticklabels(['no cloud','liquid','ice', 'mixed phase'])
    #cbar.set_label(label="cloud type",size=12)
    #cbar.ax.tick_params(labelsize=12)
    #cbar.aspect=80

    
    return (CBarray,CTarray)





# moving variance calculation function
#---------------------------------------------------------------------------------
# date :  17.01.2018
# author: Claudia Acquistapace
# goal: function that calculates the moving average of an array of values over a given window given as imput
#--------------------------------------------------------------------------------
def runningMeanVariance(x, N):
    return np.power(pd.rolling_std(x, N),2)


# variance of vertical velocity calculation
#---------------------------------------------------------------------------------
# date :  17.01.2018
# author: Claudia Acquistapace
# goal: function that calculates the variance of the vertical velocity matrix
#--------------------------------------------------------------------------------
def f_calcWvariance(Wwind,time,height,window):
    dimTime = len(time)
    dimHeight = len(height)
    variance = np.zeros((dimTime,dimHeight))
    
    for iheight in range(dimHeight):
        
        # reading array of w values at a given height
        Warray = Wwind[:,iheight]
        
        # calculating running mean over the selected array
        variance[:,iheight] = runningMeanVariance(Warray,window)
   
    
    return variance
        
    
# skewness of vertical velocity calculation
#---------------------------------------------------------------------------------
# date :  17.01.2018
# author: Claudia Acquistapace
# goal: function that calculates the variance of the vertical velocity matrix
#--------------------------------------------------------------------------------
# calculating running mean skewness matrix
def f_runningMeanSkewnessW(time, timeWindowSk, runningWindow, height, vertWind):
    
    dimTime=len(time)
    dimHeight=len(height)
    SKmatrix = np.zeros((dimTime, dimHeight))
    
    for iTime in range(0, dimTime-1, timeWindowSk):

        # for all indeces of time 
        if (iTime > runningWindow/2) & (iTime < dimTime-1-runningWindow/2):
            # generating indeces to read corresponding elements in the matrix
            timeIndeces = np.arange(iTime-runningWindow/2, iTime+runningWindow/2, dtype=np.int16)
            #print(iTime)
            for iHeight in range(dimHeight):
            
                meanW = np.mean(vertWind[timeIndeces, iHeight])            # mean of the wind array
                variance = np.var(vertWind[timeIndeces, iHeight])         # variance of the wind array
                wprime = np.subtract(vertWind[timeIndeces, iHeight],np.tile(meanW,len(timeIndeces)))
                num = np.mean(np.power(wprime,3))
                den = variance**(3./2.)
                SKmatrix[iTime:iTime+timeWindowSk-1, iHeight] = num/den
                
    return SKmatrix




# cloud mask for ice/liquid/mixed phase clouds 
#---------------------------------------------------------------------------------
# date :  19.01.2018
# author: Claudia Acquistapace
# goal: function that calculates cloud mask for the day
# input: 
# -Qc: cloud liquid content
# -Qi: ice liquid content
# -QcThreshold: Threshold value for detection of Qc
# -QiThreshold: Threshold value for detection of Qi
# output:
# -cloudmask(dimTime,dimheight) containing: 1=liquid clouds, 2=ice clouds, 3=mixed phase clouds
#--------------------------------------------------------------------------------
def f_cloudmask(time,height,Qc,Qi,QiThreshold,QcThreshold):
    
    dimTime=len(time)
    dimHeight=len(height)
    cloudMask = np.zeros((dimTime, dimHeight))
    # = 0 : no clouds
    # = 1 : liquid clouds
    # = 2 : ice clouds
    # = 3 : mixed phase clouds
    print(np.shape(Qi))

    for iTime in range(dimTime):
        for iHeight in range(dimHeight):

            # Ice and not liquid
            if (Qi[iTime, iHeight] > QiThreshold) and (Qc[iTime, iHeight] < QcThreshold): 
                cloudMask[iTime, iHeight] = 2.          # ice clouds 
            # liquid and not ice
            if (Qi[iTime, iHeight] < QiThreshold) and (Qc[iTime, iHeight] > QcThreshold): 
                cloudMask[iTime, iHeight] = 1.          # liquid clouds 
            # liquid and ice 
            if (Qi[iTime, iHeight] >= QiThreshold) and (Qc[iTime, iHeight] >= QcThreshold): 
                cloudMask[iTime, iHeight] = 3.          # mixed phase clouds 
    return cloudMask


def f_cloudmask2(time,height,Qc,Qi,Qs, Qg, QcThreshold,QiThreshold, QsThreshold, QgThreshold):    
    dimTime=len(time)
    dimHeight=len(height)
    cloudMask = np.zeros((dimTime, dimHeight))
    # = 0 : no clouds
    # = 1 : liquid clouds
    # = 2 : ice clouds
    # = 3 : mixed phase clouds
    print(np.shape(Qi))

    for iTime in range(dimTime):
        for iHeight in range(dimHeight):

            # Ice and not liquid
            if ((Qi[iTime, iHeight] > QiThreshold) or (Qs[iTime, iHeight] > QsThreshold) or (Qg[iTime, iHeight] > QgThreshold))\
            and (Qc[iTime, iHeight] < QcThreshold):
                cloudMask[iTime, iHeight] = 2.          # ice clouds 
            # liquid and not ice
            if (Qi[iTime, iHeight] < QiThreshold) and (Qc[iTime, iHeight] > QcThreshold) \
            and (Qs[iTime, iHeight] < QsThreshold) and (Qg[iTime, iHeight] < QgThreshold):
                cloudMask[iTime, iHeight] = 1.          # liquid clouds 
            # liquid and ice 
            if ((Qi[iTime, iHeight] > QiThreshold) or (Qs[iTime, iHeight] > QsThreshold) or (Qg[iTime, iHeight] > QgThreshold))\
            and (Qc[iTime, iHeight] >= QcThreshold):
                cloudMask[iTime, iHeight] = 3.          # mixed phase clouds 
    return cloudMask



def f_cloudmask_INSCAPE(time,height,Qc,QcThreshold):    
    # cloud mask selecting only liquid clouds 
    dimTime=len(time)
    dimHeight=len(height)
    cloudMask = np.zeros((dimTime, dimHeight))
    # = 0 : no clouds
    # = 1 : liquid clouds
    # = 2 : ice clouds
    # = 3 : mixed phase clouds


    for iTime in range(dimTime):
        for iHeight in range(dimHeight):

            # Ice and not liquid
            if (Qc[iTime, iHeight] < QcThreshold):
                cloudMask[iTime, iHeight] = 1.          # ice clouds 
    return cloudMask



# PBL cloud classification 
#---------------------------------------------------------------------------------
# date :  25.01.2018
# author: Claudia Acquistapace
# goal: function that derives PBL classification
# input: 
# -time: time array
# -height: height array
# -gradWindThr: threshold to detect wind shear regions
# -SigmaWThres: Threshold to detect turbulence
# -ylim: array of heights up to which classification is calculated
# -cloudMask: cloud mask indicating cloud presence
# -varianceMatrix: matrix containing variance of vertical velocity to detect turbulence
# -SKmatrix: matrix containing skewness of vertical velocity 
# -StabilityArr: array of time dimension indicating stability of the PBL close to the surface
# -connection2Surface: array of time dimension indicating that turbulence is connected to the surface or not
# -gradWindSpeed: matrix (dimHeight,dimTime) containing the intensity of the shear of the horizontal wind
# -cloudBaseHeightArr: array of time dimension containing cloud base heights
# output:
# -PBLclass(dimTime,dimheight) containing: 1=in cloud, 2=non turb, 3=cloud driven, 4=convective, 5=intermittent, 6=wind shear
#--------------------------------------------------------------------------------
def f_PBLClass(time,height,gradWindThr,SigmaWThres,ylim, cloudMask, varianceWmatrix, SKmatrix, stabilityArr, connection2Surface, gradWindspeed, cloudBaseHeightArr):
    
    dimTime=len(time)
    dimHeight=len(height)
    PBLclass=np.zeros((dimTime, dimHeight))   # defining output matrix
    shear=gradWindspeed.transpose()           # transposing gradWindspeed to conform to other matrices

    # defining flag matrices and filling them with nan values. Each flag corresponds to a check to be performed for the classification. Checks are for:  cloud/turbulence/cloud driven/unstable at surface/wind shear/surface driven
    flagCloud=np.zeros((dimTime, dimHeight))
    flagCloud.fill(np.nan)
    flagTurb=np.zeros((dimTime, dimHeight))
    flagTurb.fill(np.nan)
    flagcloudDriven=np.zeros((dimTime, dimHeight))
    flagcloudDriven.fill(np.nan)
    flagInstability=np.zeros((dimTime, dimHeight))
    flagInstability.fill(np.nan)
    flagConnection=np.zeros((dimTime, dimHeight))
    flagConnection.fill(np.nan)
    flagSurfaceDriven=np.zeros((dimTime, dimHeight))
    flagSurfaceDriven.fill(np.nan)
    flagWindShear=np.zeros((dimTime, dimHeight))
    flagWindShear.fill(np.nan)

    # loop on time and height to assign flags to each pixel: 
    for iTime in range(dimTime):
        for iHeight in range(dimHeight):

            # quitting loop on heights if height is greater than ylim.
            if height[iHeight] > ylim[iTime]:
                iHeight = dimHeight-1

             
            #------------------------------------------------------------
            # check if cloud
            #------------------------------------------------------------
            if (cloudMask[iTime, iHeight] != 0): 
                flagCloud[iTime, iHeight] = 1
            else:
                flagCloud[iTime, iHeight] = 0

            #------------------------------------------------------------
            # check if not in cloud and not turbulent
            #------------------------------------------------------------
            if (varianceWmatrix[iTime,iHeight] >= SigmaWThres): 
                flagTurb[iTime, iHeight] = 1
            else:
                flagTurb[iTime, iHeight] = 0

            #------------------------------------------------------------
            # check if cloud driven ( conditions to pose: 1) below cloud base , 2) sigma > sthr, 3) connected to the cloud
            #------------------------------------------------------------
            if np.count_nonzero(cloudMask[iTime,:]) == 0:
                flagcloudDriven[iTime, iHeight] = 0
            else: 
                indCB=f_closest(height, cloudBaseHeightArr[iTime])   # finding index of cloud base
                if iHeight <= indCB:
                    flagcloudDriven[iTime, iHeight] = 0
                else: 
                    flagConnect = []
                    # check if pixels between cloud base and height[iheight] fullfill conditions for being cloud driven
                    for ind in range(indCB, iHeight):
                        if (varianceWmatrix[iTime,ind] > SigmaWThres) & (SKmatrix[iTime,ind] < 0.):
                            flagConnect.append(1)

                    if (varianceWmatrix[iTime,iHeight] > SigmaWThres) & (SKmatrix[iTime,iHeight] < 0.) & (len(flagConnect) == iHeight - indCB):
                        flagcloudDriven[iTime, iHeight] = 1
                    else:
                        flagcloudDriven[iTime, iHeight] = 0


            #------------------------------------------------------------
            # check if unstable to the surface
            #------------------------------------------------------------        
            if (stabilityArr[iTime] == 0):
                flagInstability[iTime, iHeight] = 1
            else:
                flagInstability[iTime, iHeight] = 0
            #------------------------------------------------------------
            # check if wind shear is present
            #------------------------------------------------------------ 
            if (shear[iTime,iHeight] < gradWindThr):
                flagWindShear[iTime,iHeight] = 0   
            else:
                flagWindShear[iTime,iHeight] = 1 
            #------------------------------------------------------------
            # check if turbulence is surface driven 
            #------------------------------------------------------------ 
            if (connection2Surface[iTime] == 0):
                flagSurfaceDriven[iTime,iHeight] = 0
            else:

                # find min height where variance is bigger than the threshold
                indArray = np.where(varianceWmatrix[iTime,:] > SigmaWThres)[0]
                if len(indArray) > 0:
                    indHmin = np.max(indArray)
                    if indHmin > 200:
                        flagSurfaceDriven[iTime,iHeight] = 0
                    else:
                        counter = indHmin
                        while varianceWmatrix[iTime,counter] > SigmaWThres:
                            counter = counter - 1
                            flagSurfaceDriven[iTime,iHeight] = 1
                else:
                    flagSurfaceDriven[iTime,iHeight] = 0 



    # defining classification by posing conditions as indicated in Manninen et al, 2018, JGR
    for iTime in range(dimTime):
        for iHeight in range(dimHeight):

            # quitting loop on heights if height is greater than ylim.
            if height[iHeight] > ylim[iTime]:
                iHeight = dimHeight-1


            # defining in cloud bins
            if (flagCloud[iTime, iHeight] == 1):
                PBLclass[iTime, iHeight] = 1 

            # defining non turbulent bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 0):
                PBLclass[iTime, iHeight] = 2

            # defining cloud driven bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 1):
                PBLclass[iTime, iHeight] = 3            

            # defining convective bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 0) & (flagInstability[iTime, iHeight] == 1) & (flagSurfaceDriven[iTime, iHeight] == 1):
                PBLclass[iTime, iHeight] = 4

            # defining intermittent bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 0) & (flagInstability[iTime, iHeight] == 1) & (flagSurfaceDriven[iTime, iHeight] == 0):
                PBLclass[iTime, iHeight] = 5        

            # defining wind shear bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 0) & (flagInstability[iTime, iHeight] == 0) & (flagWindShear[iTime, iHeight] == 1):
                PBLclass[iTime, iHeight] = 6
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 0) & (flagInstability[iTime, iHeight] == 0) & (flagSurfaceDriven[iTime, iHeight] == 0) & (flagWindShear[iTime, iHeight] == 1):
                PBLclass[iTime, iHeight] = 6

            # defining intermittent bins
            if (flagCloud[iTime, iHeight] == 0) & (flagTurb[iTime, iHeight] == 1) & (flagcloudDriven[iTime, iHeight] == 0) & (flagInstability[iTime, iHeight] == 0) & (flagWindShear[iTime, iHeight] == 0):
                PBLclass[iTime, iHeight] = 5        

    return (PBLclass, flagCloud, flagTurb, flagcloudDriven, flagInstability, flagWindShear, flagSurfaceDriven, flagConnection)

# list of ordered variables : var_long_name =
#  "Pressure", 0
#  "Temperature", 1
#  "Exner pressure", 2
#  "Density", 3
#  "virtual potential temperature", 4
#  "zonal wind", 5
#  "meridional wind", 6
#  "orthogonal vertical wind", 7
#  "specific humidity", 8
#  "specific cloud water content", 9
#  "specific cloud ice content", 10
#  "rain_mixing_ratio", 11
#  "snow_mixing_ratio", 12
#  "relative humidity", 13
#  "graupel_mixing_ratio", 14
#  "graupel_mixing_ratio", 15
#  "number concentration ice", 16
#  "number concentration snow", 17
#  "number concentration rain droplet", 18 
#  "number concentration graupel", 19
#  "number concentration hail", 20
#  "number concentration cloud water", 21
#  "number concentration activated ice nuclei", 22
#  "total specific humidity (diagnostic)", 23
#  "total specific cloud water content (diagnostic)", 24
#  "total specific cloud ice content (diagnostic)", 25
#  "cloud cover", 26
#  "turbulent diffusion coefficients for momentum", 27
#  "turbulent diffusion coefficients for heat", 28
#  "Pressure on the half levels", 29
#  "soil temperature", 30
#  "total water content (ice + liquid water)", 31
#  "ice content" ; 32


# function to calculate moist adiabatic lapse rate from temperature and mixing ratio
# date: 11 July 2018, Vancouver
# author: claudia Acquistapace
# input: height dimension, time dimension, mixing ratio, temperature and pressure matrix in time, height, lcl height array.
# output: EIS height array
# notes: EIS is calculated here as indicated in eq 3 of Wood and Bretherton, 2006: On the Relationship between Stratiform Low Cloud Cover and Lower-Tropospheric Stability (Journal of Climate)
def f_calcEIS(dimHeight, dimTime, height, mr, T, P, lcl, LTS):
    
    g = 9.8 # gravitational constant [ms^-2]
    Cp = 1005.7 # specific heat at constant pressure of air [J K^-1 Kg^-1]
    Lv = 2256 # latent heat of vaporization of water [kJ/kg]
    R = 8.314472 # gas constant for dry air [J/ molK]
    epsilon = 0.622 # ratio of the gas constants of dry air and water vapor

    # calculating moist adiabatic lapse rate
    gamma_moist =  np.zeros((dimTime,dimHeight))

    for indTime in range(dimTime):
        for indHeight in range(dimHeight):
            #print( (R*T[indTime,indHeight] ))
            #print(T[indTime,indHeight])
        
            gamma_moist[indTime,indHeight] = g*( 1.+ (Lv*mr[indTime,indHeight] )/(R*T[indTime,indHeight] ) \
                                            / Cp + ((Lv**2) * epsilon * mr[indTime,indHeight] )/(R*T[indTime,indHeight]**2)  )

            
    # finding height corresponding to 700 HPa
    z_700 = []
    gamma_moist_700 = []
    gamma_lcl = []
    EIS = []

    for indTime in range(dimTime):
        ind_height_700 = f_closest(P[indTime,:], 70000.)
        ind_lcl = f_closest(height, lcl[indTime])
        z_700.append(((R*T[indTime,0])/g)*np.log(P[indTime,0]/70000.))
        gamma_moist_700.append(gamma_moist[indTime,ind_height_700])
        gamma_lcl.append(gamma_moist[indTime,ind_lcl])
    
    for indTime in range(dimTime):
        EIS.append(LTS[indTime] - (gamma_moist_700[indTime]*z_700[indTime]*0.001 + gamma_lcl[indTime]*lcl[indTime]*0.001))


    return EIS



# Version 1.0 released by David Romps on September 12, 2017.
# 
# When using this code, please cite:
# 
# @article{16lcl,
#   Title   = {Exact expression for the lifting condensation level},
#   Author  = {David M. Romps},
#   Journal = {Journal of the Atmospheric Sciences},
#   Year    = {2017},
#   Volume  = {in press},
# }
#
# This lcl function returns the height of the lifting condensation level
# (LCL) in meters.  The inputs are:
# - p in Pascals
# - T in Kelvins
# - Exactly one of rh, rhl, and rhs (dimensionless, from 0 to 1):
#    * The value of rh is interpreted to be the relative humidity with
#      respect to liquid water if T >= 273.15 K and with respect to ice if
#      T < 273.15 K. 
#    * The value of rhl is interpreted to be the relative humidity with
#      respect to liquid water
#    * The value of rhs is interpreted to be the relative humidity with
#      respect to ice
# - ldl is an optional logical flag.  If true, the lifting deposition
#   level (LDL) is returned instead of the LCL. 
# - min_lcl_ldl is an optional logical flag.  If true, the minimum of the
#   LCL and LDL is returned.
def lcl(p,T,rh=None,rhl=None,rhs=None,return_ldl=False,return_min_lcl_ldl=False):

    import math
    import scipy.special
    import numpy as numpy 
    
    # Parameters
    Ttrip = 273.16     # K
    ptrip = 611.65     # Pa
    E0v   = 2.3740e6   # J/kg
    E0s   = 0.3337e6   # J/kg
    ggr   = 9.81       # m/s^2
    rgasa = 287.04     # J/kg/K 
    rgasv = 461        # J/kg/K 
    cva   = 719        # J/kg/K
    cvv   = 1418       # J/kg/K 
    cvl   = 4119       # J/kg/K 
    cvs   = 1861       # J/kg/K 
    cpa   = cva + rgasa
    cpv   = cvv + rgasv

    # The saturation vapor pressure over liquid water
    def pvstarl(T):
        return ptrip * (T/Ttrip)**((cpv-cvl)/rgasv) * \
         math.exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) )
   
    # The saturation vapor pressure over solid ice
    def pvstars(T):
        return ptrip * (T/Ttrip)**((cpv-cvs)/rgasv) * \
         math.exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) )

    # Calculate pv from rh, rhl, or rhs
    rh_counter = 0
    if rh  is not None:
        rh_counter = rh_counter + 1
    if rhl is not None:
        rh_counter = rh_counter + 1
    if rhs is not None:
        rh_counter = rh_counter + 1
    if rh_counter != 1:
        print(rh_counter)
        exit('Error in lcl: Exactly one of rh, rhl, and rhs must be specified')
    if rh is not None:
        # The variable rh is assumed to be 
        # with respect to liquid if T > Ttrip and 
        # with respect to solid if T < Ttrip
        if T > Ttrip:
            pv = rh * pvstarl(T)
        else:
            pv = rh * pvstars(T)
        rhl = pv / pvstarl(T)
        rhs = pv / pvstars(T)
    elif rhl is not None:
        pv = rhl * pvstarl(T)
        rhs = pv / pvstars(T)
        if T > Ttrip:
            rh = rhl
        else:
            rh = rhs
    elif rhs is not None:
        pv = rhs * pvstars(T)
        rhl = pv / pvstarl(T)
        if T > Ttrip:
            rh = rhl
        else:
            rh = rhs
    if pv > p:
        return np.nan

    # Calculate lcl_liquid and lcl_solid
    qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
    
   
    rgasm = (1-qv)*rgasa + qv*rgasv
    cpm = (1-qv)*cpa + qv*cpv
    if rh == 0:
        return cpm*T/ggr
    aL = -(cpv-cvl)/rgasv + cpm/rgasm
    bL = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
    cL = pv/pvstarl(T)*math.exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
    aS = -(cpv-cvs)/rgasv + cpm/rgasm
    bS = -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
    cS = pv/pvstars(T)*math.exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
    lcl = cpm*T/ggr*( 1 - \
       bL/(aL*scipy.special.lambertw(bL/aL*cL**(1/aL),-1).real) )
    ldl = cpm*T/ggr*( 1 - \
      bS/(aS*scipy.special.lambertw(bS/aS*cS**(1/aS),-1).real) )

    # Return either lcl or ldl
    if return_ldl and return_min_lcl_ldl:
        exit('return_ldl and return_min_lcl_ldl cannot both be true')
    elif return_ldl:
        return ldl
    elif return_min_lcl_ldl:
        return min(lcl,ldl)
    else:
        return lcl
    
    
def f_histo(sample, bins):
     pylab.hist(sample, bins=50, range=(-5.5), histtype = 'step', color = 'black')
     pylab.show()
