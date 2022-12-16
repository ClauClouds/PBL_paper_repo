#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:53:19 2019

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
import math

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle


def f_processRadiosondesDay(fileList, yy, mm, dd):
    # goal: process radiosoundings and calculating quantities of interest for the analysis
    # author: claudia Acquistapace
    # date; 24 July 2019 (heat wave in Cologne)
    # contact: cacquist@meteo.uni-koeln.de
    # output: list of dictionaries. Each dictionary is for a specific hour of the day
    # and contains 'time', 'P', T', 'Td', 'Td_surf','RH','z_lcl','z_ccl','T_ccl'
    # names of the variables in the radiosondes profiles
    cols                   = ['Zeit [min:sec]','P [hPa]','T [C]','U [%]','Wind speed [m/s]','Wdir [inf]','Lange [inf]'\
                              ,'Breite [inf]','Hˆhe [m]','Geo Pot [m]', 'dew [C]', 'Tv [C]','Rs [m/s]', 'D [kg/m3]' \
                              ,'Azimut []','Elevation []','Entfernung [m]']
    Nfiles                 = len(fileList)
    RadiosondesData       = []
    for iFile in range(Nfiles):
        hourFile           = int(fileList[iFile][52:54])
        DatetimeRadiosonde = datetime.datetime(int(yy), int(mm), int(dd), hourFile, 0, 0)
        print(hourFile)

        DF                 = pd.read_csv(fileList[iFile], sep='\t', skipinitialspace=True, \
                                         encoding='latin_1', names=cols, header=0, dtype={'Hˆhe [m]':str})
        DF[DF.columns] = DF[DF.columns].apply(pd.to_numeric, errors='coerce')
        #print(iFile, fileList[iFile], DF.dtypes)
        # ---- reading variables in the file
        P                  = DF.values[:,1] # in [hPa]
        T                  = DF.values[:,2]
        U                  = DF.values[:,3]
        Hwind              = DF.values[:,4]
        Wdir               = DF.values[:,5]
        Td                 = DF.values[:,10]
        Td                 = Td + 273.15
        T                  = T + 273.15
        height             = DF.values[:,8]

# =============================================================================
#         for ind in range(len(height)):
#             #print(ind, height[ind])
#             height[ind] = height[ind].strip()
#             if height[ind] == '-----':
#                 height[ind] =np.nan
#                 height[ind] == float(height[ind])
#             else:
#                 height[ind] = float(height[ind].strip())
#
# =============================================================================
        #print(hourFile)
        Ndata              = len(DF.count(axis='columns'))
        #print(DF.values[:,4])
        #type(DF.values[:,4])

        # ---- Calculating RH for radiosondes from Td (Td = 1/To - Rv/L log(e/eo)^-1) from stull
        # RH = e/es ; es is from clausius clapeyron while e is obtained from the formula for td
        e0                 = 0.611 # Hpa
        T0                 = 273.15 # K
        cost_lrv           = 5423. # K (Lv/Rv)
        e                  = []
        es                 = []
        RH                 = []
        for indData in range(Ndata):
            e.append(e0 * np.exp(cost_lrv*(T0**(-1)-Td[indData]**(-1))))
            es.append(e0 * np.exp(cost_lrv*(T0**(-1)-T[indData]**(-1))))
            RH.append(100*(e0 * np.exp(cost_lrv*(T0**(-1)-Td[indData]**(-1))))/ \
                      (e0 * np.exp(cost_lrv*(T0**(-1)-T[indData]**(-1)))))

        # ---- finding height where RH = 100%
        ind_sat            = []
        for indData in range(1,Ndata,+1):
            if ((RH[indData-1] < 100.) and (RH[indData] > 100.)):
                ind_sat = indData
                break

        #z_sat              = height[ind_sat]
        # ---- calculating CCL from radiosoundings
        # ---- calculating the saturation mixing ratio for td at the surface (assuming RH =100%)
        # defining constants
        #cost_rvl           = np.power(5423,-1.) #K
        Rv                 = 461 # J K^-1 Kg^-1
        epsilon            = 0.622
        Td_surf            = Td[0]
        P_surf             = float(P[0])
        T_surf             = float(T[0])

        M0                 = epsilon*e0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))) / \
        (P_surf - e0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))))

        # ---- calculating mixing ratio profile for each P,T, RH using profile RH
        m                  = []
        arg_exp            = np.array((1./Rv)*(T0**(-1.)-T**(-1.)))
        for indData in range(Ndata):
            m.append((float(RH[indData])/100.)*epsilon*e0*np.exp(arg_exp[indData])\
                 / (float(P[indData]) - (float(RH[indData])/100.)*e0*np.exp(arg_exp[indData])))


        for indData in range(1,Ndata):
            if ((m[indData-1] < M0) and (m[indData] > M0)):
                ind_CCL = indData
                break
        z_ccl              = height[ind_CCL]

        Ad_rate            = -9.8 # K/Km
        T_ground_CCL       = []
            # ---- finding z(CCL) using the dry adiabatic lapse rate
        T_top              = T[ind_CCL]

        #print((z_ccl))
        #print((T_top))
        T_ground_CCL       = float(T_top) - Ad_rate* float(z_ccl)*10.**(-3)

        # ---- calculating LCL height
        #------ calculating LCL heights from tower measurements resampled
        # important: provide pressure in Pascals, T in K, RH in 70.3
        #---------------------------------------------------------------------------------
        from myFunctions import lcl
        z_lcl              = lcl(np.array(P_surf*100),np.array(T_surf),np.array(RH[0])/100.)



        # ------------------------------------------------------------------
        # calculate LTS index for lower tropospheric stability (Wood and Bretherton, 2006)
        # ------------------------------------------------------------------
        from myFunctions import f_closest
        Pthr = 700 * 100. # Pressure level of 700 Hpa used as a reference
        # calculating height of the surface
        indP700 = f_closest(P*100.,Pthr)
        Theta = []
        Cp = 1004.
        Rl = 287.
        for ind in range(len(height)):
            Theta.append(T[ind]*((100000./(P[ind]*100.))**(Rl/Cp))) # potential temperature in K

        LTS = Theta[indP700] - Theta[0]




        #------------------------------------------------------------------
        # calculate EIS index for lower tropospheric stability (Wood and Bretherton, 2006) for observations
        # ------------------------------------------------------------------
        g = 9.8 # gravitational constant [ms^-2]
        Cp = 1005.7 # specific heat at constant pressure of air [J K^-1 Kg^-1]
        Lv = 2256 # latent heat of vaporization of water [kJ/kg]
        R = 8.314472 # gas constant for dry air [J/ molK]
        epsilon = 0.622 # ratio of the gas constants of dry air and water vapor
        gamma_d = g / Cp # dry adiabatic lapse rate


        # ---- calculating mixing ratio
        mr =[]
        for indHeight in range(len(height)):
            mr.append((0.622*e[indHeight]*100.)/(P[indHeight]*100.-e[indHeight]*100.)) # water vapor mixing ratio in kg/kg


        # ---- calculating saturation mixing ratio
        ws = []#mpcalc.saturation_mixing_ratio(P, T)
        for indHeight in range(len(height)):
            ws.append(epsilon* (es[indHeight]/(P[indHeight]- es[indHeight]))) # saturation water vapor mixing ratio kg/Kg


        gamma_moist = []
        gamma_moist_atmos = []
        # calculating moist adiabatic lapse rate
        for indHeight in range(len(height)):
            gamma_moist.append(g*((1.+(Lv*mr[indHeight])/(R*T[indHeight]))/(Cp \
                                  + (mr[indHeight]*epsilon*Lv**2)/(R*T[indHeight]**2))))

        Ws_array = np.asarray(ws)
        T_array = np.asarray(T)
        gamma_moist_atmos  = atmos.equations.Gammam_from_rvs_T(Ws_array.astype(float), T_array.astype(float)) # in [k/m]
        indP700 = f_closest(P*100.,Pthr)
        gamma_moist_700 = gamma_moist[indP700]
        z_700 = height[indP700]
        #print('here')
        #print(float(z_lcl))
        #for ind in range(len(height)):
        #    print((height[ind]))
        #print(f_closest(height, float(z_lcl)))

        ind_lcl = f_closest(height, float(z_lcl))
        gamma_lcl = gamma_moist[ind_lcl]

        LFC = np.nan
        T_LFC = np.nan
# =============================================================================
#         # calculating Level of free convection from LCL
#         if (z_lcl > 0):
#
#             for indloop in range(len(height)-ind_lcl):
#                 indHeight = indloop + ind_lcl
#                 print(gamma_lcl[indHeight]*height[indHeight])
#                 print(T[indHeight])
#                 if gamma_lcl[indHeight]*height[indHeight] < T[indHeight]:
#                     print('LCL')
#                     print(z_lcl)
#                     print('height found')
#                     print(height[indHeight])
#                     print(indHeight)
#                     print(ind_lcl)
#                     LFC = height[indHeight]
#                     T_LFC = T[indHeight]
#
# =============================================================================
        # finding height corresponding to 700 HPa
        EIS = LTS - gamma_moist_700*z_700 + gamma_lcl*float(z_lcl)
        EIS_atmos = LTS - gamma_moist_atmos[indP700]*z_700 + gamma_moist_atmos[ind_lcl]*float(z_lcl)
        #print('EIS obtained from the Wood and Bretherton formula:')
        #print(EIS)


        # calculating profiles of virtual potential temperature
        Theta_v = []
        Rd = 287.058  # gas constant for dry air [Kg-1 K-1 J]
        for indHeight in range(len(height)):
            k = Rd*(1-0.23*mr[indHeight])/Cp
            Theta_v.append( (1 + 0.61 * mr[indHeight]) * T[indHeight] * (1000./P[indHeight])**k)

        # calculating EIS with the methodology of maximum deltaTheta_v
        Delta_ThetaV = [x - Theta_v[i - 1] for i, x in enumerate(Theta_v)][1:]
        # cutting profiles at 4000mt height
        indCut = f_closest(height, 3500.)
        Delta_ThetaV= Delta_ThetaV[0:indCut]
        #print('il massimo shift risulta :')
        #print(np.max(Delta_ThetaV))
        #print(np.argmax(Delta_ThetaV))
        EIS_height = height[np.argmax(Delta_ThetaV)]
        #print('e si trova ad altezza:')
        #print(EIS_height)

        EIS2 = Theta_v[np.argmax(Delta_ThetaV)]- Theta_v[0]
        #print('EIS obtained with the maximum deltaTheta virtual difference:')
        #print(EIS2)
        #print('da qui')

        # ------------------------------------------------------------------------------
        # calculating PBL height
        # ------------------------------------------------------------------------------
        dimHeight = len(height)
        g=9.8                                                # gravity constant
        Rithreshold=0.25                                     # Threshold values for Ri
        Rithreshold2=0.2
        zs=height[0]                                         # height of the surface reference
        RiMatrix=np.zeros((dimHeight))                       # Richardson number matrix
        PBLheightArr=[]
        RiCol=np.zeros((dimHeight))
        # calculating richardson number matrix
        thetaS=Theta_v[0]

        for iHeight in range(dimHeight):
            if isinstance((DF.values[iHeight,4]), float):
                den = DF.values[iHeight,4]**2
            if isinstance((DF.values[iHeight,4]), str):
                if (DF.values[iHeight,4] == '-----    '):
                    den == 0.
                else:
                    den = (float(DF.values[iHeight,4].strip()))**2
            if den == 0.:
                RiMatrix[iHeight] = 0.
            else:
                RiMatrix[iHeight] = (1/den) * (g/thetaS) * (Theta_v[iHeight]-thetaS)*(height[iHeight]-zs)


        # find index in height where Ri > Rithreshold
        RiCol=RiMatrix[:]
        #print(RiCol)
        #print(np.where(RiCol > Rithreshold2)[0][:])
        #print(len(np.where(RiCol > Rithreshold)[0][:]))
        if len(np.where(RiCol > Rithreshold)[0][:]) != 0:
            PBLheight = (height[np.where(RiCol > Rithreshold)[0][0]] - height[0])
        else:
            PBLheight = 0.
        #print('pbl height for the radiosonde:')
        #print(PBLheight)

        # ---- saving variables in dictionary
        dict_day           = {
                'time':DatetimeRadiosonde,
                'P':P,
                'T':T,
                'Td': Td,
                'Td_surf': Td[0],
                'RH':RH,
                'z_lcl':z_lcl,
                'z_ccl':z_ccl,
                'T_ccl':T_ground_CCL,
                'PBLheight':PBLheight,
                'EISWood':EIS,
                'EIS2':EIS2,
                'LTS':LTS,
                'theta_v':Theta_v,
                'surfaceTemperature':T_surf,
                'height':height,
                }
        RadiosondesData.append(dict_day)
    return(RadiosondesData)



# -----------------------------------------------------------------------------------
# ---- radiosoundings
# -----------------------------------------------------------------------------------
date              = '20130518'
yy                = date[0:4]
mm                = date[4:6]
dd                = date[6:8]
path_radiosondes = '/home/juelich/rs_hope/KIT/'
pathIn   = path_radiosondes+yy+mm+dd+'/'
fileList = glob.glob(pathIn+'*.txt')
print(fileList)
if os.listdir(pathIn):
    # folder full, process radiosondes that are in
    radiosondeList = f_processRadiosondesDay(fileList, yy, mm, dd)
