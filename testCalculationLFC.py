#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:25:38 2019

@author: cacquist
"""
import pandas as pd
import numpy as np
from myFunctions import f_closest
import matplotlib.pyplot as plt

#%%

def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure
          
    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)


def SatVap(dwpt,phase="liquid"):
    """ water vapor pressure, 
    
    Inputs: dew poitn temperature
    output: e(Pa) water vapor pressure
    
    SOURCE: bolton, Montlhy Weather review, 1980, p 1047, eq 10
    """

    return(611.2*np.exp(17.67*dwpt/(243.5+dwpt)))


def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk*(1.0+0.6*mixr)

def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """
   
    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)


def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """
    degCtoK = 273.15
    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=MixRatio(es,pres)

    if e is None:
        # assume saturated
        e=es

    w=MixRatio(e,pres)

    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma


#%%
Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218	          # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.80665          # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153e-3         # Mean molar mass of water vapor(kg/mol)
m_a=28.9644e-3        # Mean molar mass of air(kg/mol)
Rstar_a=8.31432       # Universal gas constant for air (N m /(mol K))
gammaDry = -9.8       # C/km or K/km dry adiabatic lapse rate
# -----------------------------------------------------------------------------------
# ---- radiosoundings
# -----------------------------------------------------------------------------------   
date              = '20130502'
yy                = date[0:4]
mm                = date[4:6]
dd                = date[6:8]
path_radiosondes = '/home/juelich/rs_hope/KIT/'
pathIn   = path_radiosondes+yy+mm+dd+'/'
filename = 'KIT_HOPE_2013050217.txt'
file = pathIn+filename


cols                   = ['Zeit [min:sec]','P [hPa]','T [C]','U [%]','Wind speed [m/s]','Wdir [inf]','Lange [inf]'\
                              ,'Breite [inf]','Hˆhe [m]','Geo Pot [m]', 'dew [C]', 'Tv [C]','Rs [m/s]', 'D [kg/m3]' \
                              ,'Azimut []','Elevation []','Entfernung [m]']
DF                 = pd.read_csv(file, sep='\t', skipinitialspace=True, \
                                         encoding='latin_1', names=cols, header=0, dtype={'Hˆhe [m]':str})
DF[DF.columns] = DF[DF.columns].apply(pd.to_numeric, errors='coerce')
Ndata              = len(DF.count(axis='columns'))


#print(iFile, fileList[iFile], DF.dtypes)
# ---- reading variables in the file
P                  = DF.values[:,1]  # in [hPa]
T                  = DF.values[:,2]
U                  = DF.values[:,3]
Hwind              = DF.values[:,4]
Wdir               = DF.values[:,5]
Td                 = DF.values[:,10]
Td                 = Td #+ 273.15
T                  = T #+ 273.15
testCalheight      = DF.values[:,8]
Rv                 = 461 # J K^-1 Kg^-1
epsilon            = 0.622 
Td_surf            = Td[0]
P_surf             = float(P[0])
T_surf             = float(T[0])
height             = DF.values[:,8]


#%%

# ---- Calculating RH for radiosondes from Td (Td = 1/To - Rv/L log(e/eo)^-1) from stull
# RH = e/es ; es is from clausius clapeyron while e is obtained from the formula for td 
g = 9.8 # gravitational constant [ms^-2]
Cp = 1005.7 # specific heat at constant pressure of air [J K^-1 Kg^-1]
Lv = 2256 # latent heat of vaporization of water [kJ/kg]
R = 8.314472 # gas constant for dry air [J/ molK]
epsilon = 0.622 # ratio of the gas constants of dry air and water vapor

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
        
from myFunctions import lcl
z_lcl              = lcl(np.array(P_surf*100),np.array(T_surf),np.array(RH[0])/100.)
T_lcl = (gammaDry* 0.001 * z_lcl) + T[0] # obtained following the dry adiabat from surface to lcl
P0 = P[0] # in Hpa  # 101.325 #Kpa
Hp = 7.29 * 1000  #m scale heigth for pressure 
#P_lcl = P0 * 100. * np.exp (-z_lcl/Hp) # lcl pressure in hPa
ind_lcl = f_closest(height, z_lcl)
P_lcl = P[ind_lcl] * 100.
pres = P[ind_lcl:-1]
temp = T_lcl
Tamb = T[ind_lcl:-1]
height_abovelcl = height[ind_lcl:-1]
t_out=np.zeros(pres.shape)
t_out[0]=T_lcl
z_LFCArr = np.zeros(pres.shape)
z_LFCArr.fill(np.nan)
# loop on pressure levels to calculate decrease of T due to moist adiabat and intersection with T profile
for ii in range(pres.shape[0]-1):
    # interval of pressures 
    delp=pres[ii]-pres[ii+1]
    
    # moist saturated adiabatic decrease of T
    temp=temp-100*delp*GammaW(temp+degCtoK,(pres[ii]-delp/2)*100)
    
    # storing the decrease of T
    t_out[ii+1]=temp
    
    if temp > Tamb[ii]:
        z_LFCArr[ii]= height_abovelcl[ii]
z_LFC = np.min(z_LFCArr)
print(z_LFC)

#%%

plt.plot(T, height, label='T')
plt.plot(Td, height, label='Td')
plt.plot(t_out, height_abovelcl, label = 'sat ad ' )
plt.hlines(z_LFC, -40., 20., label='z_LFC')
plt.hlines(z_lcl, -20., 20., label='z_lcl')
plt.plot((gammaDry* 0.001 * height) + T[0], height, label='dry ad')
plt.legend()
plt.xlim(-20.,20.)
plt.ylim(0., 8000.)
