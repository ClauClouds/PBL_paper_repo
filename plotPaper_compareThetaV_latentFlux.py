
"""
Created on wed May 27 17:16:28 2020
@author: cacquist
code to compare the area between the thetaV curves at t2 and t1 with the integral of the latent heat flux between t1 and t2.


"""

import numpy as np
import matplotlib
import scipy
import scipy.integrate as integrate
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


def f_calcintegralSummationInf(y,x, xtype):
    """
    date: 28/05/2020
    author: Claudia Acquistapace
    goal: calculate the integral under and arbitrary curve
    method: the principle is to sum all the infinitesimal rectangular shapes generated as yi * (xi-xi-1) over i, where i
    runs through the dimension of y. At each step the area is summed up to the total area.
    The total area is then returned
    input: - y, np.array containing the arbitrary curve to integrate
           - x, np,array of the x coordinate over which to integrate
           - xtype, string indicating the type of the xarray. Options accepted are
                - float: standard processing using float  numnbers
                - datetime: in this case, datetime increments are converted to seconds and then considered as floats.
    output: - intSum, result of the integration
    """
    SumTot = 0.
    Nsteps = len(y)

    for indInt in range(1,Nsteps):
        if xtype == 'float':
            AreaRectangular = y[indInt]* (abs(x[indInt]-x[indInt-1]))
            SumTot = SumTot + AreaRectangular
        if xtype == 'datetime':
            deltax_sec = abs(x[indInt]-x[indInt-1]).total_seconds()
            AreaRectangular = y[indInt]* deltax_sec
            SumTot = SumTot + AreaRectangular
    return(SumTot)


patch = 'patch003'


pathObs = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_' + patch + '/'
pathMod = '/work/cacquist/HDCP2_S2/statistics/iconLemProcessed_' + patch + '/'
pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/' + patch + '/figures_JAMES/debugging/'

fileListObs = sorted(glob.glob(pathObs + '*.p'))
fileListMod = sorted(glob.glob(pathMod + '*icon_lem*.nc'))
date_arr = []
for indFile in range(1):
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


    # reading model data from lem for fluxes at surface and atm indeces
    theta_v_mod  = new_dict[5]['virtualPotentialTemperature']
    height_mod   = new_dict[5]['height']
    LHSF_mod     = -new_dict[10]['LHF_iconlem']
    time_mod     = new_dict[5]['time']
    datetime_30m = new_dict[10]['datetime_30m']


# establishing T1 and T2 over which to integrate
T1 = datetime.datetime(2013,4,14,6,0,0)#datetime_30m[0]
T2 = datetime.datetime(2013,4,14,12,0,0)#datetime_30m[1]

print(type(datetime_30m))
# establish maximum height under which to integrate
Height_2500 = np.min(height_mod[height_mod > 2500])
indHeight_2500 = np.argmin(height_mod[height_mod > 2500])


# find indeces of time array for fluxes that are closest to t1 and t2
timeFlux = np.asarray(datetime_30m)
indT1flux = np.where(timeFlux >= T1)[0]
indT2flux = np.where(timeFlux >= T2)[0]
print(timeFlux[indT1flux[0]])
print(timeFlux[indT2flux[0]])

LatentHeatSerie = np.asarray(LHSF_mod[indT1flux[0]:indT2flux[0]])
TimeSerie = timeFlux[indT1flux[0]:indT2flux[0]]
print(np.shape(LatentHeatSerie))
print(np.shape(TimeSerie))

# calculate integral of latent heat flux between t1 and t2
intLatentHeat = f_calcintegralSummationInf(LatentHeatSerie, TimeSerie, 'datetime')
print('integral of latent heat ',intLatentHeat)



# find index in time_mod closest to t1 and T2
#print(np.where(time_mod >= T2))
#print(np.where(time_mod >= T1))
indT1 = np.where(time_mod >= T1)[0]
indT2 = np.where(time_mod >= T2)[0]

# selecting corresponding thetav profiles
thetaV_T1 = theta_v_mod[indT1[0],indHeight_2500:-1]
thetaV_T2 = theta_v_mod[indT2[0],indHeight_2500:-1]

# calculating integral in time of the latent heat flux


int_T1 = f_calcintegralSummationInf(thetaV_T1, height_mod[indHeight_2500:-1], 'float')
int_T2 = f_calcintegralSummationInf(thetaV_T2, height_mod[indHeight_2500:-1],'float')
print('integrals ', int_T1, int_T2)
IntCompare = int_T2 - int_T1
print('integral difference : ', IntCompare)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
matplotlib.rcParams['savefig.dpi'] = 100
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ax = plt.subplot(1,1,1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.xlim(0., height_mod[-1])
plt.ylim(0., 2500.)
plt.xlabel('theta_v', fontsize=16)
plt.ylabel('height', fontsize=16)
plt.grid(b=True, which='major', color='#666666', linestyle=':')
plt.plot
plt.plot(thetaV_T1, height_mod[indHeight_2500:-1], color='red', label='t = 6:00')
plt.plot(thetaV_T2, height_mod[indHeight_2500:-1], color='black', label='t = 12:00')
plt.fill_betweenx(height_mod[indHeight_2500:-1], thetaV_T1, thetaV_T2)
plt.legend()
plt.tight_layout()
plt.savefig(pathFig+'Thetav_integral.png', format='png')

##print(time_mod[indT1[0]])
#print(time_mod[indT2][0])



