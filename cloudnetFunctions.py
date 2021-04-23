#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 15:49:57 2019

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


def CheckIceLiquidCloudnet(array):
    """

    ;author: Claudia Acquistapace
    date of modification: 17/04/2020
    goal : read cloudnet target classification and convert it to a string for calculating cloud fraction amounts.
        MODIFICATION: deal with nans
    :return:
    """
    
    dimArr     = len(array)           # reading dimension for array in input
    stringFlag = []                   # generating output string list
    
    for ind in range(dimArr):
        
        # reading the element of the array
        if np.isnan(array[ind]):
            stringFlag.append('none')
        else:

            elem = bin(int(array[ind]))[2:]        # reading the element in binary base
            #least significative on the right end of the number and cutting first two characters '0b'
            #print array[ind]
            #print bin(array[ind])[2:]

            length=len(elem)           # counting number of digits that represent the number

            # adapting lenght of digits to the maximum lenght possible that can be found in Cloudnet (6 bits) if the length of the string is smaller
            if length < 6:
                Nbins2add = 6 - length               # number of zeros to add to the left
                elem      = '0' * Nbins2add + elem
                # print 'resized elem'
                # print elem

            # flags for bits of cloudnet that are on
            flagBin0 = int((elem)[-1])      # bit 0: small liquid droplets on
            flagBin1 = int((elem)[-2])      # bit 1: falling hydrometeors
            flagBin2 = int((elem)[-3])      # bit 2: wet bulb < 0, if bit 1 on, then phase
            flagBin3 = int((elem)[-4])      # bit 3: melting ice particles

            #print flagBin0, flagBin1, flagBin2, flagBin3
            # condition for only liquid clouds
            if ((flagBin0 == 1) and (flagBin2 == 0) and (flagBin3 == 0)):
                stringFlag.append('liquid') # cloud droplets and drizzle, cloud droplets only
                #print 'sono qui'
            if ((flagBin1 == 1) and (flagBin2 == 1)):
                stringFlag.append('ice')
            if (flagBin3 == 1):
                stringFlag.append('ice')
            if ((flagBin0 == 0) and (flagBin1 == 0) and (flagBin2 == 0) and (flagBin3 == 0)):
                stringFlag.append('none')

    return(stringFlag)



def f_calculateCloudMaskCloudnet(time, height, cloudnet):
    #---------------------------------------------------------------------------------------------------------
    # date : 13 April 2018 
    # author: Claudia Acquistapace
    # abstract : routine to calculate cloud mask matrix from cloudnet target classification
    # input :
    #   - time
    #   - height
    #   - cloudnet categorization in bits
    # output:
    #   - cloud mask
    #---------------------------------------------------------------------------------------------------------      
    dimTime = len(time)
    dimHeight = len(height)
    cloudMask = np.zeros(shape=(dimTime,dimHeight))		# matrix for cloud mask
    
    # loop on hour array for calculating cloud fraction
    for itime in range(dimTime-1):
        for iheight in range(dimHeight):
            
            
            # reading the element of the array
            elem=bin(cloudnet[itime, iheight])[2:]        # reading the element in binary base 
            #least significative on the right end of the number and cutting first two characters '0b' 
            #print array[ind]
            #print bin(array[ind])[2:]
        
            length=len(elem)           # counting number of digits that represent the number
        
            # adapting lenght of digits to the maximum lenght possible that can be found in Cloudnet (6 bits) if the length of the string is smaller
            if length < 6:
                Nbins2add= 6 - length               # number of zeros to add to the left 
                elem = '0' * Nbins2add + elem
                # print 'resized elem'
                # print elem
            
            # flags for bits of cloudnet that are on
            flagBin0 = int((elem)[-1])      # bit 0: small liquid droplets on
            flagBin1 = int((elem)[-2])      # bit 1: falling hydrometeors
            flagBin2 = int((elem)[-3])      # bit 2: wet bulb < 0, if bit 1 on, then phase
            flagBin3 = int((elem)[-4])      # bit 3: melting ice particles
            
                    # condition for only liquid clouds
            if ((flagBin0 == 1) and (flagBin2 == 0) and (flagBin3 == 0)):
                cloudMask[itime,iheight] = 1  # liquid
                #stringFlag.append('liquid') # cloud droplets and drizzle, cloud droplets only
                #print 'sono qui'
            if ((flagBin1 == 1) and (flagBin2 == 1)): 
                cloudMask[itime,iheight] = 2  # ice                
                #stringFlag.append('ice')
            if (flagBin3 == 1):
                cloudMask[itime,iheight] = 2  # ice                
                #stringFlag.append('ice')
            if ((flagBin0 == 0) and (flagBin1 == 0) and (flagBin2 == 0) and (flagBin3 == 0)):
                cloudMask[itime,iheight] = 0  # no cloud               
                #stringFlag.append('none')

    return (cloudMask)





def f_calculateCloudFractionCloudnet(cloudnet, yy, mm, dd, time, height):
    """
    date : 11 April 2018 (modified version from the previous code (see above))
    author: Claudia Acquistapace
    abstract : routine to calculate mean cloud fraction every hour for and mean profile for every six hours of the day
               for total cloud fraction, ice cloud fraction and liquid cloud fraction. This is done for a specific site
               and for the whole day. The routine also provides an output dictionary containing the cloud fraction
               profiles and their standard deviations, for the hourly and six hourly mean and for each phase.
    input :
      - date to process
      - site of the measurements
      - path to the data
    output:
      - dictionary containing data
      - plots of the daily cloud fractions in the
    """

    from cloudnetFunctions import CheckIceLiquidCloudnet
    # building time array for cloud fraction calculation
    deltaT       = datetime.timedelta(minutes=30)
    indInt       = 0
    datetime_out = []
    for itime in range(0,48):
        if indInt == 0:
            HourInf = datetime.datetime(int(yy), int(mm), int(dd), 0, 0, 0) 
        else:
            HourInf = HourInf + deltaT
        HourSup     = HourInf + deltaT
        datetime_out.append(HourInf)
        indInt      = indInt + 1

    # generating output matrices
    CFTCloudnet = np.zeros(shape=(48,len(height)))		# matrix for total cloud fraction (ice+liquid)
    CFICloudnet = np.zeros(shape=(48,len(height)))		# matrix for ice cloud fraction 
    CFLCloudnet = np.zeros(shape=(48,len(height)))		# matrix for liquid cloud fraction 

    # building a dataframe for the cloudnet target classification

    CloudnetDF = pd.DataFrame(cloudnet.transpose(),index=time, columns=height)

    # loop on hour array for calculating cloud fraction
    for itime in range(len(datetime_out)-1):
        for iheight in range(len(height)):
            
            # selecting lines corresponding to measurements within the hour
            CloudnetDFArr = CloudnetDF.loc[(CloudnetDF.index < datetime_out[itime+1]) * \
                                           (CloudnetDF.index >= datetime_out[itime]), height[iheight]]
        
            # applying the function to check for ice/liquid cloud presence:
            # returns an array of flags for the hour, where you find liquid/ice
            # clouds at that heigth level
            stringArr     = CheckIceLiquidCloudnet(CloudnetDFArr.values)
            
            # calculating cloud fraction in probabilistic way Nsel/Ntot
            Ntot          = len(CloudnetDFArr.values)                     # total number of bins in the column
            Nliquid       = stringArr.count('liquid')     # number of liquid clouds
            Nice          = stringArr.count('ice')           # number of ice clouds
            
            # calculating cloud fractions for the time interval selected at the height j
            if Ntot == 0: 
                CFTCloudnet[itime,iheight] = 0.
                CFICloudnet[itime,iheight] = 0.
                CFLCloudnet[itime,iheight] = 0.
            else:
                CFTCloudnet[itime,iheight] = float(Nice+Nliquid)/float(Ntot)
                # total cloud fraction (liquid+ic/mixed phase clouds)
                CFICloudnet[itime,iheight] = float(Nice)/float(Ntot)          # ice cloud fraction
                CFLCloudnet[itime,iheight] = float(Nliquid)/float(Ntot)            # liquid cloud fraction
                
                
            
    
    # defining dictionary containing data to have as output
    dict_CF = {}

    # filling dictionaries with data 
    dict_CF = {
            'TotalCloudFraction':CFTCloudnet,
            'LiquidCloudFraction':CFLCloudnet,
            'IceCloudFraction':CFICloudnet,
            'height':height,
            'time':datetime_out,
                }
    
    
    return (dict_CF)      
    
