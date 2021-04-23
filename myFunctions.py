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
from numpy import convolve
import datetime
import atmos
import matplotlib.dates as mdates



#"""
#Created on Wed Nov 13 10:41:35 2019
#functions to define,for each day of the dataset, which is the interval of time 
#to consider to select typical boundary layer clouds and exclude other cloud types
#which can mess up the statistics. the filter is based on human selection based 
#on visual inspection of cloudnet target classification and model cloud mask.
#As general rule: we exclude all clouds with cloud base above 2000 mt and vertical 
#extension larger than 1000 m associated with a mixed / ice phase.

#'20130501']#'20130502']#, '20130501','20130424', '20130425', '20130427', '20130429'
#@author: cacquist
#"""
import xarray as xr


def f_selectingPBLcloudWindow(date):
    """ function to select time intervals for processing cloud data in each day.
    Also, we pick which algorithm to use with each type of data.
    in particular:
        - minmax correspond to select as cloud top the max of the cloud tops found. 
        Clouds above 5000 mt are filtered out (done with the function 
        f_calcCloudBaseTopPBLclouds)
        - version2: selects boundary layer clouds with cloud base below 2500 mt and cloud tops below 
        CB+600mt.
    """
    if date == '20130414':
        timeStart   = datetime.datetime(2013,4,14,6,0,0)
        timeEnd     = datetime.datetime(2013,4,14,23,59,59)
        PBLheight   = 2500.
    if date == '20130420':
        timeStart   = datetime.datetime(2013,4,20,6,0,0)
        timeEnd     = datetime.datetime(2013,4,20,23,59,59)
        PBLheight   = 2000.
    if date == '20130424':
        timeStart   = datetime.datetime(2013,4,24,6,0,0)
        timeEnd     = datetime.datetime(2013,4,24,23,59,59)
        PBLheight   = 2000.
    if date == '20130425':
        timeStart   = datetime.datetime(2013,4,25,6,0,0)
        timeEnd     = datetime.datetime(2013,4,25,23,59,59)
        PBLheight   = 5000.
    if date == '20130426':
        timeStart   = datetime.datetime(2013,4,26,6,0,0)
        timeEnd     = datetime.datetime(2013,4,26,23,59,59)
        PBLheight   = 5000.
    if date == '20130427':
        timeStart   = datetime.datetime(2013,4,27,6,0,0)
        timeEnd     = datetime.datetime(2013,4,27,23,59,59)
        PBLheight   = 3000.
    if date == '20130428':
        timeStart   = datetime.datetime(2013,4,28,6,0,0)
        timeEnd     = datetime.datetime(2013,4,28,23,59,59)
        PBLheight   = 3500.
    if date == '20130429':
        timeStart   = datetime.datetime(2013,4,29,6,0,0)
        timeEnd     = datetime.datetime(2013,4,29,23,59,59)
        PBLheight   = 3000.
    if date == '20130430':
        timeStart   = datetime.datetime(2013,4,30,6,0,0)
        timeEnd     = datetime.datetime(2013,4,30,23,59,59)
        PBLheight   = 3000.
    if date == '20130501':
        timeStart   = datetime.datetime(2013,5,1,6,0,0)
        timeEnd     = datetime.datetime(2013,5,1,23,59,59)
        PBLheight   = 2500.
    if date == '20130502':
        timeStart   = datetime.datetime(2013,5,2,6,0,0)
        timeEnd     = datetime.datetime(2013,5,2,23,59,59)
        PBLheight   = 4000.
    if date == '20130503':
        timeStart   = datetime.datetime(2013,5,3,6,0,0)
        timeEnd     = datetime.datetime(2013,5,3,23,59,59)
        PBLheight   = 3000.
    if date == '20130504':
        timeStart   = datetime.datetime(2013,5,4,6,0,0)
        timeEnd     = datetime.datetime(2013,5,4,23,59,59)
        PBLheight   = 2500.
    if date == '20130505':
        timeStart   = datetime.datetime(2013,5,5,6,0,0)
        timeEnd     = datetime.datetime(2013,5,5,23,59,59)
        PBLheight   = 2500.
    if date == '20130506':
        timeStart   = datetime.datetime(2013,5,6,6,0,0)
        timeEnd     = datetime.datetime(2013,5,6,23,59,59)
        PBLheight   = 3000.
    if date == '20130509':
        timeStart   = datetime.datetime(2013,5,9,6,0,0)
        timeEnd     = datetime.datetime(2013,5,9,23,59,59)
        PBLheight   = 3000.
    if date == '20130510':
        timeStart   = datetime.datetime(2013,5,10,6,0,0)
        timeEnd     = datetime.datetime(2013,5,10,23,59,59)
        PBLheight   = 3000.
    if date == '20130518':
        timeStart   = datetime.datetime(2013,5,18,6,0,0)
        timeEnd     = datetime.datetime(2013,5,18,23,59,59)
        PBLheight   = 2500.
    if date == '20130524':
        timeStart   = datetime.datetime(2013,5,24,6,0,0)
        timeEnd     = datetime.datetime(2013,5,24,23,59,59)
        PBLheight   = 4500.
    if date == '20130525':
        timeStart   = datetime.datetime(2013,5,25,6,0,0)
        timeEnd     = datetime.datetime(2013,5,25,23,59,59)
        PBLheight   = 3000.
    if date == '20130527':
        timeStart   = datetime.datetime(2013,5,27,6,0,0)
        timeEnd     = datetime.datetime(2013,5,27,23,59,59)
        PBLheight   = 3000.
    if date == '20130528':
        timeStart   = datetime.datetime(2013,5,28,6,0,0)
        timeEnd     = datetime.datetime(2013,5,28,23,59,59)
        PBLheight   = 4000.

    dictOut = {'timeStart':timeStart, 'timeEnd':timeEnd, 'heightPBL':PBLheight}
#,'20130504', '20130505','20130506','20130509','20130510'
# '20130414','20130420', '20130426','20130428', '20130430','20130524','20130525','20130527', '20130528'
    return(dictOut)
    
#-------------------------------------------------------------------------------------


def f_calculateCloudBaseTopThickness(cloudMask, time, height, humanInfo):
    """
    date : wednesday 13 may 2020
    author: Claudia Acquistapace
    goal: build a function to identify all cloud base and cloud top of clouds in the vertical profile at the same time.
    Human observations for the day distinguish manually PBL from non-PBL clouds. An additional dataset of
     PBL clouds is delivered based on this information.
    Concept of the code:
    step 1: given the cloud mask, find all cloud base and cloud tops.
    step 2: build cloud database with all clouds saved as xarray dataset
    step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
    input: cloudmask,
            time,
            height,
            humanInfo (dictionary including timeStart, timeEnd, PBLheight from human obs on the day)
    output: AllCloudDataset (xarray Dataset including cloud base, cloud top, cloud thickness, level number)
            PBLcloudDataset (xarray Dataset for PBL clouds with cloud base, cloud top, cloud thickness, level number)
    """
    dimTime   = len(time)
    dimHeight = len(height)
    heightPBL = humanInfo['heightPBL']
    timeStart = humanInfo['timeStart']
    timeEnd   = humanInfo['timeEnd']

    # STEP 1: identifying all cloud bases and tops
    # ---------------------------------------------------
    # converting cloud mask to 1 / 0 matrices
    BinaryMatrix = np.zeros((dimTime, dimHeight))
    for itime in range(dimTime):
        for iH in range(dimHeight):
            if cloudMask[itime, iH] != 0.:
                BinaryMatrix[itime, iH] = 1

    # calculating gradient of binary cloud mask
    gradBinary = np.diff(BinaryMatrix, axis=1)

    # counting max number of cloud base/cloud top found
    numberCB = []
    numberCT = []
    for itime in range(dimTime):
        column = gradBinary[itime, :]
        numberCB.append(len(np.where(column == -1.)[0][:]))
        numberCT.append(len(np.where(column == 1.)[0][:]))

    NCB = max(numberCB)
    NCT = max(numberCT)

    # generating cloud base and cloud top arrays
    CBarray = np.zeros((dimTime, NCB))
    CBarray.fill(np.nan)
    CTarray = np.zeros((dimTime, NCT))
    CTarray.fill(np.nan)
    NlayersArray = np.zeros((dimTime))
    NlayersArray.fill(np.nan)

    # if no cloud bases or no cloud tops are found, then CB and CT are assigned to nan
    if (NCB == 0) or (NCT == 0):
        CBarray[iTime, :] = np.nan
        CTarray[iTime, :] = np.nan
    else:
        # if some cloud base / cloud tops are found, all the found values are stored
        # storing cloud base and cloud top arrays
        for iTime in range(dimTime):
            column = gradBinary[iTime, :]
            indCB = np.where(column == -1.)[0][:]
            NfoundCB = len(indCB)
            indCT = np.where(column == 1.)[0][:]
            NfoundCT = len(indCT)
            CBarray[iTime, 0:NfoundCB] = height[indCB]
            CTarray[iTime, 0:NfoundCT] = height[indCT]
            NlayersArray[iTime] = numberCB[iTime]

    # calculating cloud thickness based on the cloud base and tops found ( 2d array (time, Nlevels))
    cloudThicknessDatabase = CTarray - CBarray
    # generating array of levels
    levels = np.arange(NCB)

    # step 2: build cloud database with all clouds saved as xarray dataset
    clouds = xr.Dataset(
        data_vars = {'cloudBase' : (('time', 'levels'), CBarray),
                     'cloudTop'  : (('time', 'levels'), CTarray),
                     'cloudThick': (('time', 'levels'), cloudThicknessDatabase)},
        coords    = {'levels': levels,
                     'time'  : time})

    # step 3: identify cloud properties of PBL clouds using timeStart, timeEnd, MaxCTheight
    cloudsTimeWindow = clouds.sel(time=slice(timeStart, timeEnd))

    
    PBLclouds = cloudsTimeWindow.where(cloudsTimeWindow.cloudTop < heightPBL)


    return(clouds, PBLclouds)

#--------------------------------------------------------------------------
def f_calculateMinCloudBaseTop(clouds, PBLclouds, date_arr):
    """author: claudia Acquistapace
     date: 18/05/2020
     goal: function to calculate the minimum cloud base for the PBL and the corresponding cloud top
     input: clouds - list of xarray datasets of cloud properties
            PBLclouds - list of xarray datasets of PBL cloud properties
            date_arr - array of days to be processed
    output: minimum cloud base and corresponding cloud tops in matrices of dimtime, Nfiles dimensions
    """
    
    # definition of output matrices
    dimTime = 9600
    Nfiles = len(date_arr)
    CBarr_obs = np.zeros((dimTime, Nfiles))
    CTarr_obs = np.zeros((dimTime, Nfiles))
    TKarr_obs = np.zeros((dimTime, Nfiles))
    CBarr_PBL_obs = np.zeros((dimTime, Nfiles))
    CTarr_PBL_obs = np.zeros((dimTime, Nfiles))
    TKarr_PBL_obs = np.zeros((dimTime, Nfiles))
    # for each day, reading and saving minimum cloud base and corresponding cloud top
    for indFile in range(Nfiles):
        
        # readingt the date
        date = date_arr[indFile]
        yy = int(date[0:4])
        mm = int(date[4:6])
        dd = int(date[6:8])
        timeStandard = pd.date_range(start=datetime.datetime(yy,mm,dd,0,0,0), \
                                    end=datetime.datetime(yy,mm,dd,23,59,59), freq='9s')
        # reading xarray datasets of the day
        PBLcloud_dataset = PBLclouds[indFile]
        cloud_dataset = clouds[indFile]

        PBLCloudsStandard = PBLcloud_dataset.reindex({'time':timeStandard})
        meanCB_obs = np.nanmin(cloud_dataset.cloudBase.values, axis=1)
        meanCT_obs = np.nanmin(cloud_dataset.cloudTop.values, axis=1)
        meanTK_obs = np.nanmin(cloud_dataset.cloudThick.values, axis=1)

        meanCB_obs = np.nanmin(cloud_dataset.cloudBase.values, axis=1)
        meanCT_obs = np.nanmin(cloud_dataset.cloudTop.values, axis=1)
        meanTK_obs = np.nanmin(cloud_dataset.cloudThick.values, axis=1)

        meanCB_PBL_obs = np.nanmin(PBLCloudsStandard.cloudBase.values, axis=1)
        meanCT_PBL_obs = np.nanmin(PBLCloudsStandard.cloudTop.values, axis=1)
        meanTK_PBL_obs = np.nanmin(PBLCloudsStandard.cloudThick.values, axis=1)

        CBarr_obs[:, indFile] = meanCB_obs
        CTarr_obs[:, indFile] = meanCT_obs
        TKarr_obs[:, indFile] = meanTK_obs
        CBarr_PBL_obs[:, indFile] = meanCB_PBL_obs
        CTarr_PBL_obs[:, indFile] = meanCT_PBL_obs
        TKarr_PBL_obs[:, indFile] = meanTK_PBL_obs
    return (CBarr_obs, CTarr_obs, TKarr_obs, CBarr_PBL_obs)



#---------------------------------------------------------------------------------


def f_resampleArrays2StandardData(A, index, strDate):
    """
    author : Claudia Acquistapace
    date   : 10/04/2020
    goal   : resample data with some misssing times (matrices of dimT < 9600) to the standard size (9600,150)
    input  : matrix of data to resize, datetime_array, height_array
    output : ndarray of resampled matrix with nans wherever missing data are located

    """
    import numpy as np
    import pandas as pd

    DF = pd.Series(A, index=index)
    # I first construct my regular time index every 9s
    # Obviously I put a little more periods (entire day would be 9600)
    index = pd.date_range(strDate, periods=9600, freq='9s')
    # There you go, by default missing values should be NaN
    DFresampled = DF.loc[index]
    return(DFresampled.values)
#---------------------------------------------------------------------------------


def f_resample2StandardData(A, index, cols, strDate):
    """
    author : Claudia Acquistapace
    date   : 10/04/2020
    goal   : resample data with some missing times (matrices of dimT < 9600) to the standard size (9600,150)
    input  : matrix of data to resize, datetime_array, height_array
    output : ndarray of resampled matrix with nans wherever missing data are located

    """
    import numpy as np
    import pandas as pd

    DF = pd.DataFrame(A, index=index, columns=cols)
    # I first construct my regular time index every 9s
    # Obviously I put a little more periods (entire day would be 9600)
    index = pd.date_range(strDate, periods=9600, freq='9s')
    # There you go, by default missing values should be NaN
    DFresampled = DF.loc[index]
    return(DFresampled.values)

    
#aa
# closest function
#---------------------------------------------------------------------------------
# date :  16.10.2017
# author: Claudia Acquistapace
# goal: return the index of the element of the input array that in closest to the value provided to the function
def f_closest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx



def getNearestIndex(timeRef, timeStamp):
    """this function finds the nearest element of timeRef array to the value timeStamp within the given tolerance 
    and returns the index of the element found. If non is found within the given tolerance, it returns nan."""
    try:
        index = timeRef.index.get_loc(timeStamp, method='nearest')
    
    except:

        index = np.nan
    
    return index

def getIndexList(dataTable, reference):
    """this function reads the less resolved time array (dataTable) and the time array to be used as reference (reference) 
    and the tolerance. then for every value in the reference array, it finds the index of the nearest element of 
    dataTable for a fixed tolerance. It provides as output the list of indeces of dataTable corresponding 
    to the closest elements of the reference array. """
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











# function to calculate LWC from 
#---------------------------------------------------------------------------------
# date :  23.10.2019
# author: Claudia Acquistapace
# goal: calculate LWC using standard Frisch approach 
# input: linear reflectivity matrix, radar range gate resolution (assumed constant), \
#time array, height attar, LWP time serie
# output: LWC matrix (time, height)
def f_calculateLWCFrisch(Ze_lin, deltaZ, datetime_ICON, height_ICON, LWP_obs_res):
        
    LWC_Frisch = np.zeros((len(datetime_ICON), len(height_ICON)))
    LWC_Frisch.fill(np.nan)
    LWP_obs_res = np.insert(LWP_obs_res,0,np.nan)
    for indT in range(len(datetime_ICON)):
        for indH in range(len(height_ICON)):
            num = LWP_obs_res[indT] * np.sqrt(Ze_lin[indT,indH])
            den = deltaZ*np.nansum(np.sqrt(Ze_lin[indT,:]))
            #print(den)
            #print(num)
            LWC_Frisch[indT,indH] = num/den
            
            
    return(LWC_Frisch)






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
        numberCB.append(len(np.where(column == 1.)[0][:]))
        numberCT.append(len(np.where(column ==-1.)[0][:]))  

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
        indCB=np.where(column == -1.)[0][:]
        NfoundCB=len(indCB)
        indCT=np.where(column == 1.)[0][:] 
        NfoundCT=len(indCT)
        CBarray[iTime,0:NfoundCB]=height[indCB]
        CTarray[iTime,0:NfoundCT]=height[indCT]
    
    
    
    return (CBarray,CTarray)





# function to define cb and ct of boundary layer clouds in a column appearing when reading from the top
#---------------------------------------------------------------------------------
# date :  22.10.2019
# author: Claudia Acquistapace
# input: cloudMask, dimension of time array, dimension of height array, height from model/obs (cloudnet)
# output: array of: 
#  - CBarray: time array having 4 dimensions to record four different cloud base heights per time stamp
#  - CTarray: time array having 4 dimensions to record four different cloud top heights per time stamp
#  - NlayersArray: number of distinct cloud layers identified per time stamp 
#  - CB_collective: minimum cloud base identified per time stamp
#  - CT_collective: maximum cloud top identified per time stamp
# goal: 

def f_calcCloudBaseTopPBLclouds(cloudMask, dimTime, dimHeight, height, cloudTimeArray, time):
    # cloud mask for identifying cloud base and cloud top of PBL clouds

    #CloudMaskCut = cloudMask
    # filtering clouds above 5000mt
    # ind_above = np.where(height > 5000.)
    #CloudMaskCut[:, ind_above] = 0.
    
    
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
        numberCB.append(len(np.where(column == 1.)[0][:]))
        numberCT.append(len(np.where(column ==-1.)[0][:]))  

    NCB=max(numberCB)   
    NCT=max(numberCT)
    
    # generating cloud base and cloud top arrays 
    CBarray      = np.zeros((dimTime,NCB))
    CBarray.fill(np.nan)
    CTarray      = np.zeros((dimTime,NCT))
    CTarray.fill(np.nan)
    NlayersArray = np.zeros((dimTime))
    NlayersArray.fill(np.nan)
    
    # if no cloud bases or no cloud tops are found, then CB and CT are assigned to nan
    if (NCB == 0) or (NCT == 0):
        CB_collective = np.zeros((dimTime))
        CB_collective.fill(np.nan)
        CT_collective = np.zeros((dimTime))
        CT_collective.fill(np.nan)
        CB_PBL_out = np.zeros((dimTime))
        CB_PBL_out.fill(np.nan)
        CT_PBL_out = np.zeros((dimTime))
        CT_PBL_out.fill(np.nan)
    else:
    # if some cloud base / cloud tops are found, all the found values are stored
        # storing cloud base and cloud top arrays
        for iTime in range(dimTime):
            column                    = gradBinary[iTime,:]
            indCB                     = np.where(column == -1.)[0][:]
            NfoundCB                  = len(indCB)
            indCT                     = np.where(column == 1.)[0][:] 
            NfoundCT                  = len(indCT)
            CBarray[iTime,0:NfoundCB] = height[indCB]
            CTarray[iTime,0:NfoundCT] = height[indCT]
            NlayersArray[iTime]       = numberCB[iTime]
        
        # we define a collective cloud base/top to consider multilayer PBL clouds as one
        # we assign min CB and max CT for each PBl cloud found.
        CB_collective = np.asarray(CBarray[:,0])
        CT_collective = np.asarray(CTarray[:,0])
        for ind in range(dimTime):
        #    if (np.isnan(CB[ind,0]) == True):
            CB_collective[ind] = np.nanmin(CBarray[ind,:])
            CT_collective[ind] = np.nanmax(CTarray[ind,:])


        # filtering clouds in PBL using human filtering for hours
        #if  np.count_nonzero(~np.isnan(CB_collective)) != 0:
        timeStart = cloudTimeArray[0]
        timeEnd   = cloudTimeArray[1]
        CB_PBL = pd.Series(np.repeat(np.nan, len(time)), index=time)
        maskt = (CB_PBL.index > timeStart) * (CB_PBL.index < timeEnd)
        CB_PBL.loc[maskt] = CB_collective[maskt]
        CT_PBL = pd.Series(np.repeat(np.nan, len(time)), index=time)
        maskt = (CT_PBL.index > timeStart) * (CT_PBL.index < timeEnd)
        CT_PBL.loc[maskt] = CT_collective[maskt]
        CT_PBL_out = CT_PBL.values
        CB_PBL_out = CB_PBL.values

    
    return (CBarray, CTarray, NlayersArray, CB_PBL_out, CT_PBL_out, CB_collective, CT_collective)




def f_calcCloudBaseTopPBLcloudsV2(cloudMask, dimTime, dimHeight, height, cloudTimeArray, \
                                  time):
    """
    @ author: cacquist
    @  date : 10 November 2019
    @  goal : this function corresponds to the version2 processing mode. It has been 
    generated to detect PBL clouds over JOYCE and it has been tuned with 
    statistical observed mean PBL cloud properties from the site. 
    INPUT:
        - cloudMask     : matrix of 0/1 containing cloud mask
        - dimTime       : dimension of time array
        - dimHeight     : dimension of height array
        - cloudTimeArry : 
        - time          : time array 
    OUTPUTS:
        - CBarray       : array containing all cloud bases found with the gradient method
        - CTarray       : array containing all cloud tops found with the gradient method
        - NlayersArray  : number fo cloud base/top found for each time
        - CB_PBL_out    : array of boundary layer cloud bases found 
        - CT_PBL_out    : array of boundary layer cloud tops found
        - CB_collective : array of minimum cloud base found
        - CT_collective : array of maximum cloud top found
    Methodology:
    It sets the cloud base to be below 2500mt and the cloud geometrical thickness
    to be 600 mt. 
    Check for cloud base height to be below 2500 mt:
        If cloud base does not fullfill the condition, no PBL cloud 
            base and top are found and it returns nans.
        If cloud base fullfills the condition, then it checks for cloud tops.
        If maximum cloud top is found above the CB + 600 mt, lower cloud
            tops are searched among the cloud tops below that height and the 
            minimum is taken.
        If none are found cloud top nd cloud base are assigned to nan.
    """
    meanCloudThickness = 600.
    minCBheight = 2500.
    # cloud mask for identifying cloud base and cloud top of PBL clouds

    # filtering clouds above 5000mt
    #cloudMaskCut = cloudMask
    #ind_above = np.where(height > 5000.)
    #cloudMaskCut[:, ind_above] = 0.
    
    
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
        numberCB.append(len(np.where(column == 1.)[0][:]))
        numberCT.append(len(np.where(column ==-1.)[0][:]))  

    NCB=max(numberCB)   
    NCT=max(numberCT)
    
    # generating cloud base and cloud top arrays 
    CBarray      = np.zeros((dimTime,NCB))
    CBarray.fill(np.nan)
    CTarray      = np.zeros((dimTime,NCT))
    CTarray.fill(np.nan)
    NlayersArray = np.zeros((dimTime))
    NlayersArray.fill(np.nan)
    
    # if no cloud bases or no cloud tops are found, then CB and CT are assigned to nan
    if (NCB == 0) or (NCT == 0):
        CB_collective = np.zeros((dimTime))
        CB_collective.fill(np.nan)
        CT_collective = np.zeros((dimTime))
        CT_collective.fill(np.nan)
        CB_PBL_out = np.zeros((dimTime))
        CB_PBL_out.fill(np.nan)
        CT_PBL_out = np.zeros((dimTime))
        CT_PBL_out.fill(np.nan)
    else:
    # if some cloud base / cloud tops are found, all the found values are stored
        # storing cloud base and cloud top arrays
        for iTime in range(dimTime):
            column                    = gradBinary[iTime,:]
            indCB                     = np.where(column == -1.)[0][:]
            NfoundCB                  = len(indCB)
            indCT                     = np.where(column == 1.)[0][:] 
            NfoundCT                  = len(indCT)
            CBarray[iTime,0:NfoundCB] = height[indCB]
            CTarray[iTime,0:NfoundCT] = height[indCT]
            NlayersArray[iTime]       = numberCB[iTime]
        
        
        # we define a collective cloud base/top to consider multilayer PBL clouds as one
        # we assign min CB and max CT for each PBl cloud found.
        CB_collective = np.asarray(CBarray[:,0])
        CT_collective = np.asarray(CTarray[:,0])
        CB_PBL_out    = np.repeat(np.nan, len(time))
        CT_PBL_out    = np.repeat(np.nan, len(time))
        
        for ind in range(dimTime):
        #    if (np.isnan(CB[ind,0]) == True):
            CB_collective[ind] = np.nanmin(CBarray[ind,:])
            CT_collective[ind] = np.nanmax(CTarray[ind,:])
            #selecting temporal window in which cloud top and base for PBL clouds have to be calculated
            if (time[ind] > cloudTimeArray[0]) * (time[ind] < cloudTimeArray[1]):
                if (CB_collective[ind] < minCBheight):
                    # for boundary layer clouds, we can assume the lowest cloud base is correct
                    # we can also assume that from the lowest cloud base, the cloud does not extend 
                    # in the vertical for more than 1500 mt. If the max cloud top is above 1500 mt then 
                    # we select among cloud tops, those that are located withing such distance from cloud base 
                    maxCTheightPBL = np.nanmin(CBarray[ind,:]) + meanCloudThickness
                    #print('max cloud top', maxCTheightPBL)
                    if (np.nanmax(CTarray[ind,:]) > maxCTheightPBL):
                        findLowerCT = np.where(CTarray[ind,:] < maxCTheightPBL)
                        if (len(findLowerCT[0]) == 0): # no elements are found below the maximum allowed height for cloud top
                            CT_PBL_out[ind] = np.nan
                            CB_PBL_out[ind] = np.nan
                        else:
                            #print('sono qui')
                            CT_PBL_out[ind] = np.nanmin(CTarray[ind,findLowerCT]) # assigning minmum cloud top
                            CB_PBL_out[ind] = CB_collective[ind] # assigning cloud base if it is below 2500 mt

    return (CBarray, CTarray, NlayersArray, CB_PBL_out, CT_PBL_out, CB_collective, CT_collective)


#---------------------------------------------------------------------------------
# date :  28.01.2019
# author: Claudia Acquistapace
# goal: function that calculates cloud fraction over 30 minutes of time for the whole day for ICON-LEM
#  input:
#  QI_ICON_LEM, \
#  QC_ICON_LEM, 
#  datetime_ICON, 
#  height_2_ICON_LEM, 
#  QiThreshold, 
#  QcThreshold
#    
# output: 
# mean_CF_liquid_ICON,
# mean_CF_ice_ICON, 
# mean_CF_tot_ICON, 
# datetime_out
#--------------------------------------------------------------------------------
    
def f_calculateCloudFractionICON(QI, QC, yy, mm, dd, time, height, QiThreshold, QcThreshold):
    # calculation of cloud fraction for ICON_LEM
    # creating a dataframe with for Qi and Qc with time and height using pandas dataframe
    
    QI_ICON_DF = pd.DataFrame(QI, index=time, columns=height)
    QC_ICON_DF = pd.DataFrame(QC, index=time, columns=height)

    # defining mean cloud fraction matrices to contain average profile every hour for the supersite
    mean_CF_liquid_ICON = np.zeros((48,150))
    mean_CF_ice_ICON    = np.zeros((48,150))
    mean_CF_tot_ICON    = np.zeros((48,150))
    deltaT                  = datetime.timedelta(minutes=30)
    indInt = 0
    datetime_out = []
    # --- loop on hours to calculate the mean hourly profile
    for itime in range(0,48):
        if indInt == 0:
            HourInf = datetime.datetime(int(yy), int(mm), int(dd), 0, 0, 0) 
        else:
            HourInf = HourInf + deltaT
        HourSup     = HourInf + deltaT
        datetime_out.append(HourInf)
        indInt = indInt + 1
    
        Qi_sliced_t = QI_ICON_DF.loc[(QI_ICON_DF.index < HourSup) * (QI_ICON_DF.index > HourInf),:]
        Qc_sliced_t = QC_ICON_DF.loc[(QC_ICON_DF.index < HourSup) * (QC_ICON_DF.index > HourInf),:]
        # ---- loop on heights: for each height counting the number of elements 
        # larger than the threshold and 
        # calculating the cloud fraction as the ratio between this number and 
        # the number of elements counted in the hour
        #print len(DF_qi_hour[DF_qi_hour.iloc[:,0] > QiThreshold])
        #print len(DF_qi_hour.iloc[:,0])
        
        for iheight in range(len(height)-1):
        #for iheight in range(2):   
            
            # extracting array
            DF_qi_arr = Qi_sliced_t.loc[:,height[iheight]]
            DF_qc_arr = Qc_sliced_t.loc[:,height[iheight]]
            NelemTot = len(DF_qi_arr)
            
            # posing conditions on cloud fraction for liquid only
            Cond_iceClouds=np.isfinite(DF_qi_arr[DF_qi_arr > QiThreshold] * DF_qc_arr[DF_qc_arr < QcThreshold])
            Cond_iceClouds.apply(int) 
            Num_iceCloud=Cond_iceClouds.sum()
            
            Cond_LiquidClouds=np.isfinite(DF_qc_arr[DF_qi_arr < QiThreshold] * DF_qc_arr[DF_qc_arr > QcThreshold])
            Cond_LiquidClouds.apply(int) 
            Num_liquidCloud=Cond_LiquidClouds.sum() 
            #print(Num_liquidCloud)
            #print(Num_iceCloud)        
            if float(NelemTot) == 0:
                print('Houston, we have a problem!')
            else:
                mean_CF_ice_ICON[itime,iheight]=float(Num_iceCloud)/float(NelemTot)
                mean_CF_liquid_ICON[itime,iheight]=float(Num_liquidCloud)/float(NelemTot)
                mean_CF_tot_ICON[itime,iheight]=float(Num_iceCloud+Num_liquidCloud)/float(NelemTot)
            
            # defining dictionary containing data to have as output
        dict_CF = {}
        # filling dictionaries with data 
        dict_CF = {
            'TotalCloudFraction':mean_CF_tot_ICON,
            'LiquidCloudFraction':mean_CF_liquid_ICON,
            'IceCloudFraction':mean_CF_ice_ICON,
            'height':height,
            'time':datetime_out,
                }
    
    return(dict_CF)


def f_plotTest(matrix, time, height, figname):
    pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/patch003/figures_JAMES/'
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    matplotlib.rc('xtick', labelsize=10)  # sets dimension of ticks in the plots
    matplotlib.rc('ytick', labelsize=10)  # sets dimension of ticks in the plots
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    cax = ax.pcolormesh(time, height, matrix, vmin=0, vmax=3,
                        cmap=plt.cm.get_cmap("GnBu", 4))
    ax.set_ylim(0., 15000)  # limits of the y-axes
    # ax.set_xlim(timeStart, timeEnd)  # limits of the x-axes
    ax.set_title("cloud mask model", fontsize=10)
    ax.set_xlabel("time [hh:mm]", fontsize=10)
    ax.set_ylabel("height [m]", fontsize=10)
    cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical', aspect=10)
    cbar.ticks = ([0, 1, 2, 3])
    cbar.ax.set_yticklabels(['no cloud', 'liquid', 'ice', 'mixed phase'])
    cbar.set_label(label="cloud type", size=10)
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(pathFig + figname + '_cloudMask.png')
    return ()




    
""" function to derive wind speed and direction
#---------------------------------------------------------------------------------
 date :  17.12.2018
 author: Claudia Acquistapace (cacquist@meteo.uni-koeln.de)
 goal: derive wind speed and direction in form of list and matrices
 input: 
 - datetime_ICON: time array
 - height_ICON: height array
 - u_ms: zonal wind
 - v_ms: meridional wind
 output:
 - ws: list of wind speed
 - wd: list of wind directions
 - wind_abs: matrix of wind speed
 - wind_dir_trig_from_degrees: matrix of wind direction in degrees indicating
 the direction from where wind is coming
"""
#--------------------------------------------------------------------------------
def f_calcWindSpeed_Dir(datetime_ICON, height_ICON, u_ms, v_ms):
    import math
    wind_abs                   = np.sqrt(u_ms**2 + v_ms**2)
    wind_dir_trig_to           = np.zeros((len(datetime_ICON),len(height_ICON)))
    wind_dir_trig_to_degrees   = np.zeros((len(datetime_ICON),len(height_ICON)))
    wind_dir_trig_from_degrees = np.zeros((len(datetime_ICON),len(height_ICON)))
    wind_dir_cardinal          = np.zeros((len(datetime_ICON),len(height_ICON)))
    ws                         = []
    wd                         = []
    for itime in range(len(datetime_ICON)):
        for iHeight in range(len(height_ICON)):
            # wind dir in unit circle coordinates (wind_dir_trig_to), which increase counterclockwise and have a zero on the x-axis
            wind_dir_trig_to[itime, iHeight]           = math.atan2(v_ms[itime, iHeight],u_ms[itime, iHeight]) 
            # wind dir in degrees (wind_dir_trig_to_degrees) dir where wind goes
            wind_dir_trig_to_degrees[itime, iHeight]   = wind_dir_trig_to[itime, iHeight] * 180/math.pi ## -111.6 degrees  
            # wind dir in degrees (wind_dir_trig_to_degrees) dir from where wind comes
            wind_dir_trig_from_degrees[itime, iHeight] = wind_dir_trig_to_degrees[itime, iHeight] + 180 ## 68.38 degrees
            # wind dir in cardinal coordinates from the wind dir in degrees (wind_dir_trig_to_degrees) dir from where wind comes
            wind_dir_cardinal[itime, iHeight]          = 90 - wind_dir_trig_from_degrees[itime, iHeight]
            if np.isfinite(wind_dir_trig_from_degrees[itime, iHeight]) and \
                np.isfinite(wind_abs[itime, iHeight]) and \
                (wind_abs[itime, iHeight] != 0.):
                wd.append(wind_dir_trig_from_degrees[itime, iHeight])
                ws.append(wind_abs[itime, iHeight])
    WindDictionary={'windDirection':ws,
                    'windSpeed':wd,
                    }
    return(WindDictionary)




# function to plot color time heigth maps from a dictionary of initial data
#---------------------------------------------------------------------------------
# date :  14.12.2018
# author: Claudia Acquistapace (cacquist@meteo.uni-koeln.de)
# goal: function to derive pdfs of vertical and horizontal wind below cloud base
# check for vertical wind values observed below cloud base. for every time stamp. 
# methodology: for observations: 
# if there is cloud base in observations, store vertical wind values recorded in the 300m below cloud base.
# if there is no cloud, store vertical wind values in the 5 bins below mean estimated cloud base.
# input : vertical wind, horizontal wind, time, height
# output: verticalWindPDF_cloud, verticalWindPDF_nocloud, horizontalWindPDF_cloud, horizontalWindPDF_nocloud
#--------------------------------------------------------------------------------
def f_pdfsBelowCloudBase(w_ICON, Hwind, height, datetime_ICON, datetimeHourArr, height_ICON, mean_CB_arr_OBS, CB_array_OBS, timeStart, timeEnd):

    verticalWindPDF_cloud         = []
    horizontalWindPDF_cloud       = []
    verticalWindPDF_nocloud       = []
    horizontalWindPDF_nocloud     = []

    distHeight                    = 400.
    vertWind_ICON_DF              = pd.DataFrame(w_ICON, index=datetime_ICON, columns=height)
    HorWind_ICON_DF               = pd.DataFrame(Hwind, index=datetime_ICON, columns=height_ICON)
    limTimeInf                    = timeStart
    limTimeSup                    = timeEnd
    
    # establishing height below which to check for wind
    for indTime in range(len(datetime_ICON)):
        if (datetime_ICON[indTime] > limTimeInf) * (datetime_ICON[indTime] < limTimeSup):
            
            # case of no clouds, read mean cloud base height in the hour and extract height
            if np.isfinite(CB_array_OBS[indTime]) == False:
                findHourInd           = f_closest(np.asarray(datetimeHourArr), datetime_ICON[indTime])
                CBHeight              = mean_CB_arr_OBS[findHourInd]

                mask_h_vertWind           = (vertWind_ICON_DF.columns < CBHeight) * (vertWind_ICON_DF.columns > CBHeight-distHeight)
                valuesWwind               =  vertWind_ICON_DF.values[indTime, mask_h_vertWind].flatten()
                mask_h_horwind            = (HorWind_ICON_DF.columns < CBHeight) * (HorWind_ICON_DF.columns > CBHeight-distHeight)
                valuesHwind               =  HorWind_ICON_DF.values[indTime, mask_h_horwind].flatten()        
                for indValw in range(len(valuesWwind)):
                    verticalWindPDF_nocloud.append(valuesWwind[indValw])
                for indValh in range(len(valuesHwind)):
                    horizontalWindPDF_nocloud.append(valuesHwind[indValh])



            # case of clouds: read cloud base height and extract bins below.
            else:
                CBHeight              = CB_array_OBS[indTime]

                mask_h_vertWind           = (vertWind_ICON_DF.columns < CBHeight) * (vertWind_ICON_DF.columns > CBHeight-distHeight)
                valuesWwind               =  vertWind_ICON_DF.values[indTime, mask_h_vertWind].flatten()
                mask_h_horwind            = (HorWind_ICON_DF.columns < CBHeight) * (HorWind_ICON_DF.columns > CBHeight-distHeight)
                valuesHwind               =  HorWind_ICON_DF.values[indTime, mask_h_horwind].flatten() 

                for indValw in range(len(valuesWwind)):
                    verticalWindPDF_cloud.append(valuesWwind[indValw])
                for indValh in range(len(valuesHwind)):
                    horizontalWindPDF_cloud.append(valuesHwind[indValh])    
        
    return(verticalWindPDF_cloud, verticalWindPDF_nocloud, horizontalWindPDF_cloud, horizontalWindPDF_nocloud) 


def f_calcPblHeightRN(thetaV,Uwind,Vwind,height,time,device):
    """
    PBL height calculation function
    --------------------------------------------------------------------------------
    date created :  15.01.2018
    date modifed :  05.12.2019
    author: Claudia Acquistapace
    goal: calculate the boundary layer height following the richardson number
    derivation according to Seidel Et al, 2010
    #---------------------------------------------------------------------------------
    """
    g            = 9.8        # gravity constant
    Rithreshold  = 0.25       # Threshold values for Ri
    #Rithreshold2 = 0.2
    dimTime      = len(time)
    dimHeight    = len(height)
    if (device == 'mod'):
        zs       = height[149]                                        # height of the surface reference
    if (device == 'obs'):
        zs       = height[0]
    RiMatrix     = np.zeros((dimTime, dimHeight))                    # Richardson number matrix
    PBLheightArr = []
    RiCol        = np.zeros((dimHeight))
    
    # calculating richardson number matrix
    for iTime in range(dimTime):
        thetaS = thetaV[iTime,149]
        for iHeight in range(dimHeight):
            den = ((Uwind[iTime,iHeight])**2 + (Vwind[iTime,iHeight])**2)
            if den == 0.:
                RiMatrix[iTime,iHeight] = 0.
            else:
                RiMatrix[iTime,iHeight] = (1/den) * (g/thetaS) * (thetaV[iTime,iHeight]-thetaS)*(height[iHeight]-zs)
            
          
    # find index in height where Ri > Rithreshold
    for iTime in range(dimTime):
        RiCol=RiMatrix[iTime,:]
        #print(RiCol)
        #print(np.where(RiCol > Rithreshold2)[0][:])
        #print(len(np.where(RiCol > Rithreshold)[0][:]))
        if len(np.where(RiCol > Rithreshold)[0][:]) != 0:
            PBLheightArr.append(height[np.where(RiCol > Rithreshold)[0][-1]] - height[dimHeight-1])
        else:
            PBLheightArr.append(0)
    return PBLheightArr


def f_calcPblHeightTW(stdWmatrix,sigmaThreshold,height2,time, device):
    """
    PBL height calculation function based on threshold on std w method
    --------------------------------------------------------------------------------
    date created :  05.12.2019
    author: Claudia Acquistapace
    goal: calculate the boundary layer height following the method of a threshold on sigma w
    as indicated in Schween et al., 2014. The algorithm takes the maximum of the heights below 2000m
    at which the sigma values is larger than 0.4. 2000m is a conservative value
    threshold obtained from the paper from Schween et al., 2014 on MLH at JOYCE 
    #---------------------------------------------------------------------------------
    """
    dimTime     = len(time)
    PBLheightTW = np.zeros((dimTime))

    PBLheightTW.fill(np.nan)

    #std_matrix[:,height < height[142]] = 0.
    for ind in range(len(time)):
        if device == 'mod':
            column = stdWmatrix[ind,:]
            aboveThr = column > sigmaThreshold
            
            #selecting heights below 2000
            Hsel = height2[aboveThr]
            Hbelow = Hsel[Hsel < 2000.]
            if np.count_nonzero((Hbelow)) != 0:
                PBLheightTW[ind] = np.nanmax(Hbelow)
                

            
    return(PBLheightTW)
    
# function to calculate the convective condensation level height and temperature
#---------------------------------------------------------------------------------
# date :  17.05.2018
# author: Claudia Acquistapace
# goal: function that calculates the convective condensation level (CCL) height and temperature. 
# for the definition of the CCL check this: https://en.wikipedia.org/wiki/Convective_condensation_level
# input: 
# - T field (time, height)
# - RH field (time, height) ( es: 75,0)
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
    # temperature has to be provided in K
    # RH in % ( 70.14)
    # P In Kpa
    dimHeight = len(height_ICON) 
    # defining constants
    cost_rvl = np.power(5423,-1.) #K
    E0       = 0.611 # Kpa
    T0       = 273.  # K 
    Rv       = 461   # J K^-1 Kg^-1
    epsilon  = 0.622 
    Ad_rate  = -9.8  # K/Km

    # ---- substituting RH = 0. to RH = nan to avoid log(0) cases
    RH_ICON [ RH_ICON == 0.] = np.nan
    T_ICON [ T_ICON == 0.]   = np.nan

    # ---- calculating due point temperature profile for each time (The dew point is \
    # the temperature to which air must be cooled to become saturated with water vapor. )
    Td = np.power(np.power(T_ICON,-1.)-cost_rvl*np.log(RH_ICON/100.),-1.)

    # ---- calculating mixing ratio at the surface for T dew point 
    Td_surf = Td[:, dimHeight-1]
    P_surf = P_ICON[:,dimHeight-1]
    RH_surf = RH_ICON[:,dimHeight-1]
    
    for indtime in range(len(datetime_ICON)):
        if (~np.isfinite(RH_surf[indtime])):
            RH_surf[indtime] = RH_ICON[indtime,dimHeight-2]
        if (~np.isfinite(Td_surf[indtime])):
            Td_surf[indtime] = Td[indtime,dimHeight-2]
            

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
    #print(dimHeight)
    for indTime in range(len(datetime_ICON)):
        for indHeight in range(dimHeight-2,1,-1):
            #print(height_ICON[indHeight])
            if (m[indTime,indHeight] < M0[indTime] and m[indTime,indHeight-1] > M0[indTime] ):
                ind_CCL.append(indHeight)
                break
            if indHeight == 1:
                ind_CCL.append(dimHeight-1)
    print(len(ind_CCL))
    z_ccl = height_ICON[ind_CCL]

    
    # ---- finding z(CCL) using the dry adiabatic lapse rate
    T_ground_CCL = []
    for indTime in range(len(ind_CCL)):
        T_top = T_ICON[indTime, ind_CCL[indTime]]
        T_ground_CCL.append(T_top - Ad_rate* z_ccl[indTime]*10.**(-3))
    
    dict_out={'z_ccl':z_ccl, 
              'T_ccl':T_ground_CCL, 
              'Td':Td
    }
    return(dict_out)


def f_CCL_new(T, P, RH, height, time, date):
    """
    function to calculate convective condensation level (CCL). For more info on definitions of this level, read pp.250
    of Petty : A first course in atmospheric thermodynamics
    input: T: temperature , to be provided in K
              relative humidity, in % (es: 70.14)
              pressure, in Kpa
              device, string for "model" or "obs"
    procedure:
        step 1: calculate dew point T
        step 2: calculate saturation mixing ratio m0 at t=Td, P=Psurf
        step 3: calculate, for every value of P, Td(m=m0, P)
        step 4: check, for every level of P, if there's a level i  of P for which  T(P)i-1 < Td(m0,P)i < T(P)i+1.
            If the level is found, assign T_star = T(P)i and Z_ccl as the height corresponding to that pressure height.
        step 5: calculate Tc using adiabatic lapse rate to come back at the height of the surface.
    output: T_ccl, z_ccl

    """
    #pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/' + patch + '/figures_JAMES/debugging/'

    print('calculating CCL height and T_CCL')
    # defining constants
    cost_rvl = np.power(5423, -1.)  # K
    E0 = 0.611  # Kpa
    T0 = 273.  # K
    Rv = 461  # J K^-1 Kg^-1
    L = 5.6 * 10 ** 6  # J/Kg
    epsilon = 0.622
    Ad_rate = -9.8  # K/Km

    # assigning dimensions:
    dimHeight = len(height)
    dimTime = len(time)

    # step 1: calculating due point temperature profile for each time (The dew point is \
    # the temperature to which air must be cooled to become saturated with water vapor. )
    # substituting RH = 0. to RH = nan to avoid log(0) cases
    RH[RH == 0.] = np.nan
    T[T == 0.] = np.nan
    # calculating Td
    Td = np.power(np.power(T, -1.) - cost_rvl * np.log(RH / 100.), -1.)

    # step 2: calculating mixing ratio at the surface for T = Td and P=Psurf
    # finding index of height corresponding to lowest level in height
    indHmin = np.nanargmin((height))

    # reading values of P, T, RH at the corresponding height
    Td_surf = Td[:, indHmin]
    P_surf = P[:, indHmin]
    RH_surf = RH[:, indHmin]
    m0 = epsilon * E0 * np.exp((1. / Rv) * (T0 ** (-1.) - Td_surf ** (-1.))) / (
                P_surf - E0 * np.exp((1. / Rv) * (T0 ** (-1.) - Td_surf ** (-1.))))

    #print(Td_surf, P_surf, RH_surf, m0)
    # step 3: calculating Td(m=m0, P) for every P value
    z_ccl = np.zeros((dimTime))
    T_cclTop = np.zeros((dimTime))
    z_ccl.fill(np.nan)
    # indPlotCount = 0
    for indTime in range(dimTime):
        Tdm0_profile = np.zeros((dimHeight))
        Tdm0_profile.fill(np.nan)
        indCCLprofile = []
        Tm0_surface = 1 / ((1 / T0) - ((1 / L) * Rv * np.log((m0 * P[:, indHmin]) / (E0 * epsilon))))
        for indHeight in range(dimHeight - 1):
            Tdm0_profile[indHeight] = 1 / (
                        (1 / T0) - ((1 / L) * Rv * np.log((m0[indTime] * P[indTime, indHeight]) / (E0 * epsilon))))
            # print(T[indTime, indHmin])

            if (T[indTime, indHeight] < Tdm0_profile[indHeight]) and (
                    T[indTime, indHeight + 1] > Tdm0_profile[indHeight]):
                indCCLprofile.append(indHeight)
                # print(Tdm0_profile[indHmin])
                # print(T[indTime, indHmin])

        #print(indCCLprofile)
            ##fig, ax = plt.subplots(figsize=(12, 5))
            # plt.plot(Tdm0_profile, height, label='TDm0')
            # plt.plot(T[indTime, :], height, label='T')
            # plt.legend()
            # plt.plot(time, z_ccl3)
            # plt.plot(np.repeat(M0[5000],len(height)), height)
            # plt.ylim(0, 6000)
            # plt.savefig(pathFig + str(indPlotCount) + 'Tm0_profile_Check.png', format='png')
            # indPlotCount = indPlotCount +1
        # print(len(indCCLprofile))
        if len(indCCLprofile) == 0:
            z_ccl[indTime] = np.nan
            T_cclTop[indTime] = np.nan
        else:
            z_ccl[indTime] = np.nanmin(height[indCCLprofile])
            T_cclTop[indTime] = np.nanmin(T[indTime, np.nanargmin(height[indCCLprofile])])

    # fig, ax = plt.subplots(figsize=(12,5))
    # plt.plot(time, z_ccl)
    # plt.ylim(0,6000)
    # plt.savefig(pathFig+date+'_z_ccl_mod.png', format='png')
    print(z_ccl)
    # ---- finding z(CCL) using the dry adiabatic lapse rate
    T_ground_CCL = np.zeros((dimTime))
    for indTime in range(dimTime):
        T_ground_CCL[indTime] = (T_cclTop[indTime] - Ad_rate * z_ccl[indTime] * 10. ** (-3))

    # providing output as standardized xarray output format
    #DatasetOut = xr.Dataset(
    #    data_vars={'z_ccl'   : (('time'), z_ccl),
    #               't_ccltop': (('time'), T_cclTop),
    #               't_ccl'   : (('time'), T_ground_CCL),
    #               'T_dew'   : (('time', 'height'), Td)},
    #    coords={'time'  : time,
    #            'height': height})

    #return (DatasetOut)

    DatasetOut = {'time':time,
                  'height':height,
                  'z_ccl':z_ccl,
                  'T_ground_ccl':T_ground_CCL,
                  'T_top_ccl':T_cclTop,
                  'T_dew':Td}
    return (DatasetOut)


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
# input: 
# - matrix of vertical velocity
# - time array
# - height array
# - time window for the running mean (30 min for comparing to obs)
#--------------------------------------------------------------------------------
def f_calcWvariance(Wwind,time,height,window, res):
    """
    OBSOLETE FUNCTION NOT USED ANYMORE
    author: claudia acquistapace
    date: 05.12.2019
    goal: calculation of variance of vertical velocity. The function performs the
    calculation of the standard deviation of vertical velocity as in the paper 
    from Schween et al., 2014., AMT, doi:10.5194/amt-7-3685-2014
    
    input: 
        - Wwind: vertical velocity matrix (time, height)
        - time: time array
        - height: height array
        - window: time window over which to calculate the standard deviation
        - res: resolution at which calculate the standard deviation (in the 
                                                                     paper it is 5 min)
    """
    
    dimTime   = len(time)
    dimHeight = len(height)
    variance  = np.zeros((dimTime,dimHeight))
    
    for iheight in range(dimHeight):
        
        # reading array of w values at a given height
        Warray = Wwind[:,iheight]
        s = pd.Series(Warray)
        # calculating running mean variance over the selected array
        # variance[:,iheight] = runningMeanVariance(Warray,window)\
        #variance[:,iheight] = np.power(pd.rolling_std(Warray, window),2) 
        variance[:,iheight] = np.power(s.rolling(window).std(),2)
    return variance
        
    
# skewness of vertical velocity calculation
#---------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
def f_runningMeanSkewnessVarianceStd_W(time, timeWindowSk, runningWindow, height, vertWind):
    """ 
        author: Claudia Acquistapace
        date created :  17.01.2018
        date modified: 05.12.2019
        goal : calculate running mean skewness and standard deviation of vertical velocity over a time window
        given. The skewness is calculated from the data of the surroundings +- timewindow/2 min
        The mean skewness is calculated on a resolution of timeWindowSk. 
        Processing of the data follows what indicated in Schween et al., 2014, AMT
        DOI: 10.5194/amt-703685-2014
        input parameters:
            - time: time array
            - timeWindowSk: time resolution on which to calculate the skewness (5 min in the paper)
            - runningWindow: time window on which to calculate the running mean
            - height: height array
            - vertWind: matrix of w velocity (time, height)
        output: 
            - skewness_w matrix (time, height)
    """
    dimTime    = len(time)
    dimHeight  = len(height)
    SKmatrix   = np.zeros((dimTime, dimHeight))
    stdWmatrix = np.zeros((dimTime,dimHeight))
    varianceW   = np.zeros((dimTime,dimHeight))
    SKmatrix.fill(np.nan)
    stdWmatrix.fill(np.nan)
    varianceW.fill(np.nan)
    
    for iTime in range(0, dimTime-1, timeWindowSk):

        # condition to skip first and last half interval of time stamps not surrounded 
        # by values
        if (iTime > runningWindow/2) & (iTime < dimTime-1-runningWindow/2):
            
            # generating indeces to read corresponding elements in the matrix
            timeIndeces = np.arange(iTime-runningWindow/2, iTime+runningWindow/2, dtype=np.int16)
            #print(iTime)
            
            for iHeight in range(dimHeight):
            
                meanW    = np.mean(vertWind[timeIndeces, iHeight])            # mean of the wind array
                variance = np.var(vertWind[timeIndeces, iHeight])         # variance of the wind array
                
                wprime   = np.subtract(vertWind[timeIndeces, iHeight],np.tile(meanW,len(timeIndeces)))
                num      = np.mean(np.power(wprime,3))
                den      = variance**(3./2.)
                
                stdWmatrix[iTime:iTime+timeWindowSk-1, iHeight] = np.sqrt(variance)
                SKmatrix[iTime:iTime+timeWindowSk-1, iHeight]   = num/den
                varianceW[iTime:iTime+timeWindowSk-1, iHeight]  = variance
                
    
    return (varianceW, stdWmatrix, SKmatrix)




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
    #                         time, \
    #                         height2, \
    #                         gradWindThr, \
    #                         SigmaWThres, \
    #                         ylim, \
    #                         cloudMask, \
    #                         varianceWmatrix, \
    #                         SKmatrix, \
    #                         stabilityArr, \
    #                        connection2Surface, \
    #                         shear_ICON, \
    #                         CB_array_ICON)
    dimTime           = len(time)
    dimHeight         = len(height)
    PBLclass          = np.zeros((dimTime, dimHeight))   # defining output matrix
    shear             = gradWindspeed#.transpose()           # transposing gradWindspeed to conform to other matrices

    # defining flag matrices and filling them with nan values. Each flag corresponds to a check to be performed for the classification. Checks are for:  cloud/turbulence/cloud driven/unstable at surface/wind shear/surface driven
    flagCloud         = np.zeros((dimTime, dimHeight))
    flagCloud.fill(np.nan)
    flagTurb          = np.zeros((dimTime, dimHeight))
    flagTurb.fill(np.nan)
    flagcloudDriven   = np.zeros((dimTime, dimHeight))
    flagcloudDriven.fill(np.nan)
    flagInstability   = np.zeros((dimTime, dimHeight))
    flagInstability.fill(np.nan)
    flagConnection    = np.zeros((dimTime, dimHeight))
    flagConnection.fill(np.nan)
    flagSurfaceDriven = np.zeros((dimTime, dimHeight))
    flagSurfaceDriven.fill(np.nan)
    flagWindShear     = np.zeros((dimTime, dimHeight))
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
    # function to calculate lcl height providing T in K and P in pascal, RH between 0 and 1.
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



def f_resamplingfield(FieldStart,datetimeFieldStart,ICON_DF):
    # goal: function to resample the 1d time seried of the scalar field 1dFieldStart 
    # having dimension equal to the length of the array datetime1dFieldStart,
    # to the time resolution provided in the dataframe ICON_DF, which is defined
    # externally and provided to the function as input. ICON_DF represent the 
    # data resolution we want to obtain with this resampling. An example of the 
    # definition of such dataframe is:
    # ICON_DF = pd.DataFrame(cloudMask_ICON, index = datetime_ICON, columns = height_ICON)
    # data: 22 July 2019
    # author: Claudia Acquistapace
    # output: dataframe of the same time res of the desired res
    
    import pandas as pd
        
    FieldStart_DF   = pd.DataFrame(FieldStart, index=datetimeFieldStart)
    SelectedIndex   = getIndexListsuka(FieldStart_DF, ICON_DF.index)
    values          = np.arange(0, len(ICON_DF.index))
    Field_resampled = pd.DataFrame(values, index=ICON_DF.index)
    Field_resampled = getResampledDataPd(Field_resampled, FieldStart_DF, SelectedIndex)
        
    return(Field_resampled)
        
    
def f_resampling_twoD_Field(FieldStart,datetimeFieldStart,heightFieldStart, ICON_DF, ICON_DF_T):
    # goal: resample the 2d field on the grid indicated by the second field 
    # provided in ICON_DF. The code first resamples in time and then it resamples in 
    # height applying the same metholodogy but to reversed arrays
    # definition of such dataframe is:
    # ICON_DF = pd.DataFrame(cloudMask_ICON, index = datetime_ICON, columns = height_ICON)
    # data: 22 July 2019
    # author: Claudia Acquistapace 
    # output: dataframe of the same time and height desired 
    
    import pandas as pd
    
    FieldStart_DF = pd.DataFrame(FieldStart, index=datetimeFieldStart, columns=heightFieldStart)
    SelectedIndex = getIndexList(FieldStart_DF, ICON_DF.index)
    values        = np.zeros((len(ICON_DF.index),len(heightFieldStart)))
    FieldResampled = pd.DataFrame(values, index=ICON_DF.index, columns=heightFieldStart)
    FieldResampled = getResampledDataPd(FieldResampled, FieldStart_DF, SelectedIndex)
    
    FieldTransposed = FieldResampled.values.transpose()
    FieldTrans_DF   = pd.DataFrame(FieldTransposed, index=heightFieldStart, columns=ICON_DF.index)
    SelectedCols    = getIndexList(FieldTrans_DF, ICON_DF_T.index)
    valuesTrans     = np.zeros((len(ICON_DF.columns), len(ICON_DF.index)))
    FieldTranspLess = pd.DataFrame(FieldTransposed, index=heightFieldStart, columns=ICON_DF.index)
    FieldFinalRes   = pd.DataFrame(valuesTrans, index=ICON_DF.columns, columns=ICON_DF.index)
    FieldFinalRes   = getResampledDataPd(FieldFinalRes, FieldTranspLess, SelectedCols)

    return(FieldFinalRes)
    
def f_downscaleScalarfield(scalarField, datetimehighRes, ICON_DF):
    # goal: function to downscale the resolution of fields having higher 
    # resolutions compared to ICON_DF
    # date: 23 July 2019
    # author: claudia Acquistapace (cacquist@meteo.uni-koeln.de)
    # output: pandas dataframe of the same scalar field downscaled
    dimLowRes      = len(ICON_DF.index)
    datetimeLowRes = ICON_DF.index
    scalarField[scalarField == -99.] = np.nan
    scalarField_DF = pd.DataFrame(scalarField, index=datetimehighRes) 
    lowResField = []
            
    for iTime in range(dimLowRes-1):
        timeStart = datetimeLowRes[iTime]
        timeEnd   = datetimeLowRes[iTime+1]
        indSel = (scalarField_DF.index < timeEnd) * (scalarField_DF.index > timeStart)
        #print(indSel)
        lowResField.append(np.nanmean(scalarField_DF.values[indSel]))
    return(np.array(lowResField))
    
    
    
    
def f_downscalevectorfield(vectorField, datetimehighRes, heightField, ICON_DF):
    # goal: function to downscale the resolution of fields having higher 
    # resolutions compared to ICON_DF
    # date: 23 July 2019
    # author: claudia Acquistapace (cacquist@meteo.uni-koeln.de)
    # output: array of the same scalar field downscaled
    dimLowRes      = len(ICON_DF.index)
    datetimeLowRes = ICON_DF.index
    vectorField[vectorField == -99.] = np.nan
    vectorField_DF = pd.DataFrame(vectorField, index=datetimehighRes, columns=heightField ) 
    lowResVectorField = pd.DataFrame(np.zeros((len(ICON_DF.index), len(heightField))), index=ICON_DF.index,\
                                     columns=heightField ) 
            
    for iTime in range(dimLowRes-1):
        timeStart = datetimeLowRes[iTime]
        timeEnd   = datetimeLowRes[iTime+1]
        maskT = (vectorField_DF.index < timeEnd) * (vectorField_DF.index > timeStart)
        #print(indSel)
        lowResVectorField.values[iTime,:] = np.nanmean(vectorField_DF.values[maskT, :], axis=0)
    return(lowResVectorField)
    
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
    cols2                  = ['Zeit [min:sec]','P [hPa]','T [C]','U [%]','Wind speed [m/s]','Wdir [inf]','Lange [inf]'\
                              ,'Breite [inf]','Hˆhe [m]','Geo Pot [m]', 'dew [C]', 'Tv [C]','Rs [m/s]', 'D [kg/m3]']
    Nfiles                 = len(fileList)
    RadiosondesData       = []
    for iFile in range(Nfiles):
        dayFile            = fileList[iFile][44:52]
        hourFile           = int(fileList[iFile][52:54])
        DatetimeRadiosonde = datetime.datetime(int(yy), int(mm), int(dd), hourFile, 0, 0)
        print(dayFile, hourFile)

        if dayFile == '20130414':
            #print('def correct for 14 columns')
            DF                 = pd.read_csv(fileList[iFile], sep='\t', skipinitialspace=True, \
                                         encoding='latin_1', names=cols2, header=0, dtype={'Hˆhe [m]':str})
        else:
            #print('def correct for 17 columns')
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
        e0                 = 0.611 # Kpa
        T0                 = 273.15 # K
        cost_lrv           = 5423. # K (Lv/Rv)
        L                  = 5.6 * 10 ** 6  # J/Kg

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
        P_surf             = P[0]
        T_surf             = T[0]

        # step 1: calculating due point temperature profile for each time was done at the previous step
        # step 2: calculating mixing ratio M0 at the surface for T = Td and P=Psurf
        M0                 = epsilon*e0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))) / \
        (P_surf*0.1 - e0*np.exp((1./Rv)*(T0**(-1.)-Td_surf**(-1.))))


        # step 3: calculating Td(m=m0, P) for every P value and assigning Z_ccl when there's a level i  of P for which
        # T(P)i-1 < Td(m0,P)i < T(P)i+1. If the level is found, assign T_star = T(P)i and Z_ccl as the height
        # corresponding to that pressure height.
        dimHeight = len(height)
        Tdm0_profile = np.zeros((dimHeight))
        Tdm0_profile.fill(np.nan)
        indCCLprofile = []
        for indHeight in range(dimHeight - 1):
            Pkpa = P[indHeight] * 0.1
            Tdm0_profile[indHeight] = 1 / (
                            (1 / T0) - ((1 / L) * Rv * np.log((M0 * Pkpa) / (e0 * epsilon))))
            if T[0] > Tdm0_profile[indHeight]:
                if (T[indHeight] > Tdm0_profile[indHeight]) and (
                            T[indHeight + 1] < Tdm0_profile[indHeight]):
                    indCCLprofile.append(indHeight)
            else:
                if (T[indHeight] < Tdm0_profile[indHeight]) and (
                            T[indHeight + 1] > Tdm0_profile[indHeight]):
                    indCCLprofile.append(indHeight)
        #pathFig = '/work/cacquist/HDCP2_S2/statistics/figs/patch003/figures_JAMES/debugging/radiosondes/'
        #fig, ax = plt.subplots(figsize=(12, 5))
        #plt.plot(Tdm0_profile, height, label='TDm0')
        #plt.plot(T[:], height, label='T')
        #plt.legend()
        #plt.ylim(0, 6000)
        #plt.savefig(pathFig + str(dayFile) + '_' + str(hourFile) +  'radiosonde_Check.png', format='png')
        #print(len(indCCLprofile))
        if len(indCCLprofile) == 0:
            z_ccl = np.nan
            T_cclTop = np.nan
        else:
            z_ccl = np.nanmin(height[indCCLprofile])
            T_cclTop = np.nanmin(T[np.nanargmin(height[indCCLprofile])])

        Ad_rate            = -9.8 # K/Km
        T_ground_CCL       = []
        # ---- finding z(CCL) using the dry adiabatic lapse rate
        T_ground_CCL       = float(T_cclTop) - Ad_rate* float(z_ccl)*10.**(-3)
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

        # ---- calculating mixing ratio
        mr =[]
        for indHeight in range(len(height)):
            mr.append((0.622*e[indHeight]*100.)/(P[indHeight]*100.-e[indHeight]*100.)) # water vapor mixing ratio in kg/kg
       
       
        # ---- calculating saturation mixing ratio
        ws = [] #mpcalc.saturation_mixing_ratio(P, T)
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
        gamma_moist_atmos  = atmos.equations.Gammam_from_rvs_T(Ws_array.astype(float), T_array.astype(float))
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

        # averaging lowest 100 mt of Ri matrix values
        indHeightMean = np.where(height < 100.)
        RImeanval = np.nanmean(RiCol[indHeightMean])
        RiCol[indHeightMean] = np.repeat(RImeanval, len(indHeightMean))
        
        #print(RiCol)
        #print(np.where(RiCol > Rithreshold2)[0][:])
        #print(len(np.where(RiCol > Rithreshold)[0][:]))
        if len(np.where(RiCol > Rithreshold)[0][:]) != 0:
            PBLheight = (height[np.where(RiCol > Rithreshold)[0][0]] - height[0])
        else:
            PBLheight = 0.

        # ---- saving variables in dictionary: every dictionary for one hour
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
        # appending  to a list: the list is for the day
        RadiosondesData.append(dict_day)
    return(RadiosondesData)        
    
def f_calcThermodynamics(P,Q,T, LTS, time, height, Hsurf, date):
    """ 
    author: claudia Acquistapace
    date; 25 July 2019 (heat wave in Cologne)
    contact: cacquist@meteo.uni-koeln.de
    goal: derive thermodinamic quantities of interest for the analysis: 
     output: dictionary containing the following variables: 
                    'mixingRatio':r, 
                    'relativeHumidity':rh, 
                    'virtualTemperature':tv,
                    'cclHeight':result_ccl['z_ccl'],
                    'cclTemperature':result_ccl['T_ccl'],
                    'lclHeight':lclArray,
                    'surfaceTemperature':TSurf, 
                    'virtualPotentialTemperature':Theta_v,
                    'time':time, 
                    'height':height,
                    'LTS':LTS,
     input:         matrices of Pressure (pa), 
                    temperature (K), 
                    absolute humidity (Kg/Kg), 
                    time, 
                    height, 
                    height of the surface
    """
    r  = np.zeros((len(time), len(height)))
    rh = np.zeros((len(time), len(height)))
    tv = np.zeros((len(time), len(height)))
        
    # calculation of mixing ratio and relative humidity
    T0 = 273.15
    for iTime in range(len(time)):
        for iHeight in range(len(height)):
            r[iTime,iHeight]  = (Q[iTime, iHeight])/(1-Q[iTime, iHeight])
            rh[iTime,iHeight] = 0.263*P[iTime, iHeight] * \
            Q[iTime, iHeight] * (np.exp( 17.67 * (T[iTime, iHeight]-T0) / (T[iTime, iHeight] - 29.65)))**(-1)

    #print('RH', rh[0,0], rh[0,-1])
    #print('pressure ' , P[0,0]*0.001, P[0,-1]*0.001)
    #print('temp', T[0,0], T[0,-1])


    # calculation of virtual temperature
    for indH in range(len(height)-1, 0, -1):
        tv[:,indH] = T[:,indH]*(1+0.608 * Q[:,indH])
            
        
    # calculation of convective condensation level (CCL) height and temperature

    from myFunctions import f_CCL_new # input variables: T, P, RH, height, time, date
    #(provide temperaturein K,  RH in % ( 70.14), P In Kpa)
    result_ccl = f_CCL_new(T, P*0.001, rh, height, time, date)
    #print(result_ccl['z_ccl'])
    #    DatasetOut = {'time':time,
    #              'height':height,
    #              'z_ccl':z_ccl,
    #              'T_ground_ccl':T_ground_CCL,
    #              'T_top_ccl':T_cclTop,
    #              'T_dew':Td}

    
    from myFunctions import lcl # T in K and P in pascal
    indSurf = len(height)-1
    PSurf    = P[:,indSurf]
    TSurf    = T[:,indSurf]         # T in K
    rhSurf   = rh[:,indSurf-1]
    lclArray = []

    for iTime in range(len(time)):
        lclArray.append(lcl(PSurf[iTime],TSurf[iTime],rhSurf[iTime]/100.))
    
    # calculation of potential and virtual potential temperature (P in pascal)
    Rd = 287.058  # gas constant for dry air [Kg-1 K-1 J]
    Cp = 1004.
    Theta = np.zeros((len(time), len(height)))
    Theta_v = np.zeros((len(time), len(height)))
    for indTime in range(len(time)):
        for indHeight in range(len(height)):
            k_val = Rd*(1-0.23*r[indTime, indHeight])/Cp
            Theta_v[indTime, indHeight] = ( (1 + 0.61 * r[indTime, indHeight]) * \
                   T[indTime, indHeight] * (100000./P[indTime, indHeight])**k_val)
            Theta[indTime, indHeight]   = T[indTime, indHeight] * (100000./P[indTime, indHeight])**k_val
    
    
    ThermodynPar={'mixingRatio':r, 
                  'relativeHumidity':rh, 
                  'virtualTemperature':tv,
                  'cclHeight':result_ccl['z_ccl'],
                  'cclTemperature':result_ccl['T_ground_ccl'],
                  'lclHeight':lclArray,
                  'surfaceTemperature':TSurf, 
                  'virtualPotentialTemperature':Theta_v,
                  'potentialTemperature':Theta,
                  'time':time, 
                  'height':height,
                  'LTS':LTS,
                  }
    
    return(ThermodynPar)

# =============================================================================
    # OBSOLETE FUNCTION >> Calculation of variance, wind speed, with direction for model outputs
    # are done in t he f_processModelOutput.py function
# def f_calcDynamics(w,u,v,thetaV,time,height,timeWindow):
#         
#     # calculating variance of vertical velocity
#     from myFunctions import f_calcWvariance
#     #print('calculating variance of vertical velocity for observations')
#     varW = f_calcWvariance(w,time,height,timeWindow)
#         
#     #print('Calculating PBL height with Richardson number method')
#     from myFunctions import f_calcPblHeightRN
#     PBLHeightArr    = f_calcPblHeightRN(thetaV,u,v,height,time)
# 
#     # calculation of wind direction and intensity for model output
#     windData_ICON   = f_calcWindSpeed_Dir(time, height, v, u)
#     #print('wind speed and direction calculated for ICON-LEM ')
# 
#     DynPar={'varianceW':varW, 
#             'PBLHeight':PBLHeightArr, 
#             'windSpeed':windData_ICON['windSpeed'], 
#             'windDirection':windData_ICON['windDirection'], 
#             }
#     return(DynPar)
#     
# =============================================================================
    
    
def f_resamplingMatrixCloudnet(time2change, height2change, matrix2change, timeRef, heightRef, matrixRef):
    

    # resampling Cloudnet observations on ICON time/height resolution
    # ---- defining ICON data as dataframe reference
    ICON_DF             = pd.DataFrame(matrixRef, index=timeRef, columns=heightRef)
    values              = np.empty((len(timeRef), len(height2change)))

    # ---- resampling PBL classification on ICON resolution
    print('resampling CLOUDNET observations on ICON time resolution')
    ZE_DF               = pd.DataFrame(matrix2change, index=time2change, columns=height2change)
    SelectedIndex_ZE    = getIndexList(ZE_DF, ICON_DF.index)
    ZE_resampled        = pd.DataFrame(values, index=timeRef, columns=height2change)
    ZE_resampled        = getResampledDataPd(ZE_resampled, ZE_DF, SelectedIndex_ZE)

    # ---- defining ICON data as dataframe reference
    ICON_DF_T           = pd.DataFrame(matrixRef.transpose(), index=heightRef, columns=timeRef)
    ZE_values           = ZE_resampled.values.transpose()
    ZEFinal_DF          = pd.DataFrame(ZE_values, index=height2change, columns=timeRef)
    Selectedcol_ZEfinal = getIndexList(ZEFinal_DF, ICON_DF_T.index)
    # define empty values for dataframe finer resolved
    values_ZEfinal      = np.empty((len(heightRef), len(timeRef)))

    # define dataframe coarse resolution
    ZE_less             = pd.DataFrame(ZE_values, index=height2change, columns=timeRef)
    # define output dataframe with resolution as icon
    ZE                  = pd.DataFrame(values_ZEfinal, index=heightRef, columns=timeRef)
    ZE                  = getResampledDataPd(ZE, ZE_less, Selectedcol_ZEfinal)  
    matrix              = np.ma.array(ZE.values, mask=np.isnan(ZE.values))
    return(matrix)








def f_calculateAllCloudQuantities(CloudInfo, \
                                  time, \
                                  height, \
                                  LWP, \
                                  LWC, \
                                  humanInfo, \
                                  Hwind, \
                                  Wwind, \
                                  yy, \
                                  dd, \
                                  mm, \
                                  QiThreshold, \
                                  QcThreshold, \
                                  iconLemData, \
                                  device, \
                                  verboseFlag, \
                                  debuggingFlag, \
                                  pathDebugFig):
    """ 
    @ author  : claudia acquistapace
    @ date    : 30 july 2019
    @ modified: 11 November 2019
    @ contact : cacquist@meteo.uni-koeln.de
    @ goal    : code develop to calculate cloud base, cloud top, cloud fraction, when they are not 
    # previously calculated, for observations and model in the same way. 
    # In addition, it identifies cloud units and counts the amount of clouds detected during
    # the day and for obs and model, it derives duration, chord length, massFlux and cloudLWP
    
    @ input   : 
    #     - CloudInfo: this input variable is - cloudnet target categorization data matrix for obs (cloudnet_res.data)
    #                                         - cloud mask for the model output
    #     - time
    #     - height
    #     - LWP
    #     - LWC: this variable corresponds to the matrix of reflectivity in linear scale when processing the observations, \
    # and corresponds to Qc matrix (filtered with respect to the Qc threshold for model data, when processing is done for ICON LEM 
    #     - cloudTimeArray: array indicating starting and ending time for the day of the PBl cloud period
    #     - Hwin
    #     - Wwind
    #     - QiThreshold: threshold value to consider for reading Qi matrix from iconlem model output
    #     - QcThreshold: threshold value to consider for reading Qc matrix from iconlem model output
    #     - iconLemData: data structure containing iconlem extracted variables
    #     - device: string specifying if observations ('obs') or model data ('iconlem') are processed.
    @ output: dictionary containing
    #            cloudMask':cloudMask, 
    #           'cloudBase':CB_array,
    #           'cloudTop':CT_array, 
    #           'liquidCloudFraction':mean_CF_liquid (mean profiles calculated over 30 minutes)
    #           'iceCloudFraction':mean_CF_ice (mean profiles calculated over 30 minutes)
    #           'totalCloudFraction':mean_CF_tot (mean profiles calculated over 30 minutes)
    #           'datetimeCloudFraction':datetime_CF  (corresponding mean time array for cloud fraction)
    #           'heightCloudFraction':height,
    #           'duration':duration of each cloud found.
    #           'cloudLWP':cloudLWP of each cloud found.
    #           'chordLength':chordLength of each cloud found. 
    #           'massFlux':massFlux of each cloud found.
    #           'Nclouds':Nclouds,
    """
    from myFunctions import f_cloudmask
    from myFunctions import f_calcCloudBaseTopPBLcloudsV2
    from myFunctions import f_calcCloudBaseTopPBLclouds
    from myFunctions import f_closest
    from myFunctions import f_calculateCloudFractionICON
    from myFunctions import f_calculateCloudProperties
    from cloudnetFunctions import f_calculateCloudFractionCloudnet
    from cloudnetFunctions import f_calculateCloudMaskCloudnet 
    from myFunctions import f_plotTest
    date = str(yy)+str(mm)+str(dd)
        
    # checking if the input is model or observations: in this case calculation 
    # of CB, CT, cloud fraction, cloud mask, reassignement of the variable LWC 
    # to ZE_lin as it really is, filtering Ze between cloud base and cloud top 
    # and calculating corresponding LWC profile for each CB/CT identified
    if device == 'obs':
        
        stringData     = 'obs'
        # calculating cloud mask for obs
        CategoryCN_res    = CloudInfo.data
        cloudMask         = f_calculateCloudMaskCloudnet(time, height, \
                                        CategoryCN_res.transpose().astype(int))
        #PLOTDONE1 = f_plotTest(cloudMask.transpose(), time, height, 'pre_call')
        # calculating 30 min mean profiles of cloud fraction for observations
        cloudFraction     = f_calculateCloudFractionCloudnet(CategoryCN_res,\
                                                yy, mm, dd, time, height)
        Ze_lin            = LWC

        # calculating cloud base , cloud top and cloud thickness for all clouds and for pbl clouds
        clouds, PBLclouds = f_calculateCloudBaseTopThickness(cloudMask, time, height, humanInfo)

        # deriving lowest cloud base and corresponding cloud top for PBL clouds
        CBarr = np.zeros(len(time))
        CBarr.fill(np.nan)
        CTarr = np.zeros(len(time))
        CTarr.fill(np.nan)
        iPBL = 0
        for itime in range(len(time)):
            if iPBL < len(PBLclouds.time.values):
                if clouds.time.values[itime] == PBLclouds.time.values[iPBL]:
                    #print(iPBL)
                    CBarray = PBLclouds.cloudBase.values[iPBL, :]
                    if CBarray.size - np.count_nonzero(np.isnan(CBarray)) != 0:
                        minCB = np.nanmin(PBLclouds.cloudBase.values[iPBL, :])
                        CBarr[itime] = minCB
                        indexLevelMin = np.nanargmin(PBLclouds.cloudBase.values[iPBL, :])
                        CTarr[itime] = PBLclouds.cloudTop[iPBL, indexLevelMin]
                    else:
                        CBarr[itime] = np.nan
                        CTarr[itime] = np.nan
                    iPBL = iPBL + 1

        # filtering Ze linear between cloud base and cloud top 
        for indT in range(len(time)):
            ZE_lin_DF          = pd.Series(Ze_lin[indT,:], index=height)
            
            if (~np.isnan(CBarr[indT])):
                mask_CB            = (ZE_lin_DF.index < CBarr[indT])
                ZE_lin_DF[mask_CB] = np.nan
            else: # case in which cloud base is nan (excluded clouds or no clouds)
                ZE_lin_DF[:] = np.nan
            if (~np.isnan(CTarr[indT])):
                mask_CT            = (ZE_lin_DF.index > CTarr[indT])
                ZE_lin_DF[mask_CT] = np.nan            
            else: # case in which cloud top is nan (excluded clouds or no clouds)
                ZE_lin_DF[:] = np.nan
                
            # copying values from the selection in the matrix
            Ze_lin[indT,:]     = ZE_lin_DF.values        
        
        # calculating LWC matrix adopting Frisch approach using Ze and LWP
        deltaZ            = 30. # range gate resolution of JOYRAD35 in meters 
        LWC               = f_calculateLWCFrisch(Ze_lin, deltaZ, time, \
                                                 height, LWP)        
        if verboseFlag == 1:
            print('cloud base, cloud top, cloud fraction and cloud mask calculated for observation')    
            print('filtering Ze linear values between cloud base and cloud top (PBL) for each time, \
                  removed values of Ze for clouds excluded from the data')
            print('calculating LWC with Frisch approach from Ze for observations')
    
    # check if input is from model
    if device == 'iconlem':
        stringData     = 'iconlem'
        Qi             = iconLemData.groups['Temp_data'].variables['Qi'][:].copy()
        Qc             = iconLemData.groups['Temp_data'].variables['Qc'][:].copy()
        # calculating cloud mask
        cloudMask      = f_cloudmask(time,height,Qc,Qi,QiThreshold,QcThreshold)

        # calculating cloud base , cloud top and cloud thickness for all clouds and for pbl clouds
        clouds, PBLclouds = f_calculateCloudBaseTopThickness(cloudMask, time, height, humanInfo)

        # deriving lowest cloud base and corresponding cloud top for PBL clouds
        CBarr = np.zeros(len(time))
        CBarr.fill(np.nan)
        CTarr = np.zeros(len(time))
        CTarr.fill(np.nan)
        iPBL = 0
        for itime in range(len(time)):
            if iPBL < len(PBLclouds.time.values):
                if clouds.time.values[itime] == PBLclouds.time.values[iPBL]:
                    #print(iPBL)
                    CBarray = PBLclouds.cloudBase.values[iPBL, :]
                    if CBarray.size - np.count_nonzero(np.isnan(CBarray)) != 0:
                        minCB = np.nanmin(PBLclouds.cloudBase.values[iPBL, :])
                        CBarr[itime] = minCB
                        indexLevelMin = np.nanargmin(PBLclouds.cloudBase.values[iPBL, :])
                        CTarr[itime] = PBLclouds.cloudTop[iPBL, indexLevelMin]
                    else:
                        CBarr[itime] = np.nan
                        CTarr[itime] = np.nan
                    iPBL = iPBL + 1

        # calculating 30 min mean profiles of cloud fraction for ICON-LEM
        cloudFraction  = f_calculateCloudFractionICON(Qi, Qc, \
                        yy, mm, dd, time, height, QiThreshold, QcThreshold)
        LWC            = Qc
        
        if verboseFlag == 1:
            print('cloud fraction and cloud mask calculated for model output')    

    # --------------------------------------------------------------------
    #  # starting common processing for obs and model data based on definitions of the variables given above
    # --------------------------------------------------------------------   

    # calculation of cloud fraction 
    mean_CF_liquid = cloudFraction['LiquidCloudFraction']
    mean_CF_ice    = cloudFraction['IceCloudFraction'] 
    mean_CF_tot    = cloudFraction['TotalCloudFraction'] 
    datetime_CF    = cloudFraction['time']
    
    # updraft speed at cloud base
    UpdraftCB = np.zeros(len(time))
    UpdraftCB.fill(np.nan)


    for indT in range(len(time)):
        if (~np.isnan(CBarr[indT])):
            indCB               = f_closest(height, CBarr[indT])
            UpdraftCB[indT]     = Wwind[indT,indCB]

    # calculation of cloud duration, chord length, mass flux and mean LWP for each cloud unit identified
    Dict_Clouds_arr  = f_calculateCloudProperties(time, \
                                                  height, \
                                                  CBarr, \
                                                  CTarr, \
                                                  UpdraftCB, \
                                                  LWP, \
                                                  Hwind, \
                                                  Wwind, \
                                                  LWC)

    duration           = []
    chordLength        = []
    massFlux           = []
    cloudLWP           = []
    cloudTimeStart     = []
    cloudTimeEnd       = []
    meanCT             = []
    meanCB             = []
    meanCloudThickness = []
    meanUpdraftCB      = []
    Nclouds            = len(Dict_Clouds_arr)
    meanheightFromCB   = np.zeros((Nclouds, len(height)))
    cloudLWC           = np.zeros((Nclouds, len(height))) 
    cloudLWC.fill(np.nan)
    meanheightFromCB.fill(np.nan)
    
    # building arrays in which each element correspond to a cloud
    for iCloud in range(Nclouds):
        duration.append(Dict_Clouds_arr[iCloud]['duration'].total_seconds())
        chordLength.append(Dict_Clouds_arr[iCloud]['chordLength'])
        massFlux.append(Dict_Clouds_arr[iCloud]['MassFlux'])    
        cloudLWP.append(Dict_Clouds_arr[iCloud]['meanLWP'])    
        cloudLWC[iCloud,:] = Dict_Clouds_arr[iCloud]['meanLWC']
        #.append(Dict_Clouds_arr[iCloud]['meanLWC'])
        cloudTimeStart.append(Dict_Clouds_arr[iCloud]['timeStart'])
        cloudTimeEnd.append(Dict_Clouds_arr[iCloud]['timeEnd'])
        #meanheightFromCB.append(Dict_Clouds_arr[iCloud]['meanheightFromCB'])
        meanheightFromCB[iCloud,:] = Dict_Clouds_arr[iCloud]['meanheightFromCB']
        meanCT.append(Dict_Clouds_arr[iCloud]['meanCT'])
        meanCB.append(Dict_Clouds_arr[iCloud]['meanCB'])
        meanCloudThickness.append(Dict_Clouds_arr[iCloud]['cloudThickness'])
        meanUpdraftCB.append(Dict_Clouds_arr[iCloud]['meanUpdraftSpeedCB'])
        
    if verboseFlag == 1:
        print('cloud properties of duration, chord length, mass flux, mean cloud LWC and LWP calculated')


    # output of a plot is keyword is selected:
    if debuggingFlag == 1:
        
        # plot of cloud mask, cloud base and cloud top normal and filtered for PBL
        fig, ax = plt.subplots(figsize=(10,4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left() 
        ax.xaxis_date()
        cax = ax.pcolormesh(time, height, cloudMask.transpose(), vmin=0, vmax=3, cmap=plt.cm.get_cmap("RdPu", 4))
        ax.set_ylim(0,12000.)                                               # limits of the y-axes
        #ax.set_xlim(0,24)                                 # limits of the x-axes
        ax.set_title("cloud mask", fontsize=14)
        ax.set_xlabel("time ", fontsize=12)
        ax.set_ylabel("height [m]", fontsize=12)
        plt.plot(time, CBarr, color='black', label='CB PBL')
        plt.plot(time, CTarr, color='black', linestyle=':', label='CT PBL')
        plt.legend()
        cbar = fig.colorbar(cax, ticks=[0, 1, 2, 3], orientation='vertical')
        cbar.ticks=([0,1,2,3])
        cbar.ax.set_yticklabels(['no cloud','liquid','ice', 'mixed phase'])
        cbar.set_label(label="cloud type",size=12)
        cbar.ax.tick_params(labelsize=12)
        cbar.aspect=80
        plt.savefig(pathDebugFig+'cloudMask_'+stringData+'_'+date+'.png', format='png')

    # define output dictionary of data: lenght of each element of the dictionary 
    # is given by the number of clouds found in the day
    dictOut = {}
    dictOut = {'cloudMask':cloudMask, 
               'timeSerie':time,
               'liquidCloudFraction':mean_CF_liquid,
               'iceCloudFraction':mean_CF_ice, 
               'totalCloudFraction':mean_CF_tot, 
               'datetimeCloudFraction':datetime_CF, 
               'heightCloudFraction':height,
               'duration':duration,
               'cloudLWP':cloudLWP,
               'cloudLWC':cloudLWC,
               'meanheightFromCB':meanheightFromCB,
               'cloudMeanCB':meanCB,
               'cloudMeanCT':meanCT,
               'cloudTimeEnd':cloudTimeEnd,
               'chordLength':chordLength, 
               'massFlux':massFlux, 
               'timeCloudStart':cloudTimeStart, 
               'timeCloudEnd':cloudTimeEnd,
               'cloudUpdraftCB':meanUpdraftCB,
               'Nclouds':Nclouds,
               'LWPall':LWP,
               'UpdraftCB_PBL':UpdraftCB,
            }

    return(dictOut, clouds, PBLclouds)




#---------------------------------------------------------------------------------
# date :  28.01.2019
# author: Claudia Acquistapace
# goal: function that identifies cloud continuous entities and derives ( based on definitions in Laureau et al, 2018, JAS)
#  - cloud duration
#    chord lenght
#    Mass flux
#    mean cloud base height
#    vertical wind one range gate below cloud base
#    mean LWP per cloud entity
#    
# output: array of dictionaries: every dictionary correspond to an identified cloud, and contains the clodu properties
# -DictCloudPropArr
#--------------------------------------------------------------------------------
def f_calculateCloudProperties(datetime_ICON, \
                               height_ICON, \
                               CBarr, \
                               CTarr, \
                               UpdraftCB, \
                               LWP_ICON, \
                               Hwind_ICON, \
                               w_ICON, \
                               LWC):
    """ 
    @ author: cacquist
    @ date  : November 2019
    @ goal  : this function was created to process in the same way model and obs data
    It gets the inputs and calculates various cloud useful quantities, listed below.
    locaql variables are called *_ICON because originally the code was developed only 
    for icon variables, but it applies to obs and model data
    @ INPUT : 
        datetime_ICON: time array
        height_ICON  : height array
        PBLclouds    : xarray dataset containing PBl cloud bases, tops, and thicknesses
        UpdraftCB    : updraft velocity at cloud base time serie
        LWP_ICON     : LWP time serie
        Hwind_ICON   : Horizontal wind matrix (time, height)
        w_ICON       : vertical wind matrix (time, height)
        LWC          : liquid water content / Qc matrix (time, height)
    @ OUTPUT: 
        one dictionary for each cloud unit. Each dictionary contains:
        dictProp = {'timeStart':timeStart, 
                    'indStart':indStart, 
                    'timeEnd':timeEnd, 
                    'indEnd':indEnd, 
                    'meanLWP':meanLWP, 
                    'meanLWC':meanLWC, 
                    'meanheightFromCB':meanheightFromCB,
                    'meanCT':meanCT,
                    'stdLWP':stdLWP,
                    'meanCB':meanCB, 
                    'stdCB':stdCB,
                    'WwindCB':WwindCloudBase, 
                    'MassFlux':MassFlux, 
                    'duration':duration, 
                    'chordLength':chordLength }
    """
    LWC[LWC == 0]   = np.nan   # setting to nans null values for better averaging
    cloudStart      = 0
    Dict_Clouds_arr = []


    
    # assigning starting and ending time of cloudy intervals considered continuous clouds
    for itime in range(len(datetime_ICON)-1):
        if (np.isnan(CBarr[itime]) == False) * (cloudStart == 0): #cb found, cloud not started yet
            cloudStart = 1
            timeStart  = datetime_ICON[itime]
            indStart   = itime
        # if cb found and cloudstart =1 loop does not do anything: this corresponds to the cloudy part
        if ((np.isnan(CBarr[itime]) == True) * (cloudStart == 1)): # se Cb not found, \
            #and cloudstart =1 (comes from a cloud), then it saves the previous time step as timeEnd \
            # and puts to zero the cloud flag again, ready for a new cloud. Saves time indeces and values in the dictionary
            # and sets timestart and end to nans
            #print('sono qua')
            timeEnd    = datetime_ICON[itime-1]
            indEnd     = itime-1
            cloudStart = 0
            dict_cloud = {'timeStart':timeStart, 'indStart':indStart, 'timeEnd':timeEnd, 'indEnd':indEnd}
            Dict_Clouds_arr.append(dict_cloud)
            timeStart  = np.nan
            timeEnd    = np.nan
            
    # filtering LWC profiles below CB and above cloud top to nans
    for itime in range(len(datetime_ICON)):
        #if ((~np.isnan(minCB)) or (~np.isnan(maxCT))):
        if ((~np.isnan(CBarr[itime])) or (~np.isnan(CTarr[itime]))):
            #fig, ax          = plt.subplots(figsize=(4,4))
            #plt.plot(height_ICON, LWC[itime,:], label='before', color='red')
            LWC_prof_DF              = pd.Series(LWC[itime,:], index=height_ICON)
            mask_CB                  = (LWC_prof_DF.index < CBarr[itime])
            LWC_prof_DF.loc[mask_CB] = np.nan
            mask_CT                  = (LWC_prof_DF.index > CTarr[itime])
            LWC_prof_DF.loc[mask_CT] = np.nan
            LWC[itime,:]             = LWC_prof_DF.values
            #plt.plot(height_ICON, LWC[itime,:], label='after', color='blue', linestyle=':')
            #plt.axvline(x=CB_array_ICON[itime])
            #plt.axvline(x=CT_array_ICON[itime])
            #plt.xlim(0., 6000.)
            #plt.savefig('/work/cacquist/HDCP2_S2/statistics/debug/'+str(itime))
        else:
            LWC[itime,:]             = np.repeat(np.nan, len(height_ICON))
    

    # calculating cloud duration, velocity below cloud base, mean cloud base height,\
    # cloud chord lenght and corresponding mass flux as defined in (Lareau et al, 2018., JAS)
    DictCloudPropArr = []
    for iCloud in range(len(Dict_Clouds_arr)):
        
        timeStart  = Dict_Clouds_arr[iCloud]['timeStart']
        timeEnd    = Dict_Clouds_arr[iCloud]['timeEnd']
        duration   = timeEnd - timeStart
        
        # calculating corresponding mean LWP, std LWP, mean CB height, 
        iTimeStart         = Dict_Clouds_arr[iCloud]['indStart']
        iTimeEnd           = Dict_Clouds_arr[iCloud]['indEnd']
        meanLWP            = np.nanmedian(LWP_ICON[iTimeStart:iTimeEnd])
        stdLWP             = np.nanstd(LWP_ICON[iTimeStart:iTimeEnd])
        meanLWC            = np.nanmean(LWC[iTimeStart:iTimeEnd,:], axis=0)
        meanUpdraftCBspeed = np.nanmedian(UpdraftCB[iTimeStart:iTimeEnd])
        stdUpdraftCBspeed  = np.nanstd(UpdraftCB[iTimeStart:iTimeEnd])
       
        
        # finding max and min height where LWC array is non null: check the number of non nan elements
        if (np.count_nonzero(~np.isnan(meanLWC)) != 0):
            meanCB           = np.nanmean(CBarr[iTimeStart:iTimeEnd])
            meanCT           = np.nanmean(CTarr[iTimeStart:iTimeEnd])
            # calculating cloud properties 
            HwindCloudBase   = np.nanmean(Hwind_ICON[iTimeStart:iTimeEnd, f_closest(height_ICON, meanCB)])
            WwindCloudBase   = np.nanmean(w_ICON[iTimeStart:iTimeEnd, f_closest(height_ICON, meanCB)+1])
            chordLength      = HwindCloudBase * duration.total_seconds()
            meanheightFromCB = (height_ICON - np.repeat(meanCB, len(height_ICON)))/ \
            (np.repeat(meanCT, len(height_ICON))- np.repeat(meanCB, len(height_ICON)))
            MassFlux         = WwindCloudBase * chordLength 
            cloudThickness   = meanCT - meanCB
        else:
            if (len(CBarr[iTimeStart:iTimeEnd]) != 0):
                meanCB           = np.nanmean(CBarr[iTimeStart:iTimeEnd])
                meanCT           = np.nanmean(CTarr[iTimeStart:iTimeEnd])
                HwindCloudBase   = np.nanmean(Hwind_ICON[iTimeStart:iTimeEnd, f_closest(height_ICON, meanCB)])
                WwindCloudBase   = np.nanmean(w_ICON[iTimeStart:iTimeEnd, f_closest(height_ICON, meanCB)+1])
                chordLength      = HwindCloudBase * duration.total_seconds()
                meanheightFromCB = (height_ICON - np.repeat(meanCB, len(height_ICON)))/ \
                (np.repeat(meanCT, len(height_ICON))- np.repeat(meanCB, len(height_ICON)))
                MassFlux         =  WwindCloudBase * chordLength
                cloudThickness   = meanCT - meanCB

            else:
                meanCB           = np.nan
                meanCT           = np.nan
                HwindCloudBase   = np.nan
                WwindCloudBase   = np.nan
                chordLength      = np.nan
                meanheightFromCB = np.repeat(np.nan, len(height_ICON))
                MassFlux         = np.nan
                cloudThickness   = np.nan

        # storing data in dictionary array, each array is a cloudy unit
        dictProp = {'timeStart':timeStart, 
                    'indStart':indStart, 
                    'timeEnd':timeEnd, 
                    'indEnd':indEnd, 
                    'meanLWP':meanLWP, 
                    'meanLWC':meanLWC, 
                    'meanheightFromCB':meanheightFromCB,
                    'meanCT':meanCT,
                    'stdLWP':stdLWP,
                    'meanCB':meanCB, 
                    'WwindCB':WwindCloudBase,
                    'MassFlux':MassFlux, 
                    'duration':duration, 
                    'chordLength':chordLength, 
                    'cloudThickness':cloudThickness,
                    'meanUpdraftSpeedCB':meanUpdraftCBspeed,
                    'stdUpdraftSpeedCB':stdUpdraftCBspeed,
                    }
        DictCloudPropArr.append(dictProp)
        
    return(DictCloudPropArr)
    





def f_calcMeanStdVarProfiles(field, time, height, date, yy, mm, dd, NprofilesOut, timeIncrement):
    """date  : 02 aug 2019
    # author : claudia Acquistapace
    # contact: cacquist@meteo.uni-koeln.de
    # goal   : function to calculate mean profiles of a field over a given time period. 
    # input  :
    #    - field: matrix to be averaged (it has to have dimensions (dim(time), dim(height)))
    #    - time
    #    - height
    #    - date
    #    - yy (integer for year)
    #    - mm (integer for month)
    #    - dd (integer for day)
    #    - NprofilesOut: number of averaged profiles to get as outputs: it depends on the interval of time chosed 
    #      for the average Es( 48, timeIncrement= 30 min)
    #    - timeIncrement: interval of time chosed for the average (in minutes) 
    
    # output:
    # dictionary of:
    orderedDict['meanProfiles']
    orderedDict['stdProfiles']
    orderedDict['meanTime']
    orderedDict['height']
    orderedDict['date']
    """
    # defining dataframes for calculating mean    
    field_DF        = pd.DataFrame(field, index=time, columns=height)
    #print(field_DF)
    Profiles_var_DF = pd.DataFrame(np.zeros((len(height),NprofilesOut)), \
                                   columns=np.arange(0,NprofilesOut), index=height)
    Std_var_DF      = pd.DataFrame(np.zeros((len(height),NprofilesOut)), \
                                   columns=np.arange(0,NprofilesOut), index=height)



    import collections
    deltaT       = datetime.timedelta(minutes=timeIncrement)
    indInt       = 0
    datetime_out = []

    for itime in range(0,NprofilesOut):
        
        if indInt == 0:
            HourInf = datetime.datetime(int(yy), int(mm), int(dd), 0, 0, 0) 
        else:
            HourInf = HourInf + deltaT
        HourSup                       = HourInf + deltaT
        datetime_out.append(HourInf)
        indInt                        = indInt + 1
        field_sliced_t                = field_DF.loc[(field_DF.index < HourSup) * (field_DF.index >= HourInf),:]
        field_mean                    = field_sliced_t.mean(axis=0, skipna=True)
        field_std                     = field_sliced_t.std(axis=0, skipna=True)
        #print(field_mean)
        Profiles_var_DF.loc[:,indInt] = field_mean
        Std_var_DF.loc[:,indInt]      = field_std
        
    
    orderedDict                       = collections.OrderedDict()
    orderedDict['meanProfiles']       = Profiles_var_DF.values
    orderedDict['stdProfiles']        = Std_var_DF.values
    orderedDict['meanTime']           = datetime_out
    orderedDict['height']             = height
    orderedDict['date']               = date

    return(orderedDict)







def f_plotCloudFraction(datetime_CF, height, pathFig, CFmean_mod, CFmean_obs, \
                        CFmeanLiquid_mod, CFmeanLiquid_obs, CFmeanIce_mod, CFmeanIce_obs):
    import matplotlib as mpl

    import matplotlib.dates as mdates
    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmean_mod.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction icon-lem - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_tot_wholeDataset_mod.png', format='png')  
    
    
    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmean_obs.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction obs - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_tot_wholeDataset_obs.png', format='png')  
    
    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmeanLiquid_mod.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction icon-lem - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_liq_wholeDataset_mod.png', format='png')  
    
    
    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmeanLiquid_obs.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction obs - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_liq_wholeDataset_obs.png', format='png')  
    

    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmeanIce_mod.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction icon-lem - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_ice_wholeDataset_mod.png', format='png')  
    
    
    fig, ax          = plt.subplots(figsize=(12,6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left() 
    ax.xaxis_date()
    label_size = 16
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    cax1             = ax.pcolormesh(datetime_CF, height, CFmeanIce_obs.transpose(), vmin=0., vmax=0.4, cmap='BuPu')
    ax.set_ylim(400.,4000.)                                               # limits of the y-axe
    ax.set_xlim()                                                        # limits of the x-axes
    ax.set_title("cloud fraction obs - JOYCE", fontsize=16)
    ax.set_xlabel("time [hh:mm]", fontsize=16)
    ax.set_ylabel("height [m]", fontsize=16)
            #plt.plot(time_ICON, CT_array, color='black', label='cloud top')
            #plt.plot(time_ICON, CB_array, color='black',label='cloud base')
    plt.legend(loc='upper left')
    cbar = fig.colorbar(cax1, orientation='vertical')
            #cbar.ticks=([0,1,2,3])
            #cbar.ax.set_yticklabels(['no cloud','liquid','ice','mixed phase'])
    cbar.set_label(label="cloud fraction ",size=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.aspect=80
    fig.tight_layout()
    plt.savefig(pathFig+'cloudFraction_ice_wholeDataset_obs.png', format='png')  
    
def f_convertPressureToHeight(P,T):
         # pressure in Kpascal
         # T in Kelvin
    import math
    P0 = 101.325 # Kpa
    pippo = P/P0
    listLog = []
    for ind in range(len(P)):
        listLog.append(np.log10(pippo[ind]))
    g = 9.807 # m/s^2
    M = 0.02896  # Kg/mol
    R = 8.3143 # Nm/mol K
    h = - (R*T)*listLog/(M*g)

    return(h)
    




def f_plotVarianceWSingleDays(date,varWmean_obs,varWmean_mod, varWstd_obs, \
                              varWstd_mod,indHourPlotStart, height, pathFig):
    
    
    Nrows = 2
    Ncols = 5
    Nplots = Nrows*Ncols+1
    
    # ---- plotting hourly profiles of variance of vertical velocity during the day
    fig, ax       = plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(14,10))
    #matplotlib.rcParams['savefig.dpi'] = 300
    plt.gcf().subplots_adjust(bottom=0.15)
    fig.tight_layout()
    ymax          = 3000.
    ymin          = 107.
    xmax          = 1.5
    fontSizeTitle = 16
    fontSizeX     = 15
    fontSizeY     = 15
    #timeTitles = [']
    for indPlot in range(1, Nplots):
        ax        = plt.subplot(2,5,indPlot)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)  
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left() 
        #ax.text(1.8, ymax-200., 'a)', fontsize=15)
        matplotlib.rc('xtick', labelsize=15)                        # sets dimension of ticks in the plots
        matplotlib.rc('ytick', labelsize=15)                        # sets dimension of ticks in the plots
        plt.plot(varWmean_obs[:,indHourPlotStart], height, label='obs',  color='black')
    #plt.errorbar(varWmean_obs[:,8], height, xerr=varWstd_obs[:,8], color='black')
        plt.plot(varWmean_mod[:,indHourPlotStart], height, label='ICON',  color='red')
    #plt.errorbar(varWmean_mod[:,8], height, xerr=varWstd_mod[:,8], color='red')
        y1        = varWmean_obs[:,indHourPlotStart]-varWstd_obs[:,indHourPlotStart]
        y2        = varWmean_obs[:,indHourPlotStart]+varWstd_obs[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2> y1, facecolor='black', alpha=0.2)
        y1        = varWmean_mod[:,indHourPlotStart]-varWstd_mod[:,indHourPlotStart]
        y2        = varWmean_mod[:,indHourPlotStart]+varWstd_mod[:,indHourPlotStart]
        plt.fill_betweenx(height, y1, y2, where=y2> y1, facecolor='red', alpha=0.2)
        plt.legend(loc='upper right', fontsize=14)
        plt.ylim(ymin,ymax)
        plt.xlim(0.,xmax)
        #plt.title('8:00 UTC', fontsize=fontSizeTitle)
        plt.xlabel('var(w) [m/s]', fontsize=fontSizeX)
        plt.ylabel('height [m]', fontsize=fontSizeY)
        plt.tight_layout()
        
        indHourPlotStart = indHourPlotStart+1
    
    plt.savefig(pathFig+'varW_Profiles_diurnal_cycle_'+date+'.png', format='png')







# finding how many radiosondes are launched during the selected day:
def f_reshapeRadiosondes(new_dict):
    #from myFunctions import f_convertPressureToHeight
         
         
    P_radios_obs       = []
    T_radios_obs       = []
    theta_v_radios_obs = [] 
    RH_radios_obs      = []
    height_radios_obs  = []
    time_radios_obs    = []
    lengthRadiosonde   = []
    lcl_radios_obs     = []
    ccl_radios_obs     = []
    lts_radios_obs     = []
    pblHeight_radios_obs     = []    
    
    Nsoundings         = len(new_dict[0])
         
    for indSoundings in range(Nsoundings):
        P_radios_obs.append(new_dict[0][indSoundings]['P']/10.)
        T_radios_obs.append(new_dict[0][indSoundings]['T'])
        theta_v_radios_obs.append(new_dict[0][indSoundings]['theta_v'])
        RH_radios_obs.append(new_dict[0][indSoundings]['RH'])
        #height_radios_obs.append(f_convertPressureToHeight(new_dict[0][indSoundings]['P']/10.,new_dict[0][indSoundings]['T']))
        height_radios_obs.append(new_dict[0][indSoundings]['height'])
        time_radios_obs.append(new_dict[0][indSoundings]['time'])
        lengthRadiosonde.append(len(new_dict[0][indSoundings]['P']))
        lcl_radios_obs.append(new_dict[0][indSoundings]['z_lcl'])
        ccl_radios_obs.append(new_dict[0][indSoundings]['z_ccl'])
        lts_radios_obs.append(new_dict[0][indSoundings]['LTS'])
        pblHeight_radios_obs.append(new_dict[0][indSoundings]['PBLheight'])
        
        
    # find minimum length of radiosondes
    P_resized          = []
    T_resized          = []
    theta_v_resized    = []
    T_resized          = []
    RH_resized         = []
    height_resized     = []
    index_min          = np.argmin(lengthRadiosonde)
    lenghtMin          = lengthRadiosonde[index_min]

    for indSoundings in range(Nsoundings):   
        if indSoundings == index_min:
            P_resized.append(P_radios_obs[indSoundings])
            T_resized.append(T_radios_obs[indSoundings])
            theta_v_resized.append(theta_v_radios_obs[indSoundings])
            RH_resized.append(RH_radios_obs[indSoundings])
            height_resized.append(height_radios_obs[indSoundings])
        else:
            P_resized.append(P_radios_obs[indSoundings][0:lenghtMin])
            T_resized.append(T_radios_obs[indSoundings][0:lenghtMin])
            RH_resized.append(RH_radios_obs[indSoundings][0:lenghtMin])
            theta_v_resized.append(theta_v_radios_obs[indSoundings][0:lenghtMin])        
            height_resized.append(height_radios_obs[indSoundings][0:lenghtMin])                


    # building matrices
    P_radiosonde_obs       = np.reshape(P_resized, (lenghtMin, Nsoundings))
    T_radiosonde_obs       = np.reshape(T_resized, (lenghtMin, Nsoundings))
    RH_radiosonde_obs      = np.reshape(RH_resized, (lenghtMin, Nsoundings))
    theta_v_radiosonde_obs = np.reshape(theta_v_resized, (lenghtMin, Nsoundings))
    height_radiosonde_obs  = np.reshape(height_resized, (lenghtMin, Nsoundings))

         
    dict_radios = {'P':P_radiosonde_obs,
                   'T':T_radiosonde_obs,
                   'RH':RH_radiosonde_obs,
                   'theta_v':theta_v_radiosonde_obs,
                   'height':height_radiosonde_obs, 
                   'time':time_radios_obs, 
                   'lcl':lcl_radios_obs,
                   'ccl':ccl_radios_obs, 
                   'lts':lts_radios_obs, 
                   'pblHeight':pblHeight_radios_obs}
    return(dict_radios)
    
    

# function to calculate mean theta_v from the model around the hours of the radiosondes and averaging 
# =============================================================================
"""the function goes through the number of days of the dataset. For each day, it counts the number of 
radiosonde launched and reads each of them. For each hour corresponding to a radiosonde of the day, 
it calculates the mean quantities of the model around that hour (+-1 hour)
it then returns for every hour, a dictionary in which radiosonde data and 
the corresponding mean model quantities are stored togehter. 
Every dictionary associated to a given hour is appended to a list which is 
piling all hours together, independently of the day. 
each element of the list is a dictionary of an hour """
def f_calculateMeanThetaVModelProfiles(time_radiosondes, \
                                       theta_v_radiosondes,\
                                       T_radiosondes, \
                                       rh_radiosObs, \
                                       height_radiosondes, \
                                       lcl_radiosondes, \
                                       ccl_radiosondes, \
                                       lts_radiosondes, \
                                       pblHeight_radiosondes, \
                                       theta_v_mod, \
                                       T_mod, \
                                       rh_mod, \
                                       time_mod, \
                                       height_mod, \
                                       lcl_mod, \
                                       ccl_mod,\
                                       lts_mod, \
                                       pblHeight_mod):

    theta_v_dict_obs_mod_arr = []
    # time_radiosondes is a list: every element of the list corresponds to a day, 
    # and for every dat there is a different number of radiosonde launched.
    Ndays = len(time_radiosondes)
    
    # loop on the number of days: for each day there is a different number of 
    # radiosondes launched (NdarDay)
    for indDay in range(Ndays):
        
        # finding the hour of each radiosonde launched
        NradDay = len(time_radiosondes[indDay])   # number of radiosondes launched that day
        
        # loop on the number of radiosonde of the day
        for indRadDay in range(NradDay):
            #print(indRadDay)
            
            # reading the data of the selected radiosonde of the day
            radioSondeSelected = time_radiosondes[indDay][indRadDay]
            lcl_rad            = lcl_radiosondes[indDay][indRadDay]
            ccl_rad            = ccl_radiosondes[indDay][indRadDay]
            lts_rad            = lts_radiosondes[indDay][indRadDay]
            pblHeight_rad      = pblHeight_radiosondes[indDay][indRadDay]
            theta_v_rad        = theta_v_radiosondes[indDay][:,indRadDay]
            T_rad_day          = T_radiosondes[indDay][:,indRadDay]
            rh_rad_day         = rh_radiosObs[indDay][:,indRadDay]
            height_rad         = height_radiosondes[indDay][:,indRadDay]
            theta_v_day        = theta_v_mod[indDay][:,:]
            T_mod_day          = T_mod[indDay][:,:]
            rh_mod_day         = rh_mod[indDay][:,:]
            time_day           = time_mod[indDay][:]
            lcl_mod_day        = lcl_mod[indDay][:]
            ccl_mod_day        = ccl_mod[indDay][:]
            lts_mod_day        = lts_mod[indDay][:]
            pblHeight_mod_day  = pblHeight_mod[indDay][:]
            #lts_mod_day = lts_mod[:]

            # reading exact hour of the radiosounding
            hh = radioSondeSelected.hour
            dd = radioSondeSelected.day
            yy = radioSondeSelected.year
            MM = radioSondeSelected.month
            
            # defining the time interval around the hour to consider
            # for calculating the average of the model
            if hh == 23:
                hourInf = datetime.datetime(yy,MM,dd,hh-2)
                hourSup = datetime.datetime(yy,MM,dd,hh)
            if hh == 0:
                hourInf = datetime.datetime(yy,MM,dd,hh)
                hourSup = datetime.datetime(yy,MM,dd,hh+2)
            if (hh!= 23) * (hh!= 0):            
                hourInf = datetime.datetime(yy,MM,dd,hh-1)
                hourSup = datetime.datetime(yy,MM,dd,hh+1)
            #print(hh)
            
            # calculate now mean profile of the theta_v from the model output 
            # corresponding to the selected time interval
            
            # defining pandas dataframes 
            thetaV_DF            = pd.DataFrame(theta_v_day, index=time_day, columns=height_mod)
            T_mod_DF             = pd.DataFrame(T_mod_day, index=time_day, columns=height_mod)
            RH_mod_DF            = pd.DataFrame(rh_mod_day, index=time_day, columns=height_mod)
            lcl_mod_DF           = pd.DataFrame(lcl_mod_day, index=time_day)
            ccl_mod_DF           = pd.DataFrame(ccl_mod_day, index=time_day)
            lts_mod_DF           = pd.DataFrame(lts_mod_day, index=time_day)
            pblHeight_mod_DF     = pd.DataFrame(pblHeight_mod_day, index=time_day)
            
            # selecting the theta_v profiles in the time interval corresponding to the hour
            field_sliced_t       = thetaV_DF.loc[(thetaV_DF.index < hourSup) * (thetaV_DF.index >= hourInf),:]
            theta_v_mod_mean     = field_sliced_t.mean(axis=0, skipna=True)
            theta_v_mod_std      = field_sliced_t.std(axis=0, skipna=True)
            
            field_sliced_T_mod   = T_mod_DF.loc[(T_mod_DF.index < hourSup) * (T_mod_DF.index >= hourInf),:]
            T_mod_mean           = field_sliced_T_mod.mean(axis=0, skipna=True)
            T_mod_std            = field_sliced_T_mod.std(axis=0, skipna=True)
            
            field_sliced_RH_mod  = RH_mod_DF.loc[(RH_mod_DF.index < hourSup) * (RH_mod_DF.index >= hourInf),:]
            RH_mod_mean          = field_sliced_RH_mod.mean(axis=0, skipna=True)
            RH_mod_std           = field_sliced_RH_mod.std(axis=0, skipna=True)
            
            lcl_slice            = lcl_mod_DF.loc[(lcl_mod_DF.index < hourSup) * (lcl_mod_DF.index >= hourInf)]
            lcl_mod_mean         = lcl_slice.mean(skipna=True)
            lcl_mod_std          = lcl_slice.std(skipna=True)

            ccl_slice            = ccl_mod_DF.loc[(ccl_mod_DF.index < hourSup) * (ccl_mod_DF.index >= hourInf)]
            ccl_mod_mean         = ccl_slice.mean(skipna=True)
            ccl_mod_std          = ccl_slice.std(skipna=True)

            lts_slice            = lts_mod_DF.loc[(lts_mod_DF.index < hourSup) * (lts_mod_DF.index >= hourInf)]
            lts_mod_mean         = lts_slice.mean(skipna=True)
            lts_mod_std          = lts_slice.std(skipna=True)

            pblHeight_slice      = pblHeight_mod_DF.loc[(pblHeight_mod_DF.index < hourSup) * (pblHeight_mod_DF.index >= hourInf)]
            pblHeight_mod_mean   = pblHeight_slice.mean(skipna=True)
            pblHeight_mod_std    = pblHeight_slice.std(skipna=True)
            
            dict_theta = {'theta_v_radios':theta_v_rad, 
                          'T_radios':T_rad_day,
                          'RH_radios':rh_rad_day,
                          'height_rad':height_rad,
                          'theta_v_mod_mean':theta_v_mod_mean, 
                          'theta_v_mod_std':theta_v_mod_std,
                          'T_mod_mean':T_mod_mean, 
                          'T_mod_std':T_mod_std,
                          'RH_mod_mean':RH_mod_mean, 
                          'rh_mod_std':RH_mod_std, 
                          'date':radioSondeSelected,
                          'lcl_mod_mean':lcl_mod_mean,
                          'lcl_mod_std':lcl_mod_std,
                          'lcl_rad': lcl_rad,
                          'ccl_mod_mean': ccl_mod_mean,
                          'ccl_mod_std' : ccl_mod_std,
                          'ccl_rad':ccl_rad,
                          'lts_mod_mean':lts_mod_mean,
                          'lts_mod_std':lts_mod_std,
                          'lts_rad':lts_rad,
                          'pblHeight_mod_mean':pblHeight_mod_mean, 
                          'pblHeight_mod_std':pblHeight_mod_std,
                          'pblHeight_rad':pblHeight_rad,
                          'hour':hh}
            theta_v_dict_obs_mod_arr.append(dict_theta)
    return(theta_v_dict_obs_mod_arr)
    

def f_calculateMeanProfilesPlotThetaVRadiosondes(theta_v_dict_obs_mod_arr, height_mod):
    """function to derive mean profiles for each hour at which radiosondes were launched
    the first step is to sort the elements of the list of dictionaries 
    (each dictionary correspond to one radiosonde launched at a given hour) based 
    on the hour, so that we collect consequently same hours belonging to different days
    """
    import operator
    from collections import Counter
    
    # sorting dictionary based on hours: 
    theta_v_dict_obs_mod_arr.sort(key=operator.itemgetter('hour')) #list of dictionaries ordered per hour
    
    # counting how many profiles for each hour are present (on different days)
    k            = [i['hour'] for i in theta_v_dict_obs_mod_arr]
    m            = Counter(k)
    k0           = 0
    listHourDict = []
    hourPrec     = 0
    # loop on hours found
    for ind in range(len(k)):
        hourSel = k[ind]
        # if not the first hour of the loop 
        if hourSel != hourPrec:
            hourPrec          = hourSel
            #print('hoursel', hourSel)
            
            # counting how many profiles for each hour are present
            Nprofiles         = m[hourSel]
            #print('Nprofiles', Nprofiles)
            #print('k0', k0)
            
            lenghtProf_radios = []
            lenghtProf_mod    = []
            
            # loop on the number of "same" hours present to find out the 
            # lengths of the radiosonde profiles and of the modelled ones. 
            # Lengths are stored in length_mod and length_obs arrays
            for iloop in range(k0, k0+Nprofiles):
                #print('iloop', iloop)
                lenghtProf_radios.append(len(theta_v_dict_obs_mod_arr[iloop]['theta_v_radios']))
                lenghtProf_mod.append(len(theta_v_dict_obs_mod_arr[iloop]['theta_v_mod_mean']))
                
            #print('lenght of the profiles', lenghtProf_radios)
            #print('max lenght', np.max(lenghtProf_radios))
            
            # defining matrices where to store all profiles collected for the hour on which the loop is running
            # and the matrices where to calculate the mean profiles for the model data corresponding to the 
            # given hour. All matrices are filled with nans
            MatrixProfiles_radios        = np.zeros((np.max(lenghtProf_radios), Nprofiles))
            MatrixProfiles_mod           = np.zeros((np.max(lenghtProf_mod), Nprofiles))
            MatrixProfiles_T_radios      = np.zeros((np.max(lenghtProf_radios), Nprofiles))
            MatrixProfiles_T_mod         = np.zeros((np.max(lenghtProf_mod), Nprofiles))
            MatrixProfiles_RH_radios     = np.zeros((np.max(lenghtProf_radios), Nprofiles))
            MatrixProfiles_RH_mod        = np.zeros((np.max(lenghtProf_mod), Nprofiles))
            MatrixHeight_radios          = np.zeros((np.max(lenghtProf_radios), Nprofiles))
            arr_lcl_hour_radios          = np.zeros(Nprofiles)
            arr_ccl_hour_radios          = np.zeros(Nprofiles)
            arr_lts_hour_radios          = np.zeros(Nprofiles)
            arr_pblHeight_hour_radios    = np.zeros(Nprofiles)
            arr_lcl_hour_mod             = np.zeros(Nprofiles)
            arr_ccl_hour_mod             = np.zeros(Nprofiles)
            arr_lts_hour_mod             = np.zeros(Nprofiles)
            arr_pblHeight_mod            = np.zeros(Nprofiles)
            MatrixProfiles_radios.fill(np.nan) 
            MatrixProfiles_mod.fill(np.nan)
            MatrixProfiles_RH_radios.fill(np.nan) 
            MatrixProfiles_RH_mod.fill(np.nan)
            MatrixProfiles_T_radios.fill(np.nan) 
            MatrixProfiles_T_mod.fill(np.nan)            
            MatrixHeight_radios.fill(np.nan)
            arr_lcl_hour_radios.fill(np.nan)
            arr_ccl_hour_radios.fill(np.nan)
            arr_lts_hour_radios.fill(np.nan)
            arr_pblHeight_hour_radios.fill(np.nan)
            arr_lcl_hour_mod.fill(np.nan)
            arr_ccl_hour_mod.fill(np.nan)
            arr_lts_hour_mod.fill(np.nan)
            arr_pblHeight_mod.fill(np.nan)
            
            # loop on the number of profiles for the given hour to fill the data in the matrices
            for iloop in range(Nprofiles):
                #print('calculating mean for hour', hourSel)
                #print(iloop+k0)
                
                # filling mean profiles from the model (mean calculated around that hour)
                MatrixProfiles_mod[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['theta_v_mod_mean']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['theta_v_mod_mean']
                MatrixProfiles_T_mod[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['T_mod_mean']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['T_mod_mean']
                MatrixProfiles_RH_mod[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['RH_mod_mean']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['RH_mod_mean']                
                
                
                # filling profiles from radiosondes
                MatrixProfiles_radios[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['theta_v_radios']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['theta_v_radios']
                MatrixProfiles_T_radios[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['T_radios']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['T_radios']
                MatrixProfiles_RH_radios[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['RH_radios']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['RH_radios']                
                
                MatrixHeight_radios[0:len(theta_v_dict_obs_mod_arr[iloop+k0]['height_rad']),iloop]\
                = theta_v_dict_obs_mod_arr[iloop+k0]['height_rad'] 
                
                # filling values of lcl, lts and pbl from all radiosondes of the same hour from different days
                arr_lcl_hour_radios[iloop]       = theta_v_dict_obs_mod_arr[iloop+k0]['lcl_rad']
                arr_ccl_hour_radios[iloop]       = theta_v_dict_obs_mod_arr[iloop+k0]['ccl_rad']
                arr_lts_hour_radios[iloop]       = theta_v_dict_obs_mod_arr[iloop+k0]['lts_rad']
                arr_pblHeight_hour_radios[iloop] = theta_v_dict_obs_mod_arr[iloop+k0]['pblHeight_rad']
                
                # filling valued of lcl,lts and pbl from model mean already calculated
                arr_lcl_hour_mod[iloop]          = theta_v_dict_obs_mod_arr[iloop+k0]['lcl_mod_mean']
                arr_ccl_hour_mod[iloop]          = theta_v_dict_obs_mod_arr[iloop+k0]['ccl_mod_mean']
                arr_lts_hour_mod[iloop]          = theta_v_dict_obs_mod_arr[iloop+k0]['lts_mod_mean']
                arr_pblHeight_mod[iloop]         = theta_v_dict_obs_mod_arr[iloop+k0]['pblHeight_mod_mean']
                #print('lts obs', arr_lts_hour_radios)
                #print('lts mod', arr_lts_hour_mod)
            
            # incrementing Ko of the number of profiles of the given hour to be ready to process the next hour
            k0 = k0+Nprofiles
        
            # calculating mean profiles of the model profiles collected around radiosondes
            meanProfile_mod    = np.nanmean(MatrixProfiles_mod, axis=1)
            stdProfile_mod     = np.nanstd(MatrixProfiles_mod, axis=1)
            meanProfile_T_mod  = np.nanmean(MatrixProfiles_T_mod, axis=1)
            stdProfile_T_mod   = np.nanstd(MatrixProfiles_T_mod, axis=1)
            meanProfile_RH_mod = np.nanmean(MatrixProfiles_RH_mod, axis=1)
            stdProfile_RH_mod  = np.nanstd(MatrixProfiles_RH_mod, axis=1)
            
            outDict         = {'hour':hourSel, \
                               'Nprofiles':Nprofiles, \
                               'MatrixProfile_radios':MatrixProfiles_radios, \
                               'MatrixProfile_T_radios':MatrixProfiles_T_radios, \
                               'MatrixProfile_RH_radios':MatrixProfiles_RH_radios, \
                               'MatrixHeight_radios':MatrixHeight_radios, \
                               'meanProfile_mod':meanProfile_mod, \
                               'stdProfileMod':stdProfile_mod, \
                               'meanProfile_T_mod':meanProfile_T_mod, \
                               'stdProfile_T_Mod':stdProfile_T_mod, \
                               'meanProfile_RH_mod':meanProfile_RH_mod, \
                               'stdProfile_RH_Mod':stdProfile_RH_mod, \
                               'lcl_rad_hour':arr_lcl_hour_radios, \
                               'lcl_mod_hour':arr_lcl_hour_mod, \
                               'ccl_rad_hour': arr_ccl_hour_radios, \
                               'ccl_mod_hour': arr_ccl_hour_mod, \
                               'lts_rad_hour':arr_lts_hour_radios, \
                               'lts_mod_hour':arr_lts_hour_mod, \
                               'pblHeight_rad_hour':arr_pblHeight_hour_radios, \
                               'pblHeight_mod_hour':arr_pblHeight_mod, \
                               'n_lts_hour':len(arr_lts_hour_mod), \
                               'n_lts_hour_obs':len(arr_lts_hour_radios), \
                               'n_lcl_hour':len(arr_lcl_hour_mod), \
                               'n_lcl_hour_obs':len(arr_lcl_hour_radios)}
            
            listHourDict.append(outDict) 
            # list of dictionaries: each dictionary contains a matrix with profiles of theta and height for that hour
            
    ## interpolating mean profiles and standard deviation from radiosounding on ICON height grid
    gridHeight                    = height_mod[0]# np.arange(0.,8000., 5.)    # grid of 5 m resolution in height
    NgridInterp                   = len(gridHeight)
    MatrixHourMeanProfileThetaRad = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourStdProfileThetaRad  = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourMeanProfileTRad     = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourStdProfileTRad      = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourMeanProfileRHRad    = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourStdProfileRHRad     = np.zeros((NgridInterp, len(listHourDict)))
    MatrixHourMeanProfileThetaRad.fill(np.nan)
    MatrixHourStdProfileThetaRad.fill(np.nan)
    MatrixHourMeanProfileTRad.fill(np.nan)
    MatrixHourStdProfileTRad.fill(np.nan)
    MatrixHourMeanProfileRHRad.fill(np.nan)
    MatrixHourStdProfileRHRad.fill(np.nan)
    
    
    # loop on the hours: for each hour, we read a matrix of height and thetav from radiosondes
    for indHour in range(len(listHourDict)):
        MatrixHeightHour = listHourDict[indHour]['MatrixHeight_radios']
        MatrixHourTheta  = listHourDict[indHour]['MatrixProfile_radios']
        MatrixHour_T     = listHourDict[indHour]['MatrixProfile_T_radios']
        MatrixHour_RH    = listHourDict[indHour]['MatrixProfile_RH_radios']
        
        sizeMatrix       = np.shape(MatrixHourTheta)
        Nradiosondes     = sizeMatrix[1]
        NheightsRad      = sizeMatrix[0]
        MeanProfileTheta = np.zeros((NgridInterp, Nradiosondes))
        MeanProfileTheta.fill(np.nan)
        MeanProfileT     = np.zeros((NgridInterp, Nradiosondes))
        MeanProfileT.fill(np.nan)
        MeanProfileRH    = np.zeros((NgridInterp, Nradiosondes))
        MeanProfileRH.fill(np.nan)
        
        # we loop on heights in radiosondes and we average for each model 
        # heigth grid box, the values of theta foudn in the radiosonde profiles
        for indRadiosonde in range(Nradiosondes):
            for indHeight in range(len(gridHeight)-1):
                Hmax     = gridHeight[indHeight]
                Hmin     = gridHeight[indHeight+1]
                #print(Hmin,Hmax)
                indFound = np.where((MatrixHeightHour[:,indRadiosonde] >= Hmin) *\
                                    (MatrixHeightHour[:,indRadiosonde] < Hmax))
                #print(indFound)
                MeanProfileTheta[indHeight,indRadiosonde] = np.nanmean(MatrixHourTheta[indFound,indRadiosonde])
                MeanProfileT[indHeight,indRadiosonde] = np.nanmean(MatrixHour_T[indFound,indRadiosonde])
                MeanProfileRH[indHeight,indRadiosonde] = np.nanmean(MatrixHour_RH[indFound,indRadiosonde])
                
                
                #print(MatrixHeightHour[indFound,indRadiosonde])
                #print(MatrixHeightHour[:,indRadiosonde])
            # calculating mean profile of theta V for observations
            #mean_theta_V_radiosondes = np.mean()
        #if hourSel == K0:
        MatrixHourMeanProfileThetaRad[:, indHour] = np.nanmean(MeanProfileTheta, axis=1)
        MatrixHourStdProfileThetaRad[:, indHour]  = np.nanstd(MeanProfileTheta, axis=1)
        MatrixHourMeanProfileTRad[:, indHour]     = np.nanmean(MeanProfileT, axis=1)
        MatrixHourStdProfileTRad[:, indHour]      = np.nanstd(MeanProfileT, axis=1)
        MatrixHourMeanProfileRHRad[:, indHour]    = np.nanmean(MeanProfileRH, axis=1)
        MatrixHourStdProfileRHRad[:, indHour]     = np.nanstd(MeanProfileRH, axis=1)
        
    return(MatrixHourMeanProfileThetaRad, MatrixHourStdProfileThetaRad, listHourDict, \
           MatrixHourMeanProfileTRad, MatrixHourStdProfileTRad, MatrixHourMeanProfileRHRad, \
           MatrixHourStdProfileRHRad)


def f_calPBLcloudMask(PBLcloud_dataset,time,height):
    """
    author: Claudia Acquistapace
    date : Friday 19 June 2020
    goal : calculate cloud mask based on cloud base and cloud top time series provided as input. It is assumed that there
    are 8 levels for cloud base/top identification i.e. np.shape(PBLcloud_dataset) = len(time),8
    """
    # defining a new cloud mask based on the PBl cloud bases and tops
    CB_matrix = PBLcloud_dataset.cloudBase.values
    CT_matrix = PBLcloud_dataset.cloudTop.values
    PBLtime   = PBLcloud_dataset.time.values

    PBL_cloudMask = np.zeros((len(time),len(height)))
    for indTime in range(len(time)):
        for indLev in range(8):
            if (~np.isnan(CB_matrix[indTime, indLev])):
                cb_height = CB_matrix[indTime, indLev]
                ct_height = CT_matrix[indTime, indLev]
                ind = (height >= cb_height) * (height <= ct_height)
                PBL_cloudMask[indTime,ind] = 1

    return(PBL_cloudMask)


def f_calculateCloudFractionPBLclouds(PBLcloudmask,time,height,Nmin_string):
    '''
    author: Claudia Acquistapace
    date: friday 19 June 2020
    goal: calculate cloud fraction from cloud mask over a given amount of minutes Nmin
    input: PBLcloudmask: cloud mask matrix
            time: time array corresponding to the cloud mask
            height: height array corresponding to the cloud mask
            Nmin: number of minutes over which to calculate the cloud fraction
    output: CF_dataset = xr.Dataset({'CF': (['time','height'], CF_PBL)},
                        coords = {'time':datetime_CF,
                                  'height':height})
    '''
    # calculating cloud fraction every 15 minutes
    cloudMask_DF = pd.DataFrame(PBLcloudmask, index=time, columns=height)
    datetime_CF = pd.date_range(start=time[0], end=time[-1], freq=Nmin_string+'min')
    datetime_CF = datetime_CF.to_pydatetime()

    CF_PBL = np.zeros((len(datetime_CF), len(height)))
    for indTime in range(len(datetime_CF)-1):
        mask_t = (cloudMask_DF.index > datetime_CF[indTime]) * (cloudMask_DF.index < datetime_CF[indTime+1])
        Selection_cloudMask = cloudMask_DF[mask_t]
        for indHeight in range(len(height)):
            CFArray = Selection_cloudMask.loc[:,Selection_cloudMask.columns[indHeight]]
            CF_PBL[indTime,indHeight] = len(CFArray[CFArray == 1])/len(CFArray)

    CF_dataset = xr.Dataset({'CF': (['time','height'], CF_PBL)},
                        coords = {'time':datetime_CF,
                                  'height':height})
    return(CF_dataset)