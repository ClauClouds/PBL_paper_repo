#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:15:38 2019

@author: cacquist
"""

#%%
# processing of LWC profile: goal is to derive hourly mean LWC profile from the collected clouds for model and obs
LWC_All_obs              = np.concatenate(LWC_obs, axis=0)
LWC_All_mod              = np.concatenate(LWC_mod, axis=0)
duration_all_obs         = np.concatenate(duration_obs, axis=0)
duration_all_mod         = np.concatenate(duration_mod, axis=0)
CB_all_obs               = np.concatenate(LWC_obs, axis=0)
meanheightFromCB_All_obs = np.concatenate(meanheightFromCB_obs, axis=0)
meanheightFromCB_All_mod = np.concatenate(meanheightFromCB_mod, axis=0)
cloudMeanCT_All_obs      = np.concatenate(cloudMeanCB_obs, axis=0)
cloudMeanCB_All_obs      = np.concatenate(cloudMeanCB_obs, axis=0)
cloudMeanCT_All_mod      = np.concatenate(cloudMeanCT_mod, axis=0)
cloudMeanCB_All_mod      = np.concatenate(cloudMeanCB_mod, axis=0)
#duration_intervals = [0, 400, 800, 1600, 3200]
duration_intervals       = [0, 200, 800, 1600]

LWC_obs_DF              = pd.DataFrame(LWC_All_obs, index=duration_all_obs, columns=height)
LWC_mod_DF              = pd.DataFrame(LWC_All_mod, index=duration_all_mod, columns=height)
meanHeightFromCB_obs_DF = pd.DataFrame(meanheightFromCB_All_obs, index=duration_all_obs, columns=height)
meanHeightFromCB_mod_DF = pd.DataFrame(meanheightFromCB_All_mod, index=duration_all_mod, columns=height)
CB_obs_DF               = pd.Series(cloudMeanCB_All_obs, index=duration_all_obs)
CB_mod_DF               = pd.Series(cloudMeanCB_All_mod, index=duration_all_mod)
CT_obs_DF               = pd.Series(cloudMeanCT_All_obs, index=duration_all_obs)
CT_mod_DF               = pd.Series(cloudMeanCT_All_mod, index=duration_all_mod)


#%%
# calculating mean profiles: concept is to select all profiles within the \
# specific duration interval. Then, we select all heights between 0 and 1 and we loop on them:
# for every height 
height_grid = np.linspace(0., 1., num=5)
Nprofiles_durationIntervalsObs = []
Nprofiles_durationIntervalsMod = []

LWC_mean_rescaled_obs = []
LWC_mean_rescaled_mod = []
CB_obs_distr_arr = []
CB_mod_distr_arr = []
CT_obs_distr_arr = []
CT_mod_distr_arr = []
for indDuration in range(0, len(duration_intervals)):
    #print(indDuration)
    #print(indDuration-len(duration_intervals))
    if ((indDuration-len(duration_intervals)+1) == 0.):
        mask_duration_obs    = LWC_obs_DF.index > duration_intervals[indDuration]
        mask_duration_mod    = LWC_mod_DF.index > duration_intervals[indDuration]
        #print('interval sampled > '+str(duration_intervals[indDuration])+' s')
    else:
        mask_duration_obs    = (LWC_obs_DF.index > duration_intervals[indDuration]) * (LWC_obs_DF.index <= duration_intervals[indDuration+1])
        mask_duration_mod    = (LWC_mod_DF.index > duration_intervals[indDuration]) * (LWC_mod_DF.index <= duration_intervals[indDuration+1])
        #print('interval sampled'+str(duration_intervals[indDuration])+'-'+str(duration_intervals[indDuration+1])+' s')
    
    # selecting the relative heights and the LWC corresponding profiles of the selected profiles     
    relHeightCBCT_obs_durInt = meanHeightFromCB_obs_DF.loc[mask_duration_obs,:].values
    LWC_obs_durInt           = LWC_obs_DF.loc[mask_duration_obs,:].values
    relHeightCBCT_mod_durInt = meanHeightFromCB_mod_DF.loc[mask_duration_mod,:].values
    LWC_mod_durInt           = LWC_mod_DF.loc[mask_duration_mod,:].values    
    # reading the number of profiles corresponding to the selected duration 
    # interval and storing the value
    N_obs                    = np.shape(relHeightCBCT_obs_durInt)[0]
    Nprofiles_durationIntervalsObs.append(N_obs)
    N_mod                    = np.shape(relHeightCBCT_mod_durInt)[0]
    Nprofiles_durationIntervalsMod.append(N_mod)    
    
    # building array to use as index given by the duration selected in the picked interval
    durationSelected_obs     = duration_all_obs[mask_duration_obs]
    durationSelected_mod     = duration_all_mod[mask_duration_mod]
    
    # defining pandas structures for filtering and averaging the data
    RescaledLWCMatrixObs     = np.zeros((N_obs,len(height_grid)))
    RescaledLWCMatrixMod     = np.zeros((N_mod,len(height_grid)))
    
    # loop on heights: for every height and every profile, we check the heights 
    # between the two height grid points and we average them and store them in the rescaled matrix
    for indHeight in range(len(height_grid)-1):
        Hmin = height_grid[indHeight]
        Hmax = height_grid[indHeight+1]
        
        for indTimeProf in range(len(durationSelected_obs)):
            ProfSel_obs                                 = pd.Series(LWC_obs_durInt[indTimeProf,:], \
                                                index=relHeightCBCT_obs_durInt[indTimeProf,:])
            mask_h                                      = (ProfSel_obs.index >= Hmin) * (ProfSel_obs.index <= Hmax)
            RescaledLWCMatrixObs[indTimeProf,indHeight] = np.nanmean(ProfSel_obs.loc[mask_h])

        for indTimeProf in range(len(durationSelected_mod)):
            ProfSel_mod                                 = pd.Series(LWC_mod_durInt[indTimeProf,:], \
                                                index=relHeightCBCT_mod_durInt[indTimeProf,:])
            mask_h                                      = (ProfSel_mod.index >= Hmin) * (ProfSel_mod.index <= Hmax)
            RescaledLWCMatrixMod[indTimeProf,indHeight] = np.nanmean(ProfSel_mod.loc[mask_h])
    
    
    LWC_mean_rescaled_obs.append(RescaledLWCMatrixObs)
    LWC_mean_rescaled_mod.append(RescaledLWCMatrixMod)
    
    
    # calculating CB/CT distributions for each of the duration types
    CB_obs_distr_arr.append(CB_obs_DF.loc[mask_duration_obs])
    CB_mod_distr_arr.append(CB_mod_DF.loc[mask_duration_mod])
    CT_obs_distr_arr.append(CT_obs_DF.loc[mask_duration_obs])
    CT_mod_distr_arr.append(CT_mod_DF.loc[mask_duration_mod])
    
# calculating mean profiles for observations on the rescaled grid
LWC_mean_prof_obs = np.zeros((len(duration_intervals),len(height_grid)))
LWC_std_prof_obs  = np.zeros((len(duration_intervals),len(height_grid)))

for indDurInt in range(len(LWC_mean_rescaled_obs)):
    LWC_mean_prof_obs[indDurInt,:] = np.nanmean(LWC_mean_rescaled_obs[indDurInt], axis=0)
    LWC_std_prof_obs[indDurInt,:]  = np.nanstd(LWC_mean_rescaled_obs[indDurInt], axis=0)
    
LWC_mean_prof_mod = np.zeros((len(duration_intervals),len(height_grid)))
LWC_std_prof_mod  = np.zeros((len(duration_intervals),len(height_grid)))

for indDurInt in range(len(LWC_mean_rescaled_mod)):
    LWC_mean_prof_mod[indDurInt,:] = np.nanmean(LWC_mean_rescaled_mod[indDurInt], axis=0)
    LWC_std_prof_mod[indDurInt,:]  = np.nanstd(LWC_mean_rescaled_mod[indDurInt], axis=0)   
    
#%%
   
# plotting CB distributions
fig = plt.figure(figsize=(10,14))
plt.gcf().subplots_adjust(bottom=0.15)
stringplot = ['a)','b)','c)','d)']
stringmidHours = ['0-200 s', '200 - 800 s', '800-1600s', '>1600 s']
Rangearr_CB = [[600., 3000.], [600., 3000.],[600., 3000.], [600., 3000.]]
Rangearr_CT = [[1300., 5000.], [1300., 5000.],[1300.,5000.], [1300., 5000.]]
Nbins_CB = 10
Nbins_CT = 10
ymaxArr = [0.0025, 0.0009, 0.0015, 0.025]
xlabelArr = [1400., 5000., 1000., 250.]
ylabelArr = [0.007, 0.0022, 0.0037, 0.056]
colorArr = ['black', 'black', 'black', 'black']
indObs = [1,3,5,7]
indMod = [2,4,6,8]
for indPlot in range(0,4):

    ax = plt.subplot(4,2,indObs[indPlot])  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('cloud base height [m]') 
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')

    plt.hist(CB_obs_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    plt.hist(CB_mod_distr_arr[indPlot], \
                      bins=Nbins_CB, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CB[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    ax = plt.subplot(4,2,indMod[indPlot])  
    ax.set_title(stringmidHours[indPlot])
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left() 
             #ax.text(xlabelArr[indPlot-1], ylabelArr[indPlot-1], stringplot[indPlot-1], fontsize=15)
         #        matplotlib.rc('xtick', labelsize=12)                        # sets dimension of ticks in the plots
         #        matplotlib.rc('ytick', labelsize=12)                        # sets dimension of ticks in the plots
    ax.set_ylabel('occurrences')
    ax.set_xlabel('cloud top height [m]')   
    plt.hist(CT_obs_distr_arr[indPlot], \
                      bins=10, \
                      normed=True, \
                      color='black', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)    
    plt.hist(CT_mod_distr_arr[indPlot], \
                      bins=10, \
                      normed=True, \
                      color='red', \
                      range=Rangearr_CT[indPlot-1], \
                      cumulative=False, \
                      alpha=0.5)  
    ax.set_title(stringmidHours[indPlot]+\
                 ' (Nobs = '+str(Nprofiles_durationIntervalsObs[indPlot])+','+\
                 ' Nmod = '+str(Nprofiles_durationIntervalsMod[indPlot])+')')

plt.tight_layout()

#plt.savefig(pathDebugFig+'histograms_cloudProperties_'+date+'.png', format='png')
#%%
# plotting all profiles for each duration interval after regridding
fig, ax       = plt.subplots(nrows=2, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.002, 0.002, 0.0004, 0.001]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr_obs = [1,2,3,4]
indPlotArr_mod = [5,6,7,8]

for ind in range(0, len(LWC_mean_rescaled_obs)):
     xmax          = xmaxArr[ind]
     ax        = plt.subplot(2,len(duration_intervals),indPlotArr_obs[ind])  
     ax.spines["top"].set_visible(False)  
     ax.spines["right"].set_visible(False)  
     ax.get_xaxis().tick_bottom()  
     ax.get_yaxis().tick_left() 
     ax.set_ylim(ymin, ymax)
     ax.set_xlim(xmin, xmax)
     ax.set_xlabel('LWC $Kg m^{-3}$')
     ax.set_ylabel('rel. dist. from cloud base')
     LWC_plot_obs = LWC_mean_rescaled_obs[ind]
     for indProfile in range(np.shape(LWC_mean_rescaled_obs[ind])[0]):
         plt.plot(LWC_plot_obs[indProfile,:], height_grid, color='black')
         
         
     ax1        = plt.subplot(2,len(duration_intervals),indPlotArr_mod[ind])  
     ax1.spines["top"].set_visible(False)  
     ax1.spines["right"].set_visible(False)  
     ax1.get_xaxis().tick_bottom()  
     ax1.get_yaxis().tick_left() 
     ax1.set_ylim(ymin, ymax)
     ax1.set_xlim(xmin, xmax)
     ax1.set_xlabel('LWC $Kg m^{-3}$')
     ax1.set_ylabel('rel. dist. from cloud base')
     LWC_plot_mod = LWC_mean_rescaled_mod[ind]
     for indProfile in range(np.shape(LWC_mean_rescaled_mod[ind])[0]):
         plt.plot(LWC_plot_mod[indProfile,:], height_grid, color='red')        
plt.tight_layout()
#%%
# =============================================================================
fig, ax       = plt.subplots(nrows=1, ncols=len(duration_intervals), figsize=(10,5))
# #matplotlib.rcParams['savefig.dpi'] = 300
plt.gcf().subplots_adjust(bottom=0.15)
fig.tight_layout()
ymax          = 1.
ymin          = 0.
xmin          = 0.
xmaxArr       = [0.002, 0.002, 0.0004, 0.001]
fontSizeTitle = 16
fontSizeX     = 12
fontSizeY     = 12
indPlotArr = [1,2,3,4]
for ind in range(0, len(LWC_mean_rescaled_obs)):
     xmax          = xmaxArr[ind]
     ax        = plt.subplot(1,len(duration_intervals),indPlotArr[ind])  
     ax.spines["top"].set_visible(False)  
     ax.spines["right"].set_visible(False)  
     ax.get_xaxis().tick_bottom()  
     ax.get_yaxis().tick_left() 
     ax.set_ylim(ymin, ymax)
     ax.set_xlim(xmin, xmax)
     ax.set_xlabel('LWC $Kg m^{-3}$')
     ax.set_ylabel('rel. dist. from cloud base')
     plt.plot(LWC_mean_prof_obs[ind,:], height_grid, color='black', label='obs')
     y1        = LWC_mean_prof_obs[ind,:]-LWC_std_prof_obs[ind,:]
     y2        = LWC_mean_prof_obs[ind,:]+LWC_std_prof_obs[ind,:]
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='black', alpha=0.2)
     
     plt.plot(LWC_mean_prof_mod[ind,:], height_grid, color='red', label='icon-lem')
     y1        = LWC_mean_prof_mod[ind,:]-LWC_std_prof_mod[ind,:]
     y2        = LWC_mean_prof_mod[ind,:]+LWC_std_prof_mod[ind,:]
     plt.fill_betweenx(height_grid, y1, y2, where=y2>y1, facecolor='red', alpha=0.2)
# =============================================================================

    