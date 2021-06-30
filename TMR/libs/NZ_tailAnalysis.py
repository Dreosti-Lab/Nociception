# -*- coding: utf-8 -*-
"""
Library to measure tail parameters

@author: Adam Kampff & Tom Ryan & Elisa C
Compiled and edited Tom Ryan (Dreosti-Lab,UCL)
"""

## created 23/01/21 TR as library for script S2a

import numpy as np
import pandas as pd 


def computeTailCurvatures(csvXFile,csvYFile,verbose=True):
    if verbose: print('loading xData')
    x_data = np.genfromtxt(csvXFile, delimiter=',')
    if verbose: print('loading yData')
    y_data = np.genfromtxt(csvYFile, delimiter=',')
    # Determine number of frames, segments and angles
    num_frames = np.shape(x_data)[0]
    num_segments = np.shape(x_data)[1]
    num_angles = num_segments - 1
      
    # Allocate space for measurements - creates empty arrays to store data
    cumulAngles = np.zeros(num_frames)
    curvatures = np.zeros(num_frames)
    motion = np.zeros(num_frames)
    
    ## Measure tail motion, angles, and curvature ##
    # Take x and y values of first frame for each segment
    prev_xs = x_data[0, :]
    prev_ys = y_data[0, :]
    
    ########### Start of frame loop (f) #################
    for f in range(num_frames):
        if verbose: print(f)
        delta_thetas = np.zeros(num_angles) # Make an array of zeros with the same size as num angles (num seg-1)
        prev_theta = 0.0 # set first theta to zero
        
        ############### Start of segment loop (a) ################
        for a in range(num_angles):
            dx = x_data[f, a+1] - x_data[f, a] # dx between each segment for the same frame
            dy = y_data[f, a+1] - y_data[f, a] # dy between each segment for the same frame
            theta = np.arctan2(dx, dy) * 360.0 / (2.0*np.pi) # calc arctangent bt dx and dy, convert to deg
            delta_thetas[a] = theta - prev_theta
            prev_theta = theta # prev theta is set to current theta
        ############### End of angle loop (a) ################
        
        cumulAngles[f] = np.sum(delta_thetas) # sum all angles for this frame
        curvatures[f] = np.mean(np.abs(delta_thetas)) # mean of abs value of angles
    
        # Measure motion
        diff_xs = x_data[f,:] - prev_xs # difference between current x and prev x from each segment and each frame
        diff_ys = y_data[f,:] - prev_ys
        motion[f] = np.sum(np.sqrt(diff_xs*diff_xs + diff_ys*diff_ys)) # motion as sqrt of sq diff of x & y
        
        # Store previous tail
        prev_xs = x_data[f, :] # So that we don't always take the 1st frame of the movie to calculate motion
        prev_ys = y_data[f, :]
    ####### End of frame loop (f) ###############
    return cumulAngles,curvatures,motion

def currentVsPrevResponses(volts,freqResp):
    sortedVolts=volts[np.argsort(volts)]
    currentVsPrev=[]
    for i in range(1,len(volts)):
        thisVolt=volts[i]
        preVolt=volts[i-1]
        allCurrentVsPrev.append([thisVolt,preVolt])
        thisResp=freqResp[np.where(sortedVolts==thisVolt)[0][0]]
        if thisResp==1:
            prevResp=freqResp[np.where(sortedVolts==preVolt)[0][0]]
            if prevResp==1:
                currentVsPrev.append([thisVolt,preVolt])
    return currentVsPrev
                
                
def findStartEndAndCutEC(volt_data,uncut_data,preStimFrames=400,postStimFrames=4000,sort=True):
    diffvolt = np.diff(volt_data)
    starts = np.where(diffvolt > 0)[0] #because that yields a tuple
    starts = starts + 1
    volts = volt_data[starts]
    num_trials = len(starts)
    startFrame = starts - preStimFrames
    endFrame = starts + postStimFrames-1
    
    cut_data = pd.DataFrame()
    for t in range(num_trials):
        cut_data[t]=uncut_data[startFrame[t]:endFrame[t]]
    
    cut_data = cut_data.T
    cut_data['Voltage']=volts
    if sort:
        cut_data = cut_data.sort_values(by=['Voltage'])
    return cut_data

def computeEvokedBouts(Data_mot,Data_CA,k,l,preStimFrames,responseFrames,reflexFrames,interbout_interv,threshold_down,threshold_up,avg_BL,SD_BL,msPerFrame,num_stim,conf_interv=0.95):
    num_stim=np.int(num_stim) # for some strange reason
    frame_cutoff = preStimFrames + reflexFrames 
    
    timeResp_PerI = []
    reflex_timeResp_PerI = []
    veryFirst_timeResp_PerI_all = []
    NumResp_PerI = []
    freqResp_PI = []
    reflex_freqResp_PI=[]
    velocities = [[] for _ in range(num_stim)] # mke space 4 list of lists
    avgVelo_PerI = []
    firstVelo_PerI = []
    reflex_firstVelo_PerI = []
    
    # Go line by line & frame by frame; find 1st time Data[f]>threshold
    # follows Data[f]<threshold -> 1st response -> calc latency; after
    # that, Data[f] must go to BL for CI*interbout_interv frames before
    # another response can be found; use beginning/end response frames to
    # define bout, calculate mean velocity for each bout
    
    for n in range(num_stim):
        count_resp = 0
        first_resp = 0
        reflex_resp = 0
        frame = preStimFrames
        #reflex_time_resp = 0
        while frame < responseFrames:
            if (Data_mot[n,frame-1] < threshold_up and Data_mot[n,frame] > threshold_up):
                begin_resp=frame
                ready=False
                while Data_mot[n,frame] >= threshold_down and frame<responseFrames: # wait til <thresh_d
                    frame+=1          
                count_below = 0
                while ready==False and frame<responseFrames: 
                    frame_max=frame+interbout_interv
                    if frame_max > responseFrames:
                        frame_max=responseFrames  # so we don't "overrun"  
                    while frame < frame_max:
                        if Data_mot[n,frame] < threshold_down:
                            count_below+=1
                        frame+=1
                    if (count_below / interbout_interv) >= conf_interv:
                        end_resp=frame-interbout_interv
                        ready=True
                        count_resp+=1
                        bout = Data_CA[n,begin_resp:end_resp]
                        velocity_f=np.mean(np.abs(np.diff(bout)))
                        velocity_ms = velocity_f/msPerFrame
                        velocities[n].append(velocity_ms)
                        if count_resp == 1: # first response
                            if begin_resp < frame_cutoff: # reflex response
                                reflex_resp=begin_resp
                                reflex_firstVelo_PerI.append(velocity_ms)
                            else:  
                                first_resp=begin_resp # regular response
                                firstVelo_PerI.append(velocity_ms)
                                reflex_firstVelo_PerI.append(None)
                        if count_resp == 2 and reflex_resp != 0:
                            first_resp=begin_resp
                            firstVelo_PerI.append(velocity_ms)
                            # (if 1st resp is reflex, find 1st reg resp)                
            frame+=1
       
        NumResp_PerI.append(count_resp) # total # responses per intensity
        avgVelocity=np.mean(velocities[n])
        avgVelo_PerI.append(avgVelocity)
         
        if count_resp == 1 and reflex_resp != 0: # only reflex response
            time_resp = None
            firstVelo_PerI.append(None)
        
        if count_resp == 0:
            time_resp = None # latency > 10 secs doesn't count as response
            reflex_time_resp = None
            velocities[n] = [None]
            firstVelo_PerI.append(None)
            reflex_firstVelo_PerI.append(None)
            resp=0
            ref_resp=0
        else:
            time_resp = (first_resp-preStimFrames)*msPerFrame # Convert to ms, t=0 -> stim starts
            resp=1
            if reflex_resp !=0: # reflex resp found
                reflex_time_resp=(reflex_resp-preStimFrames)*msPerFrame
                ref_resp = 1
            else:
                reflex_time_resp = None
                ref_resp = 0
        
        timeResp_PerI.append(time_resp) # time 1st response per intensity
        reflex_timeResp_PerI.append(reflex_time_resp) # time reflex resp per intensity
        
        # Find min between reflex and 1st resp i.e. very 1st resp
        merge_timeResp = [time_resp,reflex_time_resp]
        merge_timeResp_ar = np.array(merge_timeResp, dtype=np.float) # convert to array because of None
        veryFirst_timeResp_PerI = np.nanmin(merge_timeResp_ar)
        veryFirst_timeResp_PerI_all.append(veryFirst_timeResp_PerI)
        
        # Freq of response per intensity
        freqResp_PI.append(resp)
        reflex_freqResp_PI.append(ref_resp)
        
    return timeResp_PerI,reflex_timeResp_PerI,veryFirst_timeResp_PerI_all,freqResp_PI,reflex_freqResp_PI,NumResp_PerI,velocities,avgVelo_PerI,firstVelo_PerI,reflex_firstVelo_PerI

def findStartAndSetBL(volt_data,Data,k=2):
    diffvolt = np.diff(volt_data)
    start = np.where(diffvolt > 0)[0] # Every time laser is switched on
    BL_start = start[0] # First time laser is switched on
    # Data = np.genfromtxt(uncut_smooth_data_path, delimiter=',')
    
    i = 0
    avg_BL = 0
    SD_BL = 0
    max_BL = 0
    threshold_low = avg_BL + k*SD_BL 
    
    # Go back in 500ms (200 frames) intervals until you find 1 w/out movement
    while max_BL >= threshold_low:
        i+=1
        BL = Data[BL_start-40*i:BL_start-40*(i-1)] # 200 frames b4 stim on
        avg_BL = np.mean(BL)
        SD_BL = np.std(BL)
        max_BL = np.amax(BL)
        threshold_low = avg_BL + k*SD_BL
                
    return avg_BL,SD_BL


def MaxTimeDF(cut_data,FrameOn,respFreq):   
    cut_data_2 = abs(cut_data.drop('Voltage',axis=1))
    maxRound = cut_data_2.max(axis=1)
    maxRound.reset_index(drop=True,inplace=True) # reset indices into 0 - 5
    frameMax = cut_data_2.idxmax(axis=1)
    frameMax.reset_index(drop=True,inplace=True) # reset indices into 0 - 5
    timeMaxRound = (frameMax-FrameOn)/400*1000
    
    timeMaxRounds = []
    maxRounds = []
    for r in range(len(respFreq)):
        if respFreq[r] == 1:
            timeMaxRounds.append(timeMaxRound[r])
            maxRounds.append(maxRound[r])
        else:
            timeMaxRounds.append(None)
            maxRounds.append(None)
    
    return maxRounds, timeMaxRounds


def MaxTimeParam(Data,num_stim,respFreq):   
    maxRound = []
    boutMaxRound = []

    for i in range(num_stim): # loop through stim
        if Data[i]==[]:Data[i]=[None]
        maxInt = np.nanmax(Data[i])
        boutMax = np.nanargmax(Data[i])
        if respFreq[i] == 1:
            maxRound.append(maxInt)
            boutMaxRound.append(boutMax)
        else:
            maxRound.append(None)
            boutMaxRound.append(None)

    return maxRound, boutMaxRound


def MaxTimeReflexDF(cut_data,FrameOn,ReflexThresh,reflexRespFreq):   
    # Find max (& TMax) CA, curv, mot in the "reflex response" interval
    cut_data_2 = abs(cut_data.drop('Voltage',axis=1))
    cut_data_3 = cut_data_2.iloc[:,FrameOn:(FrameOn+ReflexThresh+1)] # window from laser on to cutoff
    maxRound = cut_data_3.max(axis=1)
    maxRound.reset_index(drop=True,inplace=True) # reset indices into 0 - 5
    frameMax = cut_data_3.idxmax(axis=1)
    frameMax.reset_index(drop=True,inplace=True) # reset indices into 0 - 5
    timeMaxRound = (frameMax-FrameOn)/400*1000
    
    timeMaxRounds = []
    maxRounds = []
    for r in range(len(reflexRespFreq)):
        if reflexRespFreq[r] == 1:
            timeMaxRounds.append(timeMaxRound[r])
            maxRounds.append(maxRound[r])
        else:
            timeMaxRounds.append(None)
            maxRounds.append(None)
    
    return maxRounds, timeMaxRounds

