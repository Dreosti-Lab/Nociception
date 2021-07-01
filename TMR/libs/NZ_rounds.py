# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:11:12 2021

@author: Adam Kampff & Tom Ryan & Elisa C
Compiled and edited Tom Ryan (Dreosti-Lab,UCL)
"""
import numpy as np
from scipy import stats

def AverageFishDict(Files,key,num_stim,preDrug=True,diff=False):
    # Standard for raw data (mot, curv, cumulAng)
    avgAllFish=[]
    stdevAllFish=[]
    semAllFish=[]
    if preDrug:
        k='preDrug'
    else:
        k='postDrug'
        if diff:
            k='PercentDiff'
        
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each fish)
            Data = np.load(File,allow_pickle=True).item() # grab data from the file
            Data=Data[k][key]['Mean']
            Data=np.array(Data)
            preAvg.append(Data[i,:])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
   
        meanOfTrials=np.nanmean(preAvg,0)
        avgAllFish.append(meanOfTrials)   
        stdevAllFish.append(np.std(preAvg,0))
        semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return preAvg, avgAllFish, stdevAllFish, semAllFish

def AverageFishParamDict(Files,key,num_stim,preDrug=True,diff=False):   
    # For processed data e.g. latency
    avgAllFish=[]
    stdevAllFish=[]
    semAllFish=[]
    if preDrug:
        k='preDrug'
    else:
        k='postDrug'
        if diff:
            k='PercentDiff'
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each fish)
            Data = np.load(File,allow_pickle=True).item() # grab data from the file
            if diff:
                Data=Data[k][key]
            else: 
                Data=Data[k][key]['Mean']
            Data=np.array(Data)
            preAvg.append(Data[i])  # populate list with ith row from each data file (not very efficient as opens each file many many times)

        meanOfTrials=np.nanmean(preAvg,0)
        avgAllFish.append(meanOfTrials)   
        stdevAllFish.append(np.nanstd(preAvg,0))
        semAllFish.append(stats.sem(preAvg,axis=0,ddof=0,nan_policy='omit'))
       
    return preAvg, avgAllFish, stdevAllFish, semAllFish


def AverageFishSingleDict(Files,key,preDrug=True,diff=False):   
    # For fish freq resp (append "doesnt like" excel file with single row?)
    avgAllFish=[]
    stdevAllFish=[]
    semAllFish=[]
    preAvg=[]
    if preDrug:
        k='preDrug'
    else:
        k='postDrug'
        if diff:
            k='PercentDiff'
    for x,File in enumerate (Files): # loop through files (each fish)
        Data_1 = np.load(File, allow_pickle=True).item() # grab data from the file
        if diff:
            Data_1=Data_1[k][key]
        else:
            Data_1=Data_1[k][key]['Mean']   
        preAvg.append(Data_1)  # populate list with ith row from each data file (not very efficient as opens each file many many times)

    meanOfTrials=np.nanmean(preAvg,0)
    avgAllFish.append(meanOfTrials)   
    stdevAllFish.append(np.std(preAvg,0))
    semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return preAvg, avgAllFish, stdevAllFish, semAllFish

def AverageRoundDict(Files,key,num_stim):
    # Standard for raw data (mot, curv, cumulAng)
    avgAllRounds=[]
    stdevAllRounds=[]
    semAllRounds=[]
        
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each fish)
            Data = np.load(File,allow_pickle=True).item() # grab data from the file
            Data=np.array(Data[key])
            preAvg.append(Data[i,:-1])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
            
        meanOfTrials=np.nanmean(preAvg,axis=0)
        avgAllRounds.append(meanOfTrials)   
        stdevAllRounds.append(np.std(preAvg,axis=0))
        semAllRounds.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return avgAllRounds, semAllRounds

def AverageRoundParamDict(Files,key,num_stim):   
    # For processed data e.g. latency
    avgAllRounds=[]
    stdevAllRounds=[]
    semAllRounds=[]
 
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each fish)
            Data = np.load(File,allow_pickle=True).item() # grab data from the file
            Data=np.array(Data[key])
            preAvg.append(Data[i])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
        
        for l,chk in enumerate(preAvg):
        #  turn them to nans
            if chk is None: preAvg[l]=np.nan
        
        meanOfTrials=np.nanmean(preAvg,0)
        avgAllRounds.append(meanOfTrials)   
        stdevAllRounds.append(np.nanstd(preAvg,0))
        semAllRounds.append(stats.sem(preAvg,axis=0,ddof=0,nan_policy='omit'))
       
    return avgAllRounds, semAllRounds


def AverageRoundParamSingleDict(Files,key):   
    # For fish freq resp (append "doesnt like" excel file with single row?)
    avgAllRounds=[]
    stdevAllRounds=[]
    semAllRounds=[]
    
    preAvg=[]
    for x,File in enumerate (Files): # loop through files (each fish)
        Data_1 = np.load(File, allow_pickle=True).item() # grab data from the file
        Data_1=Data_1[key]
        preAvg.append(Data_1)  # populate list with ith row from each data file (not very efficient as opens each file many many times)
        
    for l,chk in enumerate(preAvg):
        #  turn them to nans
        if chk is None: preAvg[l]=np.nan
        
    meanOfTrials=np.mean(preAvg,0)
    avgAllRounds.append(meanOfTrials)   
    stdevAllRounds.append(np.std(preAvg,0))
    semAllRounds.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return avgAllRounds, semAllRounds


def AverageRoundsDict(Files):   
    # Standard for raw data from dictionary (mot, curv, cumulAng)
    # load first dictionary and find the num_stim
    averages=[]    
    sems=[]
#    ## loop through dictionaries and collect data in lists
#    ## N.B. some of these nested loops would be better and maybe faster as list comprehension sttements but these are easier to read
    Data=np.load(Files[0],allow_pickle=True).item()
    ks=list(Data.keys()) # last key (should be metadata)... can be explicit about this needs be
    ks.remove('MetaData')
    ks.remove('VelocityPerBout')
    for j in ks: # last key is very first response (only one value per round)
        try:
            Data[j].shape
            avgAllRounds, semAllRounds = AverageRoundDict(Files,j,6)
        except:
            if len(Data[j])<2:
                avgAllRounds, semAllRounds = AverageRoundParamSingleDict(Files,j)   
            else:
                avgAllRounds, semAllRounds = AverageRoundParamDict(Files,j,6)
        averages.append(avgAllRounds)
        sems.append(semAllRounds)
                        
    return ks,averages,sems
#            if isinstance(Data[j][0],list): # the way velocities is collated in step 2 is awkward and builds lists of lists of lists of single numbers, unlike the others. This is a hack until I can disentagle the original problem (if we bother)
#                prepreAvg=[]
#                for k in range(len(Data[j])):
#                    prepreAvg.append(Data[j][k][0])
#                preAvg.append(prepreAvg)
#    #            del prepreAvg
#            else:
#                preAvg.append(Data[j])
#        preAverages.append(preAvg)
#    #del Data, meta, preAvg # clear from memory
#    # Build average list
#    for i,k in enumerate(ks): # loop through the keys
#        temp=[]
#        for j in range(num_rounds):
#            temp.append(preAverages[j][i])
#            # cycle through and check for 'NoneType' or None made in pandas
#            for l,chk in enumerate(temp[j]):
#                #  turn them to nans
#                if chk is None: temp[j][l]=np.nan
#        averages.append(np.nanmean(temp,axis=0))
#        sds.append(np.nanstd(temp,axis=0))
#        
#    return ks,averages,sds
        
        
#        Data = np.genfromtxt(File, delimiter=',') # grab data from the file
#        preAvg.append(Data[i,:-1])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
#
#    meanOfTrials=np.mean(preAvg,0)
#    avgFishRound.append(meanOfTrials)   
#    stdevFishRound.append(np.std(preAvg,0))
#    volts.append(Data[i,-1])
#       
#    return volts, preAvg, avgFishRound, stdevFishRound

def AverageRounds(Files,num_stim):   
    # Standard for raw data (mot, curv, cumulAng)
    volts=[]
    avgFishRound=[]
    stdevFishRound=[]
   
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each round)
            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
            preAvg.append(Data[i,:-1])  # populate list with ith row from each data file (not very efficient as opens each file many many times)

        meanOfTrials=np.mean(preAvg,0)
        avgFishRound.append(meanOfTrials)   
        stdevFishRound.append(np.std(preAvg,0))
        volts.append(Data[i,-1])
       
    return volts, preAvg, avgFishRound, stdevFishRound

def grabCSVandAverage(path,whichRounds,boo):
    
    files=[fn for fn in glob.glob(path) if os.path.basename(fn).startswith(tuple(whichRounds))]
    if boo==0:
        _, _, avg, stdev = AverageRounds(files,6)
    elif boo==1:
        _, avg, stdev = AverageRoundsSingle(files)
    elif boo==2:
        _, avg, stdev = AverageRoundsParam(files,6)
        
    return avg,stdev