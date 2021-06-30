# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:11:12 2021

@author: Adam Kampff & Tom Ryan & Elisa C
Compiled and edited Tom Ryan (Dreosti-Lab,UCL)
"""
import numpy as np

def AverageRoundsDict(Files):   
    # Standard for raw data from dictionary (mot, curv, cumulAng)
    # load first dictionary and find the num_stim
    dict1=np.load(Files[0],allow_pickle=True).item()
    num_rounds=len(Files)
    del dict1 # clear from memory
    preAverages=[]
    averages=[]    
    sds=[]
    ## loop through dictionaries and collect data in lists
    ## N.B. some of these nested loops would be better and maybe faster as list comprehension sttements but these are easier to read
    for i,File in enumerate (Files): # loop through files (each round)
        Data=np.load(File,allow_pickle=True).item()
        preAvg=[]
        ks=list(Data.keys())[:-1] # last key (should be metadata)... can be explicit about this needs be
        for j in ks: # last key is very first response (only one value per round)
            if isinstance(Data[j][0],list): # the way velocities is collated in step 2 is awkward and builds lists of lists of lists of single numbers, unlike the others. This is a hack until I can disentagle the original problem (if we bother)
                prepreAvg=[]
                for k in range(len(Data[j])):
                    prepreAvg.append(Data[j][k][0])
                preAvg.append(prepreAvg)
    #            del prepreAvg
            else:
                preAvg.append(Data[j])
        preAverages.append(preAvg)
    #del Data, meta, preAvg # clear from memory
    # Build average list
    for i,k in enumerate(ks): # loop through the keys
        temp=[]
        for j in range(num_rounds):
            temp.append(preAverages[j][i])
            # cycle through and check for 'NoneType' or None made in pandas
            for l,chk in enumerate(temp[j]):
                #  turn them to nans
                if chk is None: temp[j][l]=np.nan
        averages.append(np.nanmean(temp,axis=0))
        sds.append(np.nanstd(temp,axis=0))
        
    return ks,averages,sds
        
        
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