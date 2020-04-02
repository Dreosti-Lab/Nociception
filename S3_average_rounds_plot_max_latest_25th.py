# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:24:14 2020

Average data from csv files within fish, determine peak and plot

@author: Elisa C & Tom Ryan
"""

import os
#import cv2
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import glob
#import timeit
#import csv as csv


### Helper Functions

def readFolderList(folderListFile):
    folderFile=open(folderListFile, 'r')
    folderList=folderFile.readlines()
    data_path=folderList[0][:-1]
    folderList=folderList[1:]
    folderNames=[]
    for i,f in enumerate(folderList):
        stringline = f[:].split()
        expFolderName=data_path + stringline[0]
        folderNames.append(expFolderName) # create a list of folders from folderlistfile
        
    return data_path,folderNames
  
       
def AverageRounds(Files,num_stim):   
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


def MaxAvgTimeRounds(Files,num_stim):   
    volts = []
    maxRound = []
    avgMaxRd = []
    stdevMaxRd = []
    timeMaxRound = []
    avgTimeMaxRd = []
    stdevTimeMaxRd = []
   
    for i in range(num_stim): # loop through stim
        temp=[]
        for x,File in enumerate (Files): # loop through files (each round)
            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
            temp.append(Data[i,:-1])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
            maxInt = np.amax(temp)
            maxRound.append(maxInt)
            frameMax = np.argmax(temp)
            timeMax = frameMax/400*1000
            timeMaxRound.append(timeMax)
           
        meanOfMaxTrials=np.mean(maxRound,0)
        avgMaxRd.append(meanOfMaxTrials)
        stdevMaxRd.append(np.std(maxRound,0))
        meanOfTimeMaxTrials=np.mean(timeMaxRound,0)
        avgTimeMaxRd.append(meanOfTimeMaxTrials)
        stdevTimeMaxRd.append(np.std(timeMaxRound,0))
        volts.append(Data[i,-1])
       
    return volts, maxRound, avgMaxRd, stdevMaxRd, timeMaxRound, avgTimeMaxRd, stdevTimeMaxRd


def GetSdTraces(avData,sdData):
    all_neg_SD = []
    all_pos_SD = []
    
    for j in range(len(avData)):
        neg_SD = np.subtract(avData[j],sdData[j])
        pos_SD = np.add(avData[j],sdData[j])
        all_neg_SD.append(neg_SD)
        all_pos_SD.append(pos_SD)
    
    return all_neg_SD, all_pos_SD


######## Plotting functions #########
def PlotRoundsPerInt(Files,ylabel,Figname,xmax):
    # Plot individual rounds per intensity
    figure=plt.figure(figsize=(10,8))
    for q in range(6): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files):
            Data = np.genfromtxt(File, delimiter=',')
            preAvg.append(Data[q,:-1])
            plt.subplot(3,2,q+1)
            for h in range(len(preAvg)):
                plt.plot(x_axis,preAvg[h])
                plt.plot(x_axis,preAvg[h],color=colours_round[h],linewidth=1)
                plt.xlim(x_min,xmax)
                plt.title(amps[q],pad=2,loc="left")
            plt.xlabel('Time (msec)')
            plt.ylabel(ylabel)
            plt.legend(rounds)
    figure.tight_layout(pad=.5)
    figname_AllInt = fname + Figname
    figpath_AllInt_path = root_path + figname_AllInt
    plt.savefig(figpath_AllInt_path)

def PlotAllParamAvg (AvgData,ylabel,fignumb,xmax,title):
    # Plot each parameter (avg) w/out sd, save as 1 png file
    plt.subplot(3,1,fignumb) # - plot avgs only for now
    for k in range(len(AvgData)):
        plt.plot(x_axis,AvgData[k],color=colours[k],linewidth=1)
        plt.xlim(x_min,xmax)
    plt.legend(amps,loc="upper left",bbox_to_anchor=(box_leg_x, box_leg_y),fontsize="x-small",ncol=2)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.title(title)
    plt.xlabel('Time (msec)')
    plt.ylabel(ylabel)

def PlotAvgPerInt(AvgData,NegSD,PosSD,ylabel,Figname,xmax):
    # Plot each avg int save as 1 png file for each parameter
    figure=plt.figure(figsize=(10,8))
    for n in range(len(AvgData)):
        plt.subplot(3,2,n+1)
        plt.plot(x_axis,AvgData[n],color=colours[n],linewidth=1)
        plt.plot(x_axis,NegSD[n],color=colours[n],linewidth=0.3)
        plt.plot(x_axis,PosSD[n],color=colours[n],linewidth=0.3)
        plt.xlim(x_min,xmax)
        plt.xlabel('Time (msec)')
        plt.ylabel(ylabel)
        plt.title(amps[n],pad=2,loc="left")
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
    figure.tight_layout(pad=.5)
    figname_AI = fname + Figname
    figpath_AI = root_path + figname_AI
    plt.savefig(figpath_AI)

def PlotAllParamAvgMaxPerInt(intens,avgMaxRd,NegSdMax,PosSdMax,fignumb,ylabel):
    # Plot max vs intensity for each parameter, save as 1 png file
    plt.subplot(3,1,fignumb)
    plt.scatter(intens,avgMaxRd,c=colours_round[fignumb-1])
    plt.scatter(intens,NegSdMax,s=11,c=colours_round[fignumb-1],marker="_")
    plt.scatter(intens,PosSdMax,s=11,c=colours_round[fignumb-1],marker="_")
    plt.xlabel('Intensity (mA)')
    plt.ylabel(ylabel)
    plt.ylim(bottom=0)
    
def PlotAllParamAvgTMaxPerInt(intens,avgTMaxRd,NegSdTMax,PosSdTMax,fignumb,ylabel):
    # Plot max vs intensity for each parameter, save as 1 png file
    plt.subplot(3,1,fignumb)
    plt.scatter(intens,avgTMaxRd,c=colours_round[fignumb-1])
    plt.scatter(intens,NegSdTMax,s=11,c=colours_round[fignumb-1],marker="_")
    plt.scatter(intens,PosSdTMax,s=11,c=colours_round[fignumb-1],marker="_")
    plt.xlabel('Intensity (mA)')
    plt.ylabel(ylabel)
    #plt.ylim(bottom=0)

### Algorithm #################################################################

# Specify data file root
folderListFile='C:\TestData\FolderListTest.txt'
data_path,folderNames=readFolderList(folderListFile)

# Loop through folders: avg data w/in fish, all int, plot w/out SD
for idF, folder in enumerate(folderNames):
    
    # Go into analysis folder, find cut files
    CumulAngCutFiles = glob.glob(folder+'/Analysis/*CumulAngCut.csv')
    CurvCutFiles = glob.glob(folder+'/Analysis/*CurvCut.csv')
    MotCutFiles = glob.glob(folder+'/Analysis/*MotCut.csv')
    
    # Average cut files from each round
    vAng, preAvgAng, avgRdAng, stdevRdAng = AverageRounds(CumulAngCutFiles,6)
    vCurv, preAvgCurv, avgRdCurv, stdevRdCurv = AverageRounds(CurvCutFiles,6)
    vMot, preAvgMot, avgRdMot, stdevRdMot = AverageRounds(MotCutFiles,6)
    
    # Get max (and time to max) from each round  and average them
    vmAng,maxRdAng,avgMaxRdAng,stdevMaxRdAng,TmaxRdAng,avgTMaxRdAng,stdevTMaxRdAng=MaxAvgTimeRounds(CumulAngCutFiles,6)
    vmCurv,maxRdCurv,avgMaxRdCurv,stdevMaxRdCurv,TmaxRdCurv,avgTMaxRdCurv,stdevTMaxRdCurv=MaxAvgTimeRounds(CurvCutFiles,6)
    vmMot,maxRdMot,avgMaxRdMot,stdevMaxRdMot,TmaxRdMot,avgTMaxRdMot,stdevTMaxRdMot=MaxAvgTimeRounds(MotCutFiles,6)
    
    
    ###### Create csv paths ######
    # Create root path and root fname and new summary folder
    root_path = folder + '/Analysis/Summary/'
    if not os.path.exists(root_path):
        os.mkdir(root_path)
    
    fname = folder[-18:] # works +/- w/ 2 digit dates & 1/2 digit fish #
    fname = fname.replace('_','')
    fname = fname.replace('\\','_')
    
    # Avg and matching SD per parameter
    avgRdAng_path = root_path + fname + '_avgRdAng.csv'
    avgRdCurv_path = root_path + fname + '_avgRdCurv.csv'
    avgRdMot_path = root_path + fname + '_avgRdMot.csv'
    stdevRdAng_path = root_path + fname + '_stdevRdAng.csv'
    stdevRdCurv_path = root_path + fname + '_stdevRdCurv.csv'
    stdevRdMot_path = root_path + fname + '_stdevRdMot.csv'
    
    # Avg max and matching SD per intensity per parameter
    avgMaxRdAng_path = root_path + fname + '_avgMaxRdAng.csv'
    avgMaxRdCurv_path = root_path + fname + '_avgMaxRdCurv.csv'
    avgMaxRdMot_path = root_path + fname + '_avgMaxRdMot.csv'
    stdevMaxRdAng_path = root_path + fname + '_stdevMaxRdAng.csv'
    stdevMaxRdCurv_path = root_path + fname + '_stdevMaxRdCurv.csv'
    stdevMaxRdMot_path = root_path + fname + '_stdevMaxRdMot.csv'
    
    # HAVEN'T DONE IT FOR TIME TO MAX
    
    ###### Save csv files ######
    # Avg and matching SD per parameter
    pd.DataFrame(avgRdAng).to_csv(avgRdAng_path, header=None, index=None)
    pd.DataFrame(avgRdCurv).to_csv(avgRdCurv_path, header=None, index=None)
    pd.DataFrame(avgRdMot).to_csv(avgRdMot_path, header=None, index=None)
    pd.DataFrame(stdevRdAng).to_csv(stdevRdAng_path, header=None, index=None)
    pd.DataFrame(stdevRdCurv).to_csv(stdevRdCurv_path, header=None, index=None)
    pd.DataFrame(stdevRdMot).to_csv(stdevRdMot_path, header=None, index=None)
    
    # Avg max and matching SD per intensity per parameter
    pd.DataFrame(avgMaxRdAng).to_csv(avgMaxRdAng_path, header=None, index=None)
    pd.DataFrame(avgMaxRdCurv).to_csv(avgMaxRdCurv_path, header=None, index=None)
    pd.DataFrame(avgMaxRdMot).to_csv(avgMaxRdMot_path, header=None, index=None)
    pd.DataFrame(stdevMaxRdAng).to_csv(stdevMaxRdAng_path, header=None, index=None)
    pd.DataFrame(stdevMaxRdCurv).to_csv(stdevMaxRdCurv_path, header=None, index=None)
    pd.DataFrame(stdevMaxRdMot).to_csv(stdevMaxRdMot_path, header=None, index=None)
    
    # HAVEN'T DONE IT FOR TIME TO MAX
    
    ###### Get SD traces ######
    # SD traces per parameter
    NegSdAng,PosSdAng = GetSdTraces(avgRdAng,stdevRdAng)
    NegSdCurv,PosSdCurv = GetSdTraces(avgRdCurv,stdevRdCurv)
    NegSdMot,PosSdMot = GetSdTraces(avgRdMot,stdevRdMot)
    
    # "Max" SD traces per intensity per parameter
    NegSdMaxAng,PosSdMaxAng = GetSdTraces(avgMaxRdAng,stdevMaxRdAng)
    NegSdMaxCurv,PosSdMaxCurv = GetSdTraces(avgMaxRdCurv,stdevMaxRdCurv)
    NegSdMaxMot,PosSdMaxMot = GetSdTraces(avgMaxRdMot,stdevMaxRdMot)
    
    # "Time to Max" SD traces per intensity per parameter
    NegSdTMaxAng,PosSdTMaxAng = GetSdTraces(avgTMaxRdAng,stdevTMaxRdAng)
    NegSdTMaxCurv,PosSdTMaxCurv = GetSdTraces(avgTMaxRdCurv,stdevTMaxRdCurv)
    NegSdTMaxMot,PosSdTMaxMot = GetSdTraces(avgTMaxRdMot,stdevTMaxRdMot)
    
    
    # Convert frames to msec
    num_frames = len(avgRdAng[0])
    frame_rate = 400 #fps
    num_msec = num_frames/frame_rate*1000
    x_ax = np.arange(0,num_msec,(num_msec/num_frames)) # Create new x values
    x_axis = x_ax - 1000 # because I had 1000 msec baseline
    
    
    ############# Plots ######################################################
    # Define plot parameters 
    box_leg_x = 1
    box_leg_y = 0.94
    x_min = -500
    colours = ["y","r","m","g","b","k"]
    amps = ['300mA','350mA','400mA','450mA','500mA','540mA']
    amps_noMA = ['300','350','400','450','500','540']
    rounds = ["R1","R2","R3"]
    colours_round =  ["b","orange","g"]
 
    #### Plotting AVG PER PARAMETER for ALL INTENSITIES same plot ####
    # Plot each parameter (avg) w/out sd, save as 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAllParamAvg(avgRdAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles Per Fish")
    PlotAllParamAvg(avgRdCurv,"Curvature (a.u.)",2,2500,"Average Curvatures Per Fish")
    PlotAllParamAvg(avgRdMot,"Motion (a.u.)",3,2500,"Average Motion Per Fish")
    figname_RdAv = fname + '_avgs_short.png'
    figpath_RdAv = root_path + figname_RdAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_RdAv)
    
    figure=plt.figure(figsize=(10,8))
    PlotAllParamAvg(avgRdAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles Per Fish")
    PlotAllParamAvg(avgRdCurv,"Curvature (a.u.)",2,10500,"Average Curvatures Per Fish")
    PlotAllParamAvg(avgRdMot,"Motion (a.u.)",3,10500,"Average Motion Per Fish")
    figname_RdAv = fname + '_avgs.png'
    figpath_RdAv = root_path + figname_RdAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_RdAv)
    
    
    #### Plotting each parameter with INTENSITIES SEPARATED OUT ####
    # Plot avg per intensity, save 1 png file per parameter
    PlotAvgPerInt(avgRdAng,NegSdAng,PosSdAng,"Cumul. Angles (deg)","AvgAngPerInt_short.png",2500)
    PlotAvgPerInt(avgRdAng,NegSdAng,PosSdAng,"Cumul. Angles (deg)","AvgAngPerInt.png",10500)
    PlotAvgPerInt(avgRdCurv,NegSdCurv,PosSdCurv,"Curvature (a.u.)","AvgCurvPerInt_short.png",2500)
    PlotAvgPerInt(avgRdCurv,NegSdCurv,PosSdCurv,"Curvature (a.u.)","AvgCurvPerInt.png",10500)
    PlotAvgPerInt(avgRdMot,NegSdMot,PosSdMot,"Motion (a.u.)","AvgMotPerInt_short.png",2500)
    PlotAvgPerInt(avgRdMot,NegSdMot,PosSdMot,"Motion (a.u.)","AvgMotPerInt.png",10500)
    
    # Plot individual rounds per intensity, save 1 png file per parameter
    PlotRoundsPerInt(CumulAngCutFiles,"Cumul. Ang. (deg)","_CumulAngPerRound_short.png",2500)
    PlotRoundsPerInt(CumulAngCutFiles,"Cumul. Ang. (deg)","_CumulAngPerRound.png",10500)
    PlotRoundsPerInt(CurvCutFiles,"Curvature (a.u.)","_CurvaturePerRound_short.png",2500)
    PlotRoundsPerInt(CurvCutFiles,"Curvature (a.u.)","_CurvaturePerRound.png",10500)
    PlotRoundsPerInt(MotCutFiles,"Motion (a.u.)","_MotionPerRound_short.png",2500)
    PlotRoundsPerInt(MotCutFiles,"Motion (a.u.)","_MotionPerRound.png",10500)
    
    
    #### Plotting PEAK values for each parameter vs intensity ####
    # Plot avg max vs intensity for each parameter, save 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAllParamAvgMaxPerInt(amps_noMA,avgMaxRdAng,NegSdMaxAng,PosSdMaxAng,1,"Peak Cumul Angle (deg)")
    PlotAllParamAvgMaxPerInt(amps_noMA,avgMaxRdCurv,NegSdMaxCurv,PosSdMaxCurv,2,"Peak Curvature (a.u.)")
    PlotAllParamAvgMaxPerInt(amps_noMA,avgMaxRdMot,NegSdMaxMot,PosSdMaxMot,3,"Peak Motion (a.u.)")
    figname_MaxRdAv = fname + '_max_per_intens.png'
    figpath_MaxRdAv = root_path + figname_MaxRdAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_MaxRdAv)
    
    # Plot "per round" max vs intensity for each parameter, save 1 png file
    

    #### Plotting TIME TO PEAK values for each parameter vs intensity ####
    # Plot avg time to max vs intensity for each parameter, save as 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAllParamAvgTMaxPerInt(amps_noMA,avgTMaxRdAng,NegSdTMaxAng,PosSdTMaxAng,1,"Time to Peak Cumul Angle (msec)")
    PlotAllParamAvgTMaxPerInt(amps_noMA,avgTMaxRdCurv,NegSdTMaxCurv,PosSdTMaxCurv,2,"Time to Peak Curvature (msec)")
    PlotAllParamAvgTMaxPerInt(amps_noMA,avgTMaxRdMot,NegSdTMaxMot,PosSdTMaxMot,3,"Time to Peak Motion (msec)")
    figname_TMaxRdAv = fname + '_time_max_per_intens.png'
    figpath_TMaxRdAv = root_path + figname_TMaxRdAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_TMaxRdAv)

    # Plot "per round" time to max vs intensity per parameter, save 1 png file
    
  
    
  
plt.close('all')

# FIN
  

    
    
    
    
    
    
    
    
    
    
    