# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:24:14 2020

Average data from csv files across fish, determine peak and plot

@author: Elisa C & Tom Ryan
"""

import os
#import cv2
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import glob
from scipy import stats
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
  
       
def AverageFish(Files,num_stim):
    volts=[]
    avgAllFish=[]
    stdevAllFish=[]
    semAllFish=[]
   
    for i in range(num_stim): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files): # loop through files (each fish)
            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
            preAvg.append(Data[i,:])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
   
        meanOfTrials=np.mean(preAvg,0)
        avgAllFish.append(meanOfTrials)   
        stdevAllFish.append(np.std(preAvg,0))
        semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return volts, preAvg, avgAllFish, stdevAllFish, semAllFish


def MaxAvgTimeFish(Files,num_stim):   
    volts = []
    maxAF = []
    avgMaxAF = []
    #stdevMaxAF = []
    TMaxAF = []
    avgTMaxAF = []
    #stdevTMaxAF = []
    semMaxAF = []
    semTMaxAF = []
   
    for i in range(num_stim): # loop through stim
        temp=[]
        for x,File in enumerate (Files): # loop through files (each round)
            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
            temp.append(Data[i,:-1])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
            maxInt = np.amax(temp)
            # timeMax = np.argmax(temp)
            frameMax = np.argmax(temp)
            timeMax = frameMax/400*1000
            maxAF.append(maxInt)
            
            TMaxAF.append(timeMax)
        print(maxAF)   
        meanOfMaxTrials=np.mean(maxAF,0)
        avgMaxAF.append(meanOfMaxTrials)
        #stdevMaxAF.append(np.std(maxAF,0))
        semMaxAF.append(stats.sem(maxAF,axis=0,ddof=0))
        meanOfTimeMaxTrials=np.mean(TMaxAF,0)
        avgTMaxAF.append(meanOfTimeMaxTrials)
        #stdevTMaxAF.append(np.std(TMaxAF,0))
        semTMaxAF.append(stats.sem(TMaxAF,axis=0,ddof=0))
        volts.append(Data[i,-1])
       
    return volts, maxAF, avgMaxAF, semMaxAF, TMaxAF, avgTMaxAF, semTMaxAF


def GetSemTraces(avData,semData):
    all_neg_SEM = []
    all_pos_SEM = []
    
    for j in range(len(avData)):
        neg_SEM = np.subtract(avData[j],semData[j])
        pos_SEM = np.add(avData[j],semData[j])
        all_neg_SEM.append(neg_SEM)
        all_pos_SEM.append(pos_SEM)
    
    return all_neg_SEM, all_pos_SEM


######## Plotting functions #########
def PlotFishPerInt(Files,ylabel,Figname,xmax):
    # Plot individual fish per intensity
    figure=plt.figure(figsize=(10,8))
    for q in range(6): # loop through stim
        preAvg=[]
        for x,File in enumerate (Files):
            Data = np.genfromtxt(File, delimiter=',')
            preAvg.append(Data[q,:])
            plt.subplot(3,2,q+1)
            for h in range(len(preAvg)):
                plt.plot(x_axis,preAvg[h],color=colours_round[h])
                plt.xlim(x_min,xmax)
            plt.title(amps[q],pad=2,loc="left")
            plt.xlabel('Time (msec)')
            plt.ylabel(ylabel)
            plt.legend(fish)
    figure.tight_layout(pad=.5)
    #figname_AllInt = fname + Figname
    AllFish_AllInt_path = analysis_path + Figname
    plt.savefig(AllFish_AllInt_path)


def PlotAllFishAllParamAvg (AvgData,ylabel,fignumb,xmax,title):
    # Plot each parameter (avg) w/out sd, save as 1 png file
    plt.subplot(3,1,fignumb) # - plot avgs only for now
    for k in range(len(AvgData)):
        plt.plot(x_axis,AvgData[k],color=colours[k])
        plt.xlim(x_min,xmax)
    plt.legend(amps,loc="upper left",bbox_to_anchor=(box_leg_x, box_leg_y),fontsize="x-small",ncol=2)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.title(title)
    plt.xlabel('Time (msec)')
    plt.ylabel(ylabel)

def PlotAllFishAvgPerInt(AvgData,NegSEM,PosSEM,ylabel,Figname,xmax):
    # Plot each avg int save as 1 png file for each parameter
    figure=plt.figure(figsize=(10,8))
    for n in range(len(AvgData)):
        plt.subplot(3,2,n+1)
        plt.plot(x_axis,AvgData[n],color=colours[n])
        plt.plot(x_axis,NegSEM[n],color=colours[n],linewidth=0.3)
        plt.plot(x_axis,PosSEM[n],color=colours[n],linewidth=0.3)
        plt.xlim(x_min,xmax)
        plt.xlabel('Time (msec)')
        plt.ylabel(ylabel)
        plt.title(amps[n],pad=2,loc="left")
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
    figure.tight_layout(pad=.5)
    figpath_AFAI = analysis_path + Figname
    plt.savefig(figpath_AFAI)

def PlotAFAllParamAvgMaxPerInt(intens,avgMaxAF,NegSdMaxAF,PosSdMaxAF,fignumb,ylabel):
    # Plot max vs intensity for each parameter, save as 1 png file
    plt.subplot(3,1,fignumb)
    plt.scatter(intens,avgMaxAF,c=colours_round[fignumb-1])
    plt.scatter(intens,NegSdMaxAF,s=11,c=colours_round[fignumb-1],marker="_")
    plt.scatter(intens,PosSdMaxAF,s=11,c=colours_round[fignumb-1],marker="_")
    plt.xlabel('Intensity (mA)')
    plt.ylabel(ylabel)
    plt.ylim(bottom=0)
    
def PlotAFAllParamAvgTMaxPerInt(intens,avgTMaxAF,NegSdTMaxAF,PosSdTMaxAF,fignumb,ylabel):
    # Plot max vs intensity for each parameter, save as 1 png file
    plt.subplot(3,1,fignumb)
    plt.scatter(intens,avgTMaxAF,c=colours_round[fignumb-1])
    plt.scatter(intens,NegSdTMaxAF,s=11,c=colours_round[fignumb-1],marker="_")
    plt.scatter(intens,PosSdTMaxAF,s=11,c=colours_round[fignumb-1],marker="_")
    plt.xlabel('Intensity (mA)')
    plt.ylabel(ylabel)
    #plt.ylim(bottom=0)
     
### Algorithm #################################################################

# Specify data file root
folderListFile='C:\TestData\FolderListTest.txt'
data_path,folderNames=readFolderList(folderListFile)

AvgCumulAngFiles = []
AvgCurvFiles = []
AvgMotFiles = []

# Loop through folders: avg data w/in fish, all int, plot w/out SD
for idF, folder in enumerate(folderNames):
    
    # Go into analysis folder, find cut files
    AvgCumulAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*avgRdAng.csv'))
    AvgCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*avgRdCurv.csv'))
    AvgMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*avgRdMot.csv'))
    
    # Average avg files from each fish
    vAng, preAvgAng, avgAllAng, stdevAllAng, semAllAng = AverageFish(AvgCumulAngFiles,6)
    vCurv, preAvgCurv, avgAllCurv, stdevAllCurv, semAllCurv = AverageFish(AvgCurvFiles,6)
    vMot, preAvgMot, avgAllMot, stdevAllMot, semAllMot = AverageFish(AvgMotFiles,6)
    
    # Get max (and time to max) from each fish and average them
    vmAng,maxAFAng,avgMaxAFAng,semMaxAFAng,TMaxAFAng,avgTMaxAFAng,semTMaxAFAng=MaxAvgTimeFish(AvgCumulAngFiles,6)
    vmCurv,maxAFCurv,avgMaxAFCurv,semMaxAFCurv,TMaxAFCurv,avgTMaxAFCurv,semTMaxAFCurv=MaxAvgTimeFish(AvgCurvFiles,6)
    vmMot,maxAFMot,avgMaxAFMot,semMaxAFMot,TMaxAFMot,avgTMaxAFMot,semTMaxAFMot=MaxAvgTimeFish(AvgMotFiles,6)
    
    
    ###### Create csv paths ######
    # Create root path and root fname and new summary folder
    folderFile = open(folderListFile, 'r')
    folderList = folderFile.readlines()
    analysis_path = folderList[0][:-1] + "/Analysis_Fish/"
    if not os.path.exists(analysis_path):
        os.mkdir(analysis_path)
    
    # Avg and matching SEM per parameter
    avgAllAng_path = analysis_path + 'avgAllAng.csv'
    avgAllCurv_path = analysis_path + 'avgAllCurv.csv'
    avgAllMot_path = analysis_path + 'avgAllMot.csv'
    #stdevAllAng_path = analysis_path + 'stdevAllAng.csv'
    #stdevAllCurv_path = analysis_path + 'stdevAllCurv.csv'
    #stdevAllMot_path = analysis_path + 'stdevAllMot.csv'
    semAllAng_path = analysis_path + 'semAllAng.csv'
    semAllCurv_path = analysis_path + 'semAllCurv.csv'
    semAllMot_path = analysis_path + 'semAllMot.csv'
    
    # Avg max and matching SD per intensity per parameter
    avgMaxAFAng_path = analysis_path + '_avgMaxAFAng.csv'
    avgMaxAFCurv_path = analysis_path + '_avgMaxAFCurv.csv'
    avgMaxAFMot_path = analysis_path + '_avgMaxAFMot.csv'
    semMaxAFAng_path = analysis_path + '_semMaxAFAng.csv'
    semMaxAFCurv_path = analysis_path + '_semMaxAFCurv.csv'
    semMaxAFMot_path = analysis_path + '_ssemMaxAFMot.csv'
    
    # HAVEN'T DONE IT FOR TIME TO MAX
    
    
    ###### Save csv files ######
    # Avg and matching SEM per parameter
    pd.DataFrame(avgAllAng).to_csv(avgAllAng_path, header=None, index=None)
    pd.DataFrame(avgAllCurv).to_csv(avgAllCurv_path, header=None, index=None)
    pd.DataFrame(avgAllMot).to_csv(avgAllMot_path, header=None, index=None)
    #pd.DataFrame(stdevAllAng).to_csv(stdevAllAng_path, header=None, index=None)
    #pd.DataFrame(stdevAllCurv).to_csv(stdevAllCurv_path, header=None, index=None)
    #pd.DataFrame(stdevAllMot).to_csv(stdevAllMot_path, header=None, index=None)
    pd.DataFrame(semAllAng).to_csv(semAllAng_path, header=None, index=None)
    pd.DataFrame(semAllCurv).to_csv(semAllCurv_path, header=None, index=None)
    pd.DataFrame(semAllMot).to_csv(semAllMot_path, header=None, index=None)
    
    # Avg max and matching SD per intensity per parameter
    pd.DataFrame(avgMaxAFAng).to_csv(avgMaxAFAng_path, header=None, index=None)
    pd.DataFrame(avgMaxAFCurv).to_csv(avgMaxAFCurv_path, header=None, index=None)
    pd.DataFrame(avgMaxAFMot).to_csv(avgMaxAFMot_path, header=None, index=None)
    pd.DataFrame(semMaxAFAng).to_csv(semMaxAFAng_path, header=None, index=None)
    pd.DataFrame(semMaxAFCurv).to_csv(semMaxAFCurv_path, header=None, index=None)
    pd.DataFrame(semMaxAFMot).to_csv(semMaxAFMot_path, header=None, index=None)
    
    # HAVEN'T DONE IT FOR TIME TO MAX
    
    
    ###### Get SEM traces ######
    # SEM traces per parameter
    NegSemAng,PosSemAng = GetSemTraces(avgAllAng,semAllAng)
    NegSemCurv,PosSemCurv = GetSemTraces(avgAllCurv,semAllCurv)
    NegSemMot,PosSemMot = GetSemTraces(avgAllMot,semAllMot)
    
    # "Max" SEM traces per intensity per parameter
    NegSemMaxAng,PosSemMaxAng = GetSemTraces(avgMaxAFAng,semMaxAFAng)
    NegSemMaxCurv,PosSemMaxCurv = GetSemTraces(avgMaxAFCurv,semMaxAFCurv)
    NegSemMaxMot,PosSemMaxMot = GetSemTraces(avgMaxAFMot,semMaxAFMot)
    
    # "Time to Max" SEM traces per intensity per parameter
    NegSemTMaxAng,PosSemTMaxAng = GetSemTraces(avgTMaxAFAng,semTMaxAFAng)
    NegSemTMaxCurv,PosSemTMaxCurv = GetSemTraces(avgTMaxAFCurv,semTMaxAFCurv)
    NegSemTMaxMot,PosSemTMaxMot = GetSemTraces(avgTMaxAFMot,semTMaxAFMot)
    
    
    # Convert frames to msec
    num_frames = len(avgAllAng[0])
    frame_rate = 400 #fps
    num_msec = num_frames/frame_rate*1000
    x_ax = np.arange(0,num_msec,(num_msec/num_frames)) # Create new x values
    x_axis = x_ax - 1000 # because I had 1000 msec baseline
    
    
    ############# Plots ######################################################
    # Define plot parameters 
    box_leg_x = 1
    box_leg_y = 0.94
    x_min = -500
    #x_max = 2500
    colours = ["y","r","m","g","b","k"]
    amps = ['300mA','350mA','400mA','450mA','500mA','540mA']
    amps_noMA = ['300','350','400','450','500','540']
    fish = ["F1","F2","F3"]
    #colours_round = ["brown","orange","blue"]
    #colours_round =  ["y","r","g"]
    colours_round =  ["b","orange","g"]
    
    #### Plotting AVG PER PARAMETER for ALL INTENSITIES same plot ####
    # Plot each parameter (avg) w/out sd, save 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish")
    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (a.u.)",2,2500,"Average Curvatures For All Fish")
    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish")
    figname_AllAvSh = 'fish_avgs_short.png'
    figpath_AllAvSh = analysis_path + figname_AllAvSh
    figure.tight_layout(pad=1)
    plt.savefig(figpath_AllAvSh)
    
    figure=plt.figure(figsize=(10,8))
    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish")
    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (a.u.)",2,10500,"Average Curvatures For All Fish")
    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish")
    figname_AllAv = 'fish_avgs.png'
    figpath_AllAv = analysis_path + figname_AllAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_AllAv)
    
    
    #### Plotting each parameter with INTENSITIES SEPARATED OUT ####
    # Plot each avg int save 1 png file per parameter
    PlotAllFishAvgPerInt(avgAllAng,NegSemAng,PosSemAng,"Cumul. Angles (deg)","AllFish_AvgAngPerInt_short.png",2500)
    PlotAllFishAvgPerInt(avgAllAng,NegSemAng,PosSemAng,"Cumul. Angles (deg)","AllFish_AvgAngPerInt.png",10500)
    PlotAllFishAvgPerInt(avgAllCurv,NegSemCurv,PosSemCurv,"Curvature (a.u.)","AllFish_AvgCurvPerInt_short.png",2500)
    PlotAllFishAvgPerInt(avgAllCurv,NegSemCurv,PosSemCurv,"Curvature (a.u.)","AllFish_AvgCurvPerInt.png",10500)
    PlotAllFishAvgPerInt(avgAllMot,NegSemMot,PosSemMot,"Motion (a.u.)","AllFish_AvgMotPerInt_short.png",2500)
    PlotAllFishAvgPerInt(avgAllMot,NegSemMot,PosSemMot,"Motion (a.u.)","AllFish_AvgMotPerInt.png",10500)
    
    # Plot individual fish per intensity, save 1 png file per parameter
    PlotFishPerInt(AvgCumulAngFiles,"Cumul. Ang. (deg)","CumulAngPerFish_short.png",2500)
    PlotFishPerInt(AvgCumulAngFiles,"Cumul. Ang. (deg)","CumulAngPerFish.png",10500)
    PlotFishPerInt(AvgCurvFiles,"Curvature","CurvaturePerFish_short.png",2500)
    PlotFishPerInt(AvgCurvFiles,"Curvature","CurvaturePerFish.png",10500)
    PlotFishPerInt(AvgMotFiles,"Motion","MotionPerFish_short.png",2500)
    PlotFishPerInt(AvgMotFiles,"Motion","MotionPerFish.png",10500)
 
    
    #### Plotting PEAK values for each parameter vs intensity ####
    # Plot avg max vs intensity for each parameter, save 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAFAllParamAvgMaxPerInt(amps_noMA,avgMaxAFAng,NegSemMaxAng,PosSemMaxAng,1,"Peak Cumul Angle (deg)")
    PlotAFAllParamAvgMaxPerInt(amps_noMA,avgMaxAFCurv,NegSemMaxCurv,PosSemMaxCurv,2,"Peak Curvature (a.u.)")
    PlotAFAllParamAvgMaxPerInt(amps_noMA,avgMaxAFMot,NegSemMaxMot,PosSemMaxMot,3,"Peak Motion (a.u.)")
    figname_MaxAFAv = '_max_per_intens.png'
    figpath_MaxAFAv = analysis_path + figname_MaxAFAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_MaxAFAv)
    
    # Plot "per fish" max vs intensity for each parameter, save 1 png file
    

    #### Plotting TIME TO PEAK values for each parameter vs intensity ####
    # Plot avg time to max vs intensity for each parameter, save as 1 png file
    figure=plt.figure(figsize=(10,8))
    PlotAFAllParamAvgTMaxPerInt(amps_noMA,avgTMaxAFAng,NegSemTMaxAng,PosSemTMaxAng,1,"Time to Peak Cumul Angle (msec)")
    PlotAFAllParamAvgTMaxPerInt(amps_noMA,avgTMaxAFCurv,NegSemTMaxCurv,PosSemTMaxCurv,2,"Time to Peak Curvature (msec)")
    PlotAFAllParamAvgTMaxPerInt(amps_noMA,avgTMaxAFMot,NegSemTMaxMot,PosSemTMaxMot,3,"Time to Peak Motion (msec)")
    figname_TMaxAFAv = '_time_max_per_intens.png'
    figpath_TMaxAFAv = analysis_path + figname_TMaxAFAv
    figure.tight_layout(pad=1)
    plt.savefig(figpath_TMaxAFAv)

    # Plot "per fish" time to max vs intensity per parameter, save 1 png file
    
  
    
plt.close('all')











