# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:24:14 2020

Average data from csv files across fish, determine peak and plot

@author: Elisa C & Tom Ryan
"""
############ Set default inputs and options DO NOT CHANGE ##########
folderListFile=[]
Dose=''
#####################################################
# Set "Library Path" - Nociceptive Zebrafish Repo
lib_path='C:/Users/thoma/OneDrive/Documents/GitHub/Nociceptive_Zebrafish/libs'
import sys
sys.path.append(lib_path)

# Generic Libraries
import glob
import numpy as np
import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
# Custom libraries
import NZ_rounds as NZR
import NZ_utilities as NZU
import NZ_dictionary as NZD
import NZ_figures as NZF
################################## INPUTS #####################################
############ Directory info and options ##############
root='D:/Elisa/'
folderListFileN='PairedLido5mg.txt'
folderListFile=root+folderListFileN
outFolder='/Analysis/Summary/'
paired=True
grabDict=True # if False, will search for relevant csvs instead
saveCSV=False
group=True # group rounds together with averages into a single dictionary. Will delete the single round dictionaries created in step 2
saveDict=True
deleteIndDicts=False
ServerExchange = True 
plot=False
numStim=6
set_dpi=600

######################################################
if grabDict: 
    print('grabbing dictionaries')
else: print('grabbing csvs')
if saveCSV and saveDict: print('saving as both csv files and dictionaries')
if saveCSV and not saveDict: print('saving as csv files')
if saveDict and not saveCSV: print('saving as dictionaries')

data_path,folderNames=NZU.readFolderList(folderListFile)

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
            preAvg.append(Data[i,:])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
   
        meanOfTrials=np.mean(preAvg,0)
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
            try:                
                preAvg.append(Data[i])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
            except:
                print('debug')

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

    meanOfTrials=np.mean(preAvg,0)
    avgAllFish.append(meanOfTrials)   
    stdevAllFish.append(np.std(preAvg,0))
    semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
       
    return preAvg, avgAllFish, stdevAllFish, semAllFish


# start of fish Loop ##
DictList=[]
for idF, folder in enumerate(folderNames):
    if ServerExchange: 
        folder_serv=folder
        fpath=folder.split(sep='\\',maxsplit=5)[-1]
        folder=root+'Data/'+fpath
    # load average dictionaries per fish
    DictFolder=folder+'\\Dictionaries\\'
    ff=glob.glob(DictFolder+'*_AverageDict.npy')
    DictList.append(ff[0])
    Dict=np.load(ff[0],allow_pickle=True).item()
    
    # for each key, normalise preDrug to 1, and normalise postDrug by preDrug
    preDrug=Dict['preDrug']
    postDrug=Dict['postDrug']
    Keys=list(preDrug.keys())[3:-3]
    
    diffDict={}
    for i,thisKey in enumerate(Keys):
        thisPreDrugParam=preDrug[thisKey]['Mean']
        thisPostDrugParam=postDrug[thisKey]['Mean']
        
        percent_difference=np.divide(thisPostDrugParam,thisPreDrugParam)
        diffDict.update({thisKey:percent_difference})
    Dict.update({'PercentDiff':diffDict})
    if saveDict:
        np.save(ff[0],Dict)

# type is 0 for rounds, 1 for single and 2 for parameters
typeOfMetric=[0,0,0,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
Keys=list(preDrug.keys())[0:-3]
preDrugAvgS=[]
preDrugSemS=[]
postDrugAvgS=[]
postDrugSemS=[]
diffAvgS=[]
diffSemS=[]
for i,thisKey in enumerate(Keys):
    diffavg=[]
    diffsem=[]
    if typeOfMetric[i]==0:
        _, preavg, _, presem = AverageFishDict(DictList,thisKey,6,preDrug=True,diff=False)
        _, postavg, _, postsem = AverageFishDict(DictList,thisKey,6,preDrug=False,diff=False)
#        _, diffavg, _, diffsem = AverageFishDict(DictList,thisKey,6,preDrug=False,diff=True)
    elif typeOfMetric[i]==1:
        _, preavg, _, presem = AverageFishSingleDict(DictList,thisKey,preDrug=True,diff=False)
        _, postavg, _, postsem = AverageFishSingleDict(DictList,thisKey,preDrug=False,diff=False)
        _, diffavg, _, diffsem = AverageFishSingleDict(DictList,thisKey,preDrug=False,diff=True)
    elif typeOfMetric[i]==2:
        _, preavg, _, presem = AverageFishParamDict(DictList,thisKey,6,preDrug=True,diff=False)
        _, postavg, _, postsem = AverageFishParamDict(DictList,thisKey,6,preDrug=False,diff=False)
        _, diffavg, _, diffsem = AverageFishParamDict(DictList,thisKey,6,preDrug=False,diff=True)
    preDrugAvgS.append(preavg)
    preDrugSemS.append(presem)
    postDrugAvgS.append(postavg)
    postDrugSemS.append(postsem)
    diffAvgS.append(diffavg)
    diffSemS.append(diffsem)
    diffavg=[]
    diffsem=[]

## FIN
#import os
##import cv2
#import numpy as np
#import pandas as pd 
#import matplotlib.pyplot as plt
#from matplotlib import cm
#import glob
#from scipy import stats
#from statsmodels.stats.multicomp import pairwise_tukeyhsd
#from statsmodels.stats.multicomp import MultiComparison
#from scikit_posthocs import posthoc_dunn
#import timeit
#from matplotlib.lines import Line2D
##import csv as csv
#
#
#### Helper Functions
#
#def readFolderList(folderListFile):
#    folderFile=open(folderListFile, 'r')
#    folderList=folderFile.readlines()
#    data_path=folderList[0][:-1]
#    folderList=folderList[1:]
#    folderNames=[]
#    for i,f in enumerate(folderList):
#        stringline = f[:].split()
#        expFolderName=data_path + stringline[0]
#        folderNames.append(expFolderName) # create a list of folders from folderlistfile
#        
#    return data_path,folderNames
#  
#       
#def AverageFish(Files,num_stim):
#    # Standard for raw data (mot, curv, cumulAng)
#    avgAllFish=[]
#    stdevAllFish=[]
#    semAllFish=[]
#   
#    for i in range(num_stim): # loop through stim
#        preAvg=[]
#        for x,File in enumerate (Files): # loop through files (each fish)
#            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
#            if len(Data[i]) == 4601:
#                preAvg.append(Data[i][:-1])
#            elif len(Data[i]) == 4600:
#                preAvg.append(Data[i,:])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
#   
#        meanOfTrials=np.mean(preAvg,0)
#        avgAllFish.append(meanOfTrials)   
#        stdevAllFish.append(np.std(preAvg,0))
#        semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
#       
#    return preAvg, avgAllFish, stdevAllFish, semAllFish
#
#
#def AverageFishParam(Files,num_stim):   
#    # For processed data e.g. latency
#    avgAllFish=[]
#    stdevAllFish=[]
#    semAllFish=[]
#   
#    for i in range(num_stim): # loop through stim
#        preAvg=[]
#        for x,File in enumerate (Files): # loop through files (each fish)
#            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
#            preAvg.append(Data[i])  # populate list with ith row from each data file (not very efficient as opens each file many many times)
#
#        meanOfTrials=np.nanmean(preAvg,0)
#        avgAllFish.append(meanOfTrials)   
#        stdevAllFish.append(np.nanstd(preAvg,0))
#        semAllFish.append(stats.sem(preAvg,axis=0,ddof=0,nan_policy='omit'))
#       
#    return preAvg, avgAllFish, stdevAllFish, semAllFish
#
#
#def AverageFishSingle(Files):   
#    # For fish freq resp (append "doesnt like" excel file with single row?)
#    avgAllFish=[]
#    stdevAllFish=[]
#    semAllFish=[]
#    preAvg=[]
#    
#    for x,File in enumerate (Files): # loop through files (each fish)
#        Data_1 = np.genfromtxt(File, delimiter=',') # grab data from the file
#        preAvg.append(Data_1)  # populate list with ith row from each data file (not very efficient as opens each file many many times)
#
#    meanOfTrials=np.mean(preAvg,0)
#    avgAllFish.append(meanOfTrials)   
#    stdevAllFish.append(np.std(preAvg,0))
#    semAllFish.append(stats.sem(preAvg,axis=0,ddof=0))
#       
#    return preAvg, avgAllFish, stdevAllFish, semAllFish
#
#
#def GetSemTraces(avData,semData):
#    all_neg_SEM = []
#    all_pos_SEM = []
#    
#    for j in range(len(avData)):
#        neg_SEM = np.subtract(avData[j],semData[j])
#        pos_SEM = np.add(avData[j],semData[j])
#        all_neg_SEM.append(neg_SEM)
#        all_pos_SEM.append(pos_SEM)
#    
#    return all_neg_SEM, all_pos_SEM
#
#
#def MyStats(Files):
#    
#    AN=[0 for i in range(6)]
#    for i in range(6): # loop through stim
#        nonan=[]
#        for x,File in enumerate (Files): # loop through files (each fish)
#            Data = np.genfromtxt(File, delimiter=',') # grab data from the file
#            y = Data[i][np.logical_not(np.isnan(Data[i]))]
#            nonan.append(y)
#        AN[i]=np.concatenate(nonan)
#    if len(np.concatenate(AN)) >=3:
#        normal = stats.shapiro(np.concatenate(AN))
#        if normal[1] > 0.05:
#            final_merged=[]
#            for i in range(6):
#                full=np.full(len(AN[i]),amps_noMA[i])
#                AN_int=AN[i]
#                merged_list = [(full[j], AN_int[j]) for j in range(0, len(full))]
#                final_merged.extend(merged_list)
#            
#            data = np.rec.array(final_merged, dtype = [('Laser driver current','|U5'),('Data', '<i8')])
#            
#            f, p = stats.f_oneway(data[data['Laser driver current'] == '300'].Data,
#                                  data[data['Laser driver current'] == '350'].Data,
#                                  data[data['Laser driver current'] == '400'].Data,
#            	                  data[data['Laser driver current'] == '450'].Data,
#                                  data[data['Laser driver current'] == '500'].Data,
#                                  data[data['Laser driver current'] == '540'].Data)
#            mc = MultiComparison(data['Data'], data['Laser driver current'])
#            result = mc.tukeyhsd()
#            df = pd.DataFrame(data=result._results_table.data[1:], columns=result._results_table.data[0])
#            fp_arr=[f,p]
#        else:
#            f,p=stats.kruskal(AN[0],AN[1],AN[2],AN[3],AN[4],AN[5])
#            dunn_list=[AN[0],AN[1],AN[2],AN[3],AN[4],AN[5]]
#            df=posthoc_dunn(dunn_list)
#            fp_arr=[f,p]
#    else:
#        fp_arr = []
#        df = []
#    
#    return fp_arr, df
#      
#
######### Plotting functions #########
#def PlotAllFishAllParamAvg (exptAvgData,ylabel,fignumb,xmax,title,AvgData):
#    # Plot each parameter (avg) w/out sd, save as 1 tiff file
#    plots = plt.subplot(3,1,fignumb)
#    plots.spines['right'].set_visible(False)
#    plots.spines['top'].set_visible(False)
#    tops_1 = []
#    bottoms_1 = []
#    tops_2 = []
#    bottoms_2 = []    
#
#    for k in range(len(AvgData)):
#        temp,=plt.plot(x_axis,AvgData[k]) 
#        temp.remove()
#        bottom_1,top_1 = plt.ylim()
#        tops_1.append(top_1)
#        bottoms_1.append(bottom_1)
#        # suboptimal but couldn't get it to work otherwise: simply finding min
#        # in AvgData[k] yielded slitghly =/= results than getting ylim of plt 
#    ymax_1 = np.amax(tops_1)
#    ymin_1 = np.amin(bottoms_1)
#     
#    for k in range(len(exptAvgData)):
#        plt.plot(x_axis,exptAvgData[k],color=cmap(rgb[k]),linewidth=1)
#        plt.xlim(x_min,xmax)
#        bottom_2,top_2 = plt.ylim()
#        tops_2.append(top_2)
#        bottoms_2.append(bottom_2)
#    ymax_2 = np.amax(tops_2)
#    ymin_2 = np.amin(bottoms_2)  
#
#    ymax = np.amax((ymax_1,ymax_2))
#    ymin = np.amin((ymin_1,ymin_2))    
#    
#    plt.legend(amps,loc="upper left",bbox_to_anchor=(box_leg_x, box_leg_y),fontsize="x-small")
#    plt.xticks(fontsize=7)
#    plt.yticks(fontsize=7)
#    plt.title(title)
#    plt.ylim(ymin,ymax)
#    plt.xlabel('Time (ms)')
#    plt.ylabel(ylabel)
#
#        
#def PlotAllFishAllParamAvgSm (exptAvgData,ylabel,fignumb,xmax,title,AvgData):
#    # Plot each parameter (avg) w/out sd, save as 1 tiff file
#    plots = plt.subplot(3,1,fignumb)
#    plots.spines['right'].set_visible(False)
#    plots.spines['top'].set_visible(False)
#    tops_1 = []
#    bottoms_1 = []
#    tops_2 = []
#    bottoms_2 = []    
#    sizeSm = 5
#    xxx = -(sizeSm-1)
#
#    for k in range(len(AvgData)):
#        temp,= plt.plot(x_axis,AvgData[k])
#        temp.remove()
#        bottom_1,top_1 = plt.ylim()
#        tops_1.append(top_1)
#        bottoms_1.append(bottom_1)
#    ymax_1 = np.amax(tops_1)
#    ymin_1 = np.amin(bottoms_1) 
#    
#    for k in range(len(exptAvgData)):
#        smooth_data=np.convolve(exptAvgData[k], np.ones((int(sizeSm),))/int(sizeSm), mode='valid')
#        plt.plot(x_axis[0:xxx],smooth_data,color=cmap(rgb[k]),linewidth=1)
#        plt.xlim(x_min,xmax)
#        bottom_2,top_2 = plt.ylim()
#        tops_2.append(top_2)
#        bottoms_2.append(bottom_2)
#    ymax_2 = np.amax(tops_2)
#    ymin_2 = np.amin(bottoms_2)  
#
#    ymax = np.amax((ymax_1,ymax_2))
#    ymin = np.amin((ymin_1,ymin_2))
#         
#    plt.legend(amps,loc="upper left",bbox_to_anchor=(box_leg_x, box_leg_y),fontsize="x-small")
#    plt.xticks(fontsize=7)
#    plt.yticks(fontsize=7)
#    plt.title(title)
#    plt.ylim(ymin,ymax)
#    plt.xlabel('Time (ms)')
#    plt.ylabel(ylabel)
#    
#
#### MAX ###
#def PlotAFAllParamAvgMaxPerInt(Files,intens,avgMaxAF,SdMax,EXPTFiles,EXPTavgMaxAF,EXPTSdMax,fignumb,ylabel):
#    # Plot max+TMax vs Laser driver current for each parameter, save as 1 tiff file
#    if exptType == 'AITC':
#        myCol = 'r'
#        experimentLabel = 'AITC' 
#    if exptType == 'lidocaine':
#        myCol = 'b'
#        experimentLabel = 'Lidocaine' 
#
#    legend_elements = [Line2D([0], [0], color=myCol, lw=2, label=experimentLabel),Line2D([0], [0], color='k', lw=2, label='Control')]
#        
#    plots = plt.subplot(3,2,fignumb)
#    plots.spines['right'].set_visible(False)
#    plots.spines['top'].set_visible(False)
#    for p in range(len(avgMaxAF)):
#        plots.errorbar(intens[p],avgMaxAF[p],SdMax[p],None,c='k',linestyle='None',marker='.',markersize=9,capsize=3)
#        plots.errorbar(intens[p],EXPTavgMaxAF[p],EXPTSdMax[p],None,c=myCol,linestyle='None',marker='^',markersize=5,capsize=3)
#    plots.legend(handles=legend_elements)
#    #plt.legend(legd_expt,loc= 'upper left',frameon=False,fontsize="x-small")
#    # handles, labels = plots.get_legend_handles_labels()    
#    # handles = [h[0] for h in handles] # remove errorbars
#    # plots.legend(handles, labels) # use them in the legend
#    #plots.legendHandles.set_color('black')
#    plt.xlabel('Laser driver current (mA)')
#    plt.ylabel(ylabel)
#    bottom, top = plt.ylim()
#    new_bottom = min(bottom,0)
#    new_top = top + 20
#    plt.ylim(new_bottom,new_top)
#
#### OTHER PARAMETERS ###
#def PlotAFAvgLatFreqBoutVelo(intens,avgRd,SD,EXPTavgRd,EXPTSD,fignumb,ylabel):
#    # Plot avg latency, resp freq, bout #, velocity vs Laser driver current
#    if exptType == 'AITC':
#        myCol = 'r'
#        experimentLabel = 'AITC' 
#    if exptType == 'lidocaine':
#        myCol = 'b'
#        experimentLabel = 'Lidocaine' 
#    
#    legend_elements = [Line2D([0], [0], color=myCol, lw=2, label=experimentLabel),Line2D([0], [0], color='k', lw=2, label='Control')]
#    plots = plt.subplot(3,2,fignumb)
#    plots.spines['right'].set_visible(False)
#    plots.spines['top'].set_visible(False)
#    for p in range(len(avgRd)):
#        plt.errorbar(intens[p],avgRd[p],SD[p],None,c='k',linestyle='None',marker='.',capsize=3)
#        plt.errorbar(intens[p],EXPTavgRd[p],EXPTSD[p],None,c=myCol,linestyle='None',marker='.',capsize=3)
#    plots.legend(handles=legend_elements)
#    plt.xlabel('Laser driver current (mA)')
#    plt.ylabel(ylabel)
#    if avgRd == avgAllFreqPI:
#        bottomm, topp=plots.get_ylim()
#        ymax = 1
#        ymin = bottomm
#        plt.ylim(ymin,ymax)
#        
#    
#
#### Algorithm #################################################################
#
## NB: If I add to S4a the option to save SEM as csv files, I can change S4aCombo
## such that it just imports AVG & SEM rather than calculating them (faster);
## but then this script would always have to be run after S4a. As it is now, it
## works as a stand alone script, except I'm not saving avgs so I should always
## run S4a as well.
#
#
## Start timer
#tic=timeit.default_timer()
#
#
## Get list of files: WT & control / lidocaine
#if not 'folderListFile_WT' in locals():
#    folderListFile_WT = input('Path for R4(CONTROL)/WT/WT folder list txt file: ').strip("'")
#
#if not 'folderListFile_expt' in locals():
#    folderListFile_expt = input('Path for CONTROL/LIDOC/AITC folder list txt file: ').strip("'")
#
#
## Define expt type: 'WT' or 'control' (laser offset) or 'lidocaine'
#if not 'exptType' in locals():
#    exptType = input('Experiment type: control, lidocaine or AITC? ')
#
#if exptType == 'control':
#    expt = 'control_'
#elif exptType == 'lidocaine':
#    expt = 'lidocaine_'
#elif exptType == 'AITC':
#    expt = 'AITC_'
#
#
## Specify data file root: WT & control / lidocaine
#data_path_WT,folderNames_WT=readFolderList(folderListFile_WT)
#data_path_expt,folderNames_expt=readFolderList(folderListFile_expt)
#
#
################################### 'WT' #######################################
#
#AvgCumulAngFiles = []
#AvgCurvFiles = []
#AvgMotFiles = []
#AvgLatencyFiles = []
#AvgRespFreqPerIntFiles = []
#AvgBoutPerIntFiles = []
#AvgRespFreqFiles = []
#AvgFirstVeloPerIntFiles = []
#AvgAvgVeloPerIFiles = []
#
#AvgmaxRdAngFiles = []
#AvgmaxRdCurvFiles = []
#AvgmaxRdMotFiles = []
#AvgmaxRdVeloFiles = []
#AvgTmaxRdAngFiles = []
#AvgTmaxRdCurvFiles = []
#AvgTmaxRdMotFiles = []
#AvgTmaxRdVeloFiles = []
#
#RF_AvgLatencyFiles = []
#RF_AvgRespFreqPerIntFiles = []
#RF_AvgFirstVeloPerIntFiles = []
#RF_AvgmaxRdAngFiles = []
#RF_AvgmaxRdCurvFiles = []
#RF_AvgmaxRdMotFiles = []
#RF_AvgTmaxRdAngFiles = []
#RF_AvgTmaxRdCurvFiles = []
#RF_AvgTmaxRdMotFiles = []
#  
#  
#if exptType == 'lidocaine' or 'AITC':
#  
#    for idF, folder in enumerate(folderNames_WT):
#        
#        # Go into analysis folder, find avg files
#        AvgCumulAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdAng.csv'))
#        AvgCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdCurv.csv'))
#        AvgMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdMot.csv'))
#        
#        AvgLatencyFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdLat.csv'))
#        AvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFreqPI.csv'))
#        AvgBoutPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdBout.csv'))
#        AvgRespFreqFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFreq.csv'))
#        AvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFirstVelo.csv'))
#        AvgAvgVeloPerIFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdAvVelo.csv'))
#        
#        AvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdAng.csv'))
#        AvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdCurv.csv'))
#        AvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdMot.csv'))
#        AvgmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdVelo.csv'))
#        
#        AvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdAng.csv'))
#        AvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdCurv.csv'))
#        AvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdMot.csv'))
#        AvgTmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgBMaxRdVelo.csv'))
#        
#        RF_AvgLatencyFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdLat.csv'))
#        RF_AvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdFreqPI.csv'))
#        RF_AvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdFirstVelo.csv'))
#        RF_AvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdAng.csv'))
#        RF_AvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdCurv.csv'))
#        RF_AvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdMot.csv'))
#        RF_AvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdAng.csv'))
#        RF_AvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdCurv.csv'))
#        RF_AvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdMot.csv'))    
#     
#        
#    ## Average cut, max and time to max files from each fish ##
#    _, avgAllAng, _, semAllAng = AverageFish(AvgCumulAngFiles,6)
#    _, avgAllCurv, _, semAllCurv = AverageFish(AvgCurvFiles,6)
#    _, avgAllMot, _, semAllMot = AverageFish(AvgMotFiles,6)
#    
#    _, avgAllLat, _, semAllLat = AverageFishParam(AvgLatencyFiles,6)
#    _, avgAllFreqPI, _, semAllFreqPI = AverageFishParam(AvgRespFreqPerIntFiles,6)
#    _, avgAllBout, _, semAllBout = AverageFishParam(AvgBoutPerIntFiles,6)
#    _, avgAllFreq, _, semAllFreq = AverageFishSingle(AvgRespFreqFiles)
#    _, avgAllFVelo, _, semAllFVelo = AverageFishParam(AvgFirstVeloPerIntFiles,6)
#    _, avgAllAvVelo, _, semAllAvVelo = AverageFishParam(AvgAvgVeloPerIFiles,6)
#    
#    _, avgAllmaxAng, _, semAllmaxAng = AverageFishParam(AvgmaxRdAngFiles,6)
#    _, avgAllmaxCurv, _, semAllmaxCurv = AverageFishParam(AvgmaxRdCurvFiles,6)
#    _, avgAllmaxMot, _, semAllmaxMot = AverageFishParam(AvgmaxRdMotFiles,6)
#    _, avgAllmaxVelo, _, semAllmaxVelo = AverageFishParam(AvgmaxRdVeloFiles,6)
#    
#    _, avgAllTmaxAng, _, semAllTmaxAng = AverageFishParam(AvgTmaxRdAngFiles,6)
#    _, avgAllTmaxCurv, _, semAllTmaxCurv = AverageFishParam(AvgTmaxRdCurvFiles,6)
#    _, avgAllTmaxMot, _, semAllTmaxMot = AverageFishParam(AvgTmaxRdMotFiles,6)
#    _, avgAllBmaxVelo, _, semAllBmaxVelo = AverageFishParam(AvgTmaxRdVeloFiles,6)
#    
#    _, RF_avgAllLat, _, RF_semAllLat = AverageFishParam(RF_AvgLatencyFiles,6)
#    _, RF_avgAllFreqPI, _, RF_semAllFreqPI = AverageFishParam(RF_AvgRespFreqPerIntFiles,6)
#    _, RF_avgAllFVelo, _, RF_semAllFVelo = AverageFishParam(RF_AvgFirstVeloPerIntFiles,6)
#    _, RF_avgAllmaxAng, _, RF_semAllmaxAng = AverageFishParam(RF_AvgmaxRdAngFiles,6)
#    _, RF_avgAllmaxCurv, _, RF_semAllmaxCurv = AverageFishParam(RF_AvgmaxRdCurvFiles,6)
#    _, RF_avgAllmaxMot, _, RF_semAllmaxMot = AverageFishParam(RF_AvgmaxRdMotFiles,6)
#    _, RF_avgAllTmaxAng, _, RF_semAllTmaxAng = AverageFishParam(RF_AvgTmaxRdAngFiles,6)
#    _, RF_avgAllTmaxCurv, _, RF_semAllTmaxCurv = AverageFishParam(RF_AvgTmaxRdCurvFiles,6)
#    _, RF_avgAllTmaxMot, _, RF_semAllTmaxMot = AverageFishParam(RF_AvgTmaxRdMotFiles,6)
#    
#    # SEM traces per parameter
#    NegSemAng,PosSemAng = GetSemTraces(avgAllAng,semAllAng)
#    NegSemCurv,PosSemCurv = GetSemTraces(avgAllCurv,semAllCurv)
#    NegSemMot,PosSemMot = GetSemTraces(avgAllMot,semAllMot)
#
#
#if exptType == 'control':
#    
#    for idF, folder in enumerate(folderNames_WT):
#    
#        AvgCumulAngFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_CumulAngCut.csv'))
#        AvgCurvFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_CurvCut.csv'))
#        AvgMotFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_MotCut.csv'))
#        
#        AvgLatencyFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_LatencyPerInt.csv'))
#        AvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_RespFreqPerInt.csv'))
#        AvgBoutPerIntFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_BoutsPerInt.csv'))
#        AvgRespFreqFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_RdRespFreq.csv'))
#        AvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_FirstVelocityPerInt.csv'))
#        AvgAvgVeloPerIFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_AvgVelocityPerInt.csv'))
#        
#        AvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_maxRdAng.csv'))
#        AvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_maxRdCurv.csv'))
#        AvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_maxRdMot.csv'))
#        AvgmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_maxRdVelo.csv'))
#        
#        AvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_TimeTOmaxRdAng.csv'))
#        AvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_TimeTOmaxRdCurv.csv'))
#        AvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_TimeTOmaxRdMot.csv'))
#        AvgTmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_BoutOFmaxRdVelo.csv'))
#        
#        RF_AvgLatencyFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexLatencyPerInt.csv'))
#        RF_AvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexRespFreqPerInt.csv'))
#        RF_AvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexFirstVelocityPerInt.csv'))
#        RF_AvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexmaxRdAng.csv'))
#        RF_AvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexmaxRdCurv.csv'))
#        RF_AvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexmaxRdMot.csv'))
#        RF_AvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexTimeTOmaxRdAng.csv'))
#        RF_AvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexTimeTOmaxRdCurv.csv'))
#        RF_AvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/round_4' + '*_ReflexTimeTOmaxRdMot.csv'))
#    
#    _, avgAllAng, _, semAllAng = AverageFish(AvgCumulAngFiles,6)
#    _, avgAllCurv, _, semAllCurv = AverageFish(AvgCurvFiles,6)
#    _, avgAllMot, _, semAllMot = AverageFish(AvgMotFiles,6)
#    
#    _, avgAllLat, _, semAllLat = AverageFishParam(AvgLatencyFiles,6)
#    _, avgAllFreqPI, _, semAllFreqPI = AverageFishParam(AvgRespFreqPerIntFiles,6)
#    _, avgAllBout, _, semAllBout = AverageFishParam(AvgBoutPerIntFiles,6)
#    _, avgAllFreq, _, semAllFreq = AverageFishSingle(AvgRespFreqFiles)
#    _, avgAllFVelo, _, semAllFVelo = AverageFishParam(AvgFirstVeloPerIntFiles,6)
#    _, avgAllAvVelo, _, semAllAvVelo = AverageFishParam(AvgAvgVeloPerIFiles,6)
#    
#    _, avgAllmaxAng, _, semAllmaxAng = AverageFishParam(AvgmaxRdAngFiles,6)
#    _, avgAllmaxCurv, _, semAllmaxCurv = AverageFishParam(AvgmaxRdCurvFiles,6)
#    _, avgAllmaxMot, _, semAllmaxMot = AverageFishParam(AvgmaxRdMotFiles,6)
#    _, avgAllmaxVelo, _, semAllmaxVelo = AverageFishParam(AvgmaxRdVeloFiles,6)
#    
#    _, avgAllTmaxAng, _, semAllTmaxAng = AverageFishParam(AvgTmaxRdAngFiles,6)
#    _, avgAllTmaxCurv, _, semAllTmaxCurv = AverageFishParam(AvgTmaxRdCurvFiles,6)
#    _, avgAllTmaxMot, _, semAllTmaxMot = AverageFishParam(AvgTmaxRdMotFiles,6)
#    _, avgAllBmaxVelo, _, semAllBmaxVelo = AverageFishParam(AvgTmaxRdVeloFiles,6)
#    
#    _, RF_avgAllLat, _, RF_semAllLat = AverageFishParam(RF_AvgLatencyFiles,6)
#    _, RF_avgAllFreqPI, _, RF_semAllFreqPI = AverageFishParam(RF_AvgRespFreqPerIntFiles,6)
#    _, RF_avgAllFVelo, _, RF_semAllFVelo = AverageFishParam(RF_AvgFirstVeloPerIntFiles,6)
#    _, RF_avgAllmaxAng, _, RF_semAllmaxAng = AverageFishParam(RF_AvgmaxRdAngFiles,6)
#    _, RF_avgAllmaxCurv, _, RF_semAllmaxCurv = AverageFishParam(RF_AvgmaxRdCurvFiles,6)
#    _, RF_avgAllmaxMot, _, RF_semAllmaxMot = AverageFishParam(RF_AvgmaxRdMotFiles,6)
#    _, RF_avgAllTmaxAng, _, RF_semAllTmaxAng = AverageFishParam(RF_AvgTmaxRdAngFiles,6)
#    _, RF_avgAllTmaxCurv, _, RF_semAllTmaxCurv = AverageFishParam(RF_AvgTmaxRdCurvFiles,6)
#    _, RF_avgAllTmaxMot, _, RF_semAllTmaxMot = AverageFishParam(RF_AvgTmaxRdMotFiles,6)
#    
#    # SEM traces per parameter
#    NegSemAng,PosSemAng = GetSemTraces(avgAllAng,semAllAng)
#    NegSemCurv,PosSemCurv = GetSemTraces(avgAllCurv,semAllCurv)
#    NegSemMot,PosSemMot = GetSemTraces(avgAllMot,semAllMot)
#
#    
#
#################################### EXPT ######################################
#
#EXPTAvgCumulAngFiles = []
#EXPTAvgCurvFiles = []
#EXPTAvgMotFiles = []
#EXPTAvgLatencyFiles = []
#EXPTAvgRespFreqPerIntFiles = []
#EXPTAvgBoutPerIntFiles = []
#EXPTAvgRespFreqFiles = []
#EXPTAvgFirstVeloPerIntFiles = []
#EXPTAvgAvgVeloPerIFiles = []
#
#EXPTAvgmaxRdAngFiles = []
#EXPTAvgmaxRdCurvFiles = []
#EXPTAvgmaxRdMotFiles = []
#EXPTAvgmaxRdVeloFiles = []
#EXPTAvgTmaxRdAngFiles = []
#EXPTAvgTmaxRdCurvFiles = []
#EXPTAvgTmaxRdMotFiles = []
#EXPTAvgTmaxRdVeloFiles = []
#
#EXPTRF_AvgLatencyFiles = []
#EXPTRF_AvgRespFreqPerIntFiles = []
#EXPTRF_AvgFirstVeloPerIntFiles = []
#EXPTRF_AvgmaxRdAngFiles = []
#EXPTRF_AvgmaxRdCurvFiles = []
#EXPTRF_AvgmaxRdMotFiles = []
#EXPTRF_AvgTmaxRdAngFiles = []
#EXPTRF_AvgTmaxRdCurvFiles = []
#EXPTRF_AvgTmaxRdMotFiles = []
#
#
#for idF, folder in enumerate(folderNames_expt):
#    
#        EXPTAvgCumulAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdAng.csv'))
#        EXPTAvgCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdCurv.csv'))
#        EXPTAvgMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdMot.csv'))
#        
#        EXPTAvgLatencyFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdLat.csv'))
#        EXPTAvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFreqPI.csv'))
#        EXPTAvgBoutPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdBout.csv'))
#        EXPTAvgRespFreqFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFreq.csv'))
#        EXPTAvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdFirstVelo.csv'))
#        EXPTAvgAvgVeloPerIFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgRdAvVelo.csv'))
#        
#        EXPTAvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdAng.csv'))
#        EXPTAvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdCurv.csv'))
#        EXPTAvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdMot.csv'))
#        EXPTAvgmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgMaxRdVelo.csv'))
#        
#        EXPTAvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdAng.csv'))
#        EXPTAvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdCurv.csv'))
#        EXPTAvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgTMaxRdMot.csv'))
#        EXPTAvgTmaxRdVeloFiles.extend(glob.glob(folder+'/Analysis/Summary/*_avgBMaxRdVelo.csv'))
#        
#        EXPTRF_AvgLatencyFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdLat.csv'))
#        EXPTRF_AvgRespFreqPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdFreqPI.csv'))
#        EXPTRF_AvgFirstVeloPerIntFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgRdFirstVelo.csv'))
#        EXPTRF_AvgmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdAng.csv'))
#        EXPTRF_AvgmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdCurv.csv'))
#        EXPTRF_AvgmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgMaxRdMot.csv'))
#        EXPTRF_AvgTmaxRdAngFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdAng.csv'))
#        EXPTRF_AvgTmaxRdCurvFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdCurv.csv'))
#        EXPTRF_AvgTmaxRdMotFiles.extend(glob.glob(folder+'/Analysis/Summary/*_REFLEXavgTMaxRdMot.csv'))    
# 
#    
### Average cut, max and time to max files from each fish ##
#_, EXPTavgAllAng, _, EXPTsemAllAng = AverageFish(EXPTAvgCumulAngFiles,6)
#_, EXPTavgAllCurv, _, EXPTsemAllCurv = AverageFish(EXPTAvgCurvFiles,6)
#_, EXPTavgAllMot, _, EXPTsemAllMot = AverageFish(EXPTAvgMotFiles,6)
#
#_, EXPTavgAllLat, _, EXPTsemAllLat = AverageFishParam(EXPTAvgLatencyFiles,6)
#_, EXPTavgAllFreqPI, _, EXPTsemAllFreqPI = AverageFishParam(EXPTAvgRespFreqPerIntFiles,6)
#_, EXPTavgAllBout, _, EXPTsemAllBout = AverageFishParam(EXPTAvgBoutPerIntFiles,6)
#_, EXPTavgAllFreq, _, EXPTsemAllFreq = AverageFishSingle(EXPTAvgRespFreqFiles)
#_, EXPTavgAllFVelo, _, EXPTsemAllFVelo = AverageFishParam(EXPTAvgFirstVeloPerIntFiles,6)
#_, EXPTavgAllAvVelo, _, EXPTsemAllAvVelo = AverageFishParam(EXPTAvgAvgVeloPerIFiles,6)
#
#_, EXPTavgAllmaxAng, _, EXPTsemAllmaxAng = AverageFishParam(EXPTAvgmaxRdAngFiles,6)
#_, EXPTavgAllmaxCurv, _, EXPTsemAllmaxCurv = AverageFishParam(EXPTAvgmaxRdCurvFiles,6)
#_, EXPTavgAllmaxMot, _, EXPTsemAllmaxMot = AverageFishParam(EXPTAvgmaxRdMotFiles,6)
#_, EXPTavgAllmaxVelo, _, EXPTsemAllmaxVelo = AverageFishParam(EXPTAvgmaxRdVeloFiles,6)
#
#_, EXPTavgAllTmaxAng, _, EXPTsemAllTmaxAng = AverageFishParam(EXPTAvgTmaxRdAngFiles,6)
#_, EXPTavgAllTmaxCurv, _, EXPTsemAllTmaxCurv = AverageFishParam(EXPTAvgTmaxRdCurvFiles,6)
#_, EXPTavgAllTmaxMot, _, EXPTsemAllTmaxMot = AverageFishParam(EXPTAvgTmaxRdMotFiles,6)
#_, EXPTavgAllBmaxVelo, _, EXPTsemAllBmaxVelo = AverageFishParam(EXPTAvgTmaxRdVeloFiles,6)
#
#_, EXPTRF_avgAllLat, _, EXPTRF_semAllLat = AverageFishParam(EXPTRF_AvgLatencyFiles,6)
#_, EXPTRF_avgAllFreqPI, _, EXPTRF_semAllFreqPI = AverageFishParam(EXPTRF_AvgRespFreqPerIntFiles,6)
#_, EXPTRF_avgAllFVelo, _, EXPTRF_semAllFVelo = AverageFishParam(EXPTRF_AvgFirstVeloPerIntFiles,6)
#_, EXPTRF_avgAllmaxAng, _, EXPTRF_semAllmaxAng = AverageFishParam(EXPTRF_AvgmaxRdAngFiles,6)
#_, EXPTRF_avgAllmaxCurv, _, EXPTRF_semAllmaxCurv = AverageFishParam(EXPTRF_AvgmaxRdCurvFiles,6)
#_, EXPTRF_avgAllmaxMot, _, EXPTRF_semAllmaxMot = AverageFishParam(EXPTRF_AvgmaxRdMotFiles,6)
#_, EXPTRF_avgAllTmaxAng, _, EXPTRF_semAllTmaxAng = AverageFishParam(EXPTRF_AvgTmaxRdAngFiles,6)
#_, EXPTRF_avgAllTmaxCurv, _, EXPTRF_semAllTmaxCurv = AverageFishParam(EXPTRF_AvgTmaxRdCurvFiles,6)
#_, EXPTRF_avgAllTmaxMot, _, EXPTRF_semAllTmaxMot = AverageFishParam(EXPTRF_AvgTmaxRdMotFiles,6)
#
#
## SEM traces per parameter
#EXPTNegSemAng,EXPTPosSemAng = GetSemTraces(EXPTavgAllAng,EXPTsemAllAng)
#EXPTNegSemCurv,EXPTPosSemCurv = GetSemTraces(EXPTavgAllCurv,EXPTsemAllCurv)
#EXPTNegSemMot,EXPTPosSemMot = GetSemTraces(EXPTavgAllMot,EXPTsemAllMot)
#
#
#
################ Plots #########################################################
#
## Save it in analysis folder of expt (not WT)
#folderFile = open(folderListFile_expt, 'r')
#folderList = folderFile.readlines()
#analysisFolder_path = folderList[0][:-1] + "/Analysis_Fish/"
#if not os.path.exists(analysisFolder_path):
#    os.mkdir(analysisFolder_path)
#
#analysis_path = analysisFolder_path + expt + 'combo_'
#
#
## Convert frames to msec
#num_frames = len(avgAllAng[0])
#frame_rate = 400 #fps
#num_msec = num_frames/frame_rate*1000
#x_ax = np.arange(0,num_msec,(num_msec/num_frames)) # Create new x values
#x_axis = x_ax - 1000 # because I had 1000 msec baseline
#
## Define plot parameters 
#box_leg_x = 1
#box_leg_y = 0.94
#x_min = -500
#amps = ['300mA','350mA','400mA','450mA','500mA','540mA']
#amps_noMA = [300,350,400,450,500,540]
#fish = ["F1","F2","F3"]
#colours_round =  ["b","orange","g"]
#
#box_legT_x = 0.1
#box_legT_y = 1.1
#cels = u"\N{DEGREE SIGN}" + 'C'
#celsSec = cels + '/s'
#temps = ['35'+celsSec+'(41'+cels+')','43'+celsSec+'(46'+cels+')','50'+celsSec+'(48'+cels+')','54'+celsSec+'(51'+cels+')','59'+celsSec+'(54'+cels+')','72'+celsSec+'(61'+cels+')']
#
#if exptType == 'control':
#    legd_expt = ['Off-target','On-target']
#    
#elif exptType == 'lidocaine':
#    legd_expt = ['Control','Lidocaine']
#
#
#
## Using colourmaps to generate colour gradient
#cmap=cm.get_cmap('hot')
#rgb=np.linspace(0.7,0,len(amps))
## NB: alternative attempts to generate gradient can be found in older versions
#
## Change image resolution and font to Arial
#plt.rcParams["font.family"] = "Arial"
#set_dpi = 300
#
#
##### Plotting AVG PER PARAMETER for ALL INTENSITIES same plot ####
## Plot each parameter (avg) w/out sd, save 1 tiff file
#figure=plt.figure(figsize=(10,8))
#PlotAllFishAllParamAvg(EXPTavgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",avgAllAng)
#PlotAllFishAllParamAvg(EXPTavgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",avgAllCurv)
#PlotAllFishAllParamAvg(EXPTavgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",avgAllMot)
#figname = 'fish_avgs_short.tiff'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi)
#
#figure=plt.figure(figsize=(10,8))
#PlotAllFishAllParamAvg(EXPTavgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",avgAllAng)
#PlotAllFishAllParamAvg(EXPTavgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",avgAllCurv)
#PlotAllFishAllParamAvg(EXPTavgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",avgAllMot)
#figname = 'fish_avgs.tiff'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi)
#
#
## Smooth
#
## Plot each parameter (avg) w/out sd, save 1 tiff file
#figure=plt.figure(figsize=(10,8))
#PlotAllFishAllParamAvgSm(EXPTavgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",avgAllAng)
#PlotAllFishAllParamAvgSm(EXPTavgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",avgAllCurv)
#PlotAllFishAllParamAvgSm(EXPTavgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",avgAllMot)
#figname = 'smooth_fish_avgs_short.tiff'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi)
#
#figure=plt.figure(figsize=(10,8))
#PlotAllFishAllParamAvgSm(EXPTavgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",avgAllAng)
#PlotAllFishAllParamAvgSm(EXPTavgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",avgAllCurv)
#PlotAllFishAllParamAvgSm(EXPTavgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",avgAllMot)
#figname = 'smooth_fish_avgs.tiff'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi)
#
#
#
#if exptType == 'lidocaine':
#    # if lidocaine, plot WT but reversing what the expt data and the WT data
#    # is -> takes ymin and ymax from lidocaine and makes WT bigger
#    
#    #### Plotting AVG PER PARAMETER for ALL INTENSITIES same plot ####
#    # Plot each parameter (avg) w/out sd, save 1 tiff file
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",EXPTavgAllAng)
#    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",EXPTavgAllCurv)
#    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",EXPTavgAllMot)    
#    figname = 'WT_fish_avgs_short.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",EXPTavgAllAng)
#    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",EXPTavgAllCurv)
#    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",EXPTavgAllMot)    
#    figname = 'WT_fish_avgs.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    
#    # Smooth
#    
#    # Plot each parameter (avg) w/out sd, save 1 tiff file
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvgSm(avgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",EXPTavgAllAng)
#    PlotAllFishAllParamAvgSm(avgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",EXPTavgAllCurv)
#    PlotAllFishAllParamAvgSm(avgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",EXPTavgAllMot)    
#    figname = 'WT_smooth_fish_avgs_short.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvgSm(avgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",EXPTavgAllAng)
#    PlotAllFishAllParamAvgSm(avgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",EXPTavgAllCurv)
#    PlotAllFishAllParamAvgSm(avgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",EXPTavgAllMot) 
#    figname = 'WT_smooth_fish_avgs.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#
#
#if exptType == 'control':
#    # if control, plot R4 using ylim from R4
#    
#    #### Plotting AVG PER PARAMETER for ALL INTENSITIES same plot ####
#    # Plot each parameter (avg) w/out sd, save 1 tiff file
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",avgAllAng)
#    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",avgAllCurv)
#    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",avgAllMot)    
#    figname = 'R4_fish_avgs_short.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvg(avgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",avgAllAng)
#    PlotAllFishAllParamAvg(avgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",avgAllCurv)
#    PlotAllFishAllParamAvg(avgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",avgAllMot)    
#    figname = 'R4_fish_avgs.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    
#    # Smooth
#    
#    # Plot each parameter (avg) w/out sd, save 1 tiff file
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvgSm(avgAllAng,"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles For All Fish",avgAllAng)
#    PlotAllFishAllParamAvgSm(avgAllCurv,"Curvature (deg)",2,2500,"Average Curvatures For All Fish",avgAllCurv)
#    PlotAllFishAllParamAvgSm(avgAllMot,"Motion (a.u.)",3,2500,"Average Motion For All Fish",avgAllMot)    
#    figname = 'R4_fish_avgs_short.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#    
#    figure=plt.figure(figsize=(10,8))
#    PlotAllFishAllParamAvgSm(avgAllAng,"Cumul. Angles (deg)",1,10500,"Average Cumulative Angles For All Fish",avgAllAng)
#    PlotAllFishAllParamAvgSm(avgAllCurv,"Curvature (deg)",2,10500,"Average Curvatures For All Fish",avgAllCurv)
#    PlotAllFishAllParamAvgSm(avgAllMot,"Motion (a.u.)",3,10500,"Average Motion For All Fish",avgAllMot)    
#    figname = 'R4_fish_avgs.tiff'
#    figpath = analysis_path + figname
#    figure.tight_layout(pad=1)
#    if os.path.isfile(figpath):
#        os.remove(figpath)
#    plt.savefig(figpath,dpi=set_dpi)
#
#
#
#
######################################################################################################################### STATS AND PLOT LEGEND
#
## Plot latency, resp freq, bout #, velocity vs Laser driver current
#figure=plt.figure(figsize=(10,8))
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllLat,semAllLat,EXPTavgAllLat,EXPTsemAllLat,1,"First response latency (ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllFreqPI,semAllFreqPI,EXPTavgAllFreqPI,EXPTsemAllFreqPI,2,"Response probability")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllBout,semAllBout,EXPTavgAllBout,EXPTsemAllBout,3,"Total number of bouts")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllFVelo,semAllFVelo,EXPTavgAllFVelo,EXPTsemAllFVelo,4,"First response vigour (deg/ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllAvVelo,semAllAvVelo,EXPTavgAllAvVelo,EXPTsemAllAvVelo,5,"Average response vigour (deg/ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,avgAllmaxVelo,semAllmaxVelo,EXPTavgAllmaxVelo,EXPTsemAllmaxVelo,6,"Peak response vigour (deg/ms)")
#figname = 'avg_parameters_per_intens.pdf'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi,transparent=True)
#
## Plot latency, resp freq, bout #, velocity vs Laser driver current
#figure=plt.figure(figsize=(10,8))
#PlotAFAvgLatFreqBoutVelo(amps_noMA,RF_avgAllLat,RF_semAllLat,EXPTRF_avgAllLat,EXPTRF_semAllLat,1,"SL response latency (ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,RF_avgAllFreqPI,RF_semAllFreqPI,EXPTRF_avgAllFreqPI,EXPTRF_semAllFreqPI,3,"SL Response probability")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,RF_avgAllFVelo,RF_semAllFVelo,EXPTRF_avgAllFVelo,EXPTRF_semAllFVelo,5,"SL response vigour (deg/ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,[],[],[],[],2,"First response latency (ms)")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,[],[],[],[],4,"Response probability")
#PlotAFAvgLatFreqBoutVelo(amps_noMA,[],[],[],[],6,"First response vigour (deg/ms)")
#figname = 'REFLEX_avg_parameters_per_intens.pdf'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi,transparent=True)
#
#
## Plot avg max + time to max vs Laser driver current for each parameter, save 1 tiff
#figure=plt.figure(figsize=(10,8))
#PlotAFAllParamAvgMaxPerInt(AvgmaxRdAngFiles,amps_noMA,avgAllmaxAng,semAllmaxAng,EXPTAvgmaxRdAngFiles,EXPTavgAllmaxAng,EXPTsemAllmaxAng,1,"Peak Cumul Angle (deg)")
#PlotAFAllParamAvgMaxPerInt(AvgmaxRdCurvFiles,amps_noMA,avgAllmaxCurv,semAllmaxCurv,EXPTAvgmaxRdCurvFiles,EXPTavgAllmaxCurv,EXPTsemAllmaxCurv,3,"Peak Curvature (deg)")
#PlotAFAllParamAvgMaxPerInt(AvgmaxRdMotFiles,amps_noMA,avgAllmaxMot,semAllmaxMot,EXPTAvgmaxRdMotFiles,EXPTavgAllmaxMot,EXPTsemAllmaxMot,5,"Peak Motion (a.u.)")
#PlotAFAllParamAvgMaxPerInt(AvgTmaxRdAngFiles,amps_noMA,avgAllTmaxAng,semAllTmaxAng,EXPTAvgTmaxRdAngFiles,EXPTavgAllTmaxAng,EXPTsemAllTmaxAng,2,"Time to Peak Cumul Angle (ms)")
#PlotAFAllParamAvgMaxPerInt(AvgTmaxRdCurvFiles,amps_noMA,avgAllTmaxCurv,semAllTmaxCurv,EXPTAvgTmaxRdCurvFiles,EXPTavgAllTmaxCurv,EXPTsemAllTmaxCurv,4,"Time to Peak Curvature (ms)")
#PlotAFAllParamAvgMaxPerInt(AvgTmaxRdMotFiles,amps_noMA,avgAllTmaxMot,semAllTmaxMot,EXPTAvgTmaxRdMotFiles,EXPTavgAllTmaxMot,EXPTsemAllTmaxMot,6,"Time to Peak Motion (ms)")
#figname = 'max_per_intens.pdf'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi,transparent=True)
#
## Plot avg max + time to max vs Laser driver current for each parameter, save 1 tiff
#figure=plt.figure(figsize=(10,8))
#PlotAFAllParamAvgMaxPerInt(RF_AvgmaxRdAngFiles,amps_noMA,RF_avgAllmaxAng,RF_semAllmaxAng,EXPTRF_AvgmaxRdAngFiles,EXPTRF_avgAllmaxAng,EXPTRF_semAllmaxAng,1,"SL Peak Cumul Angle (deg)")
#PlotAFAllParamAvgMaxPerInt(RF_AvgmaxRdCurvFiles,amps_noMA,RF_avgAllmaxCurv,RF_semAllmaxCurv,EXPTRF_AvgmaxRdCurvFiles,EXPTRF_avgAllmaxCurv,EXPTRF_semAllmaxCurv,3,"SL Peak Curvature (deg)")
#PlotAFAllParamAvgMaxPerInt(RF_AvgmaxRdMotFiles,amps_noMA,RF_avgAllmaxMot,RF_semAllmaxMot,EXPTRF_AvgmaxRdMotFiles,EXPTRF_avgAllmaxMot,EXPTRF_semAllmaxMot,5,"SL Peak Motion (a.u.)")
#PlotAFAllParamAvgMaxPerInt(RF_AvgTmaxRdAngFiles,amps_noMA,RF_avgAllTmaxAng,RF_semAllTmaxAng,EXPTRF_AvgTmaxRdAngFiles,EXPTRF_avgAllTmaxAng,EXPTRF_semAllTmaxAng,2,"SL Time to Peak Cumul Angle (ms)")
#PlotAFAllParamAvgMaxPerInt(RF_AvgTmaxRdCurvFiles,amps_noMA,RF_avgAllTmaxCurv,RF_semAllTmaxCurv,EXPTRF_AvgTmaxRdCurvFiles,EXPTRF_avgAllTmaxCurv,EXPTRF_semAllTmaxCurv,4,"SL Time to Peak Curvature (ms)")
#PlotAFAllParamAvgMaxPerInt(RF_AvgTmaxRdMotFiles,amps_noMA,RF_avgAllTmaxMot,RF_semAllTmaxMot,EXPTRF_AvgTmaxRdMotFiles,EXPTRF_avgAllTmaxMot,EXPTRF_semAllTmaxMot,6,"SL Time to Peak Motion (ms)")
#figname = 'REFLEX_max_per_intens.pdf'
#figpath = analysis_path + figname
#figure.tight_layout(pad=1)
#if os.path.isfile(figpath):
#    os.remove(figpath)
#plt.savefig(figpath,dpi=set_dpi,transparent=True)
#
#
#plt.close('all')
#
#
#
#
###################################################################################################### STATS ######################################
#
#whichFiles = [AvgLatencyFiles,
#              AvgRespFreqPerIntFiles,
#              AvgBoutPerIntFiles,
#              AvgFirstVeloPerIntFiles,
#              AvgAvgVeloPerIFiles,
#              AvgmaxRdAngFiles,
#              AvgmaxRdCurvFiles,
#              AvgmaxRdMotFiles,
#              AvgmaxRdVeloFiles,
#              AvgTmaxRdAngFiles,
#              AvgTmaxRdCurvFiles,
#              AvgTmaxRdMotFiles,
#              RF_AvgLatencyFiles,
#              RF_AvgRespFreqPerIntFiles,
#              RF_AvgFirstVeloPerIntFiles,
#              RF_AvgmaxRdAngFiles,
#              RF_AvgmaxRdCurvFiles,
#              RF_AvgmaxRdMotFiles,
#              RF_AvgTmaxRdAngFiles,
#              RF_AvgTmaxRdCurvFiles,
#              RF_AvgTmaxRdMotFiles]
#
#whichFP = ['f_p_Lat.csv',
#           'f_p_Freq.csv',
#           'f_p_BoutPI.csv',
#           'f_p_FVelo.csv',
#           'f_p_AvVelo.csv',
#           'f_p_maxAng.csv',
#           'f_p_maxCurv.csv',
#           'f_p_maxMot.csv',
#           'f_p_maxVelo.csv',
#           'f_p_TmaxAng.csv',
#           'f_p_TmaxCurv.csv',
#           'f_p_TmaxMot.csv',
#           'REFLEX_f_p_Lat.csv',
#           'REFLEX_f_p_Freq.csv',
#           'REFLEX_f_p_FVelo.csv',
#           'REFLEX_f_p_maxAng.csv',
#           'REFLEX_f_p_maxCurv.csv',
#           'REFLEX_f_p_maxMot.csv',
#           'REFLEX_f_p_TmaxAng.csv',
#           'REFLEX_f_p_TmaxCurv.csv',
#           'REFLEX_f_p_TmaxMot.csv']
#
#whichDF = ['df_Lat.csv',
#           'df_Freq.csv',
#           'df_BoutPI.csv',
#           'df_FVelo.csv',
#           'df_AvVelo.csv',
#           'df_maxAng.csv',
#           'df_maxCurv.csv',
#           'df_maxMot.csv',
#           'df_maxVelo.csv',
#           'df_TmaxAng.csv',
#           'df_TmaxCurv.csv',
#           'df_TmaxMot.csv',
#           'REFLEX_df_Lat.csv',
#           'REFLEX_df_Freq.csv',
#           'REFLEX_df_FVelo.csv',
#           'REFLEX_df_maxAng.csv',
#           'REFLEX_df_maxCurv.csv',
#           'REFLEX_df_maxMot.csv',
#           'REFLEX_df_TmaxAng.csv',
#           'REFLEX_df_TmaxCurv.csv',
#           'REFLEX_df_TmaxMot.csv']
#
#for parameter in range(len(whichFiles)):
#    fp,df=MyStats(whichFiles[parameter])
#    arr_path = analysis_path  + whichFP[parameter]
#    df_path = analysis_path  + whichDF[parameter]
#    pd.DataFrame(fp).to_csv(arr_path, header=None, index=None)
#    pd.DataFrame(df).to_csv(df_path, header=None, index=None)
#
#
## NB: there was a mistake here: i was saving fpTmaxAng & dfTmaxAng, rather than
## RF_fpTmaxAng and RF_dfTmaxAng
## RF_fpTmaxAng,RF_dfTmaxAng=MyStats(RF_AvgTmaxRdAngFiles)
## RF_arrTmaxAng_path = analysis_path  + 'REFLEX_f_p_TmaxAng.csv'
## RF_dfTmaxAng_path = analysis_path  + 'REFLEX_df_TmaxAng.csv'
## pd.DataFrame(fpTmaxAng).to_csv(RF_arrTmaxAng_path, header=None, index=None)
## pd.DataFrame(dfTmaxAng).to_csv(RF_dfTmaxAng_path, header=None, index=None)
#
#
## For stats comparing reflex to normal response - go to S5
#
## End timer
#toc = timeit.default_timer()
#time_elapsed_sec = toc-tic
#time_elapsed_min = time_elapsed_sec/60
#time_elapsed_hr = time_elapsed_min/60
#
#time_elapsed_sec_r = round(time_elapsed_sec)
#time_elapsed_min_r = round(time_elapsed_min)
#time_elapsed_hr_r = round(time_elapsed_hr,1)
#
#print('Finished. Time elapsed:' + str(time_elapsed_sec_r) + ' sec, i.e. ~' + str(time_elapsed_min_r) + ' mins, i.e. ~' + str(time_elapsed_hr_r) + 'h')
#
#
## FIN #
#
