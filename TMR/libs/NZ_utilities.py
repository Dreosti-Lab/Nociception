# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 15:53:27 2021

@author: Tom Ryan (Dreosti-Lab, UCL)
"""
## Collection of utility scripts for Nociceptive_Zebrafish repo
import os
import shutil
import glob 
import numpy as np
import pandas as pd
import cv2

def saveCSV(data,path,suff=''):
    if path[-3:]=='csv' and suff!='':
        path=path[0:-4]
    if path[-3:]!='csv' and (suff=='' or suff[-4:]!='.csv'):
        print('Warning! csv extension not set. Adding it automatically but you should check the file name')
        if suff[-4:]=='.csv':
            suff=suff[0:-4]
    path=path+suff+'.csv'
    pd.DataFrame(data).to_csv(path,header=None,index=None)
    return path

def checkInputs(folderListFile=[],exptType=[],dose=[],paired=[]):
    # Specify data file root
    if len(folderListFile)==0:
        folderListFile = input('Path for folder list txt file:').strip("'")    
    
    # Define expt type: 'WT' or 'control' (laser offset), 'lidocaine'
    if len(exptType)==0:
        exptType = input('Experiment type: WT, control, lidocaine or AITC? ')
        
    if exptType == 'WT':
        expt = ''
    
    elif exptType == 'control':
        expt = '_control'
    
    elif exptType == 'lidocaine':
        if len(dose)==0:
            dose = input('Dose? water, xmgL ')
        expt = '_lidocaine_' + dose
    
    elif exptType == 'AITC':
        if len(dose)==0:
            dose = input('Dose? water, DMSO, xuM ')
        expt = '_AITC_' + dose
    
    if (paired!=True and paired!=False) and (paired!='Yes' and paired!='No'):
        if len(paired)==0:
            paired = input('Was the experiment paired?: Yes / No ')
    numRounds = 3
    if paired==True or paired == 'Yes':
        numRounds = 6
        paired=True
    else: paired=False
    return folderListFile,expt,numRounds,paired

def readFolderList(folderListFile):
    folderFile=open(folderListFile, 'r')
    folderList=folderFile.readlines()
    data_path=folderList[0][:-1]
    folderList=folderList[1:]
    folderNames=[]
    for i,f in enumerate(folderList):
        stringline = f[:].split()
        expFolderName=data_path + stringline[0]
        folderNames.append(expFolderName) # create list of folders from folderlistfile
        
    return data_path,folderNames

def readIndFolderList(folderListFile):
    folderFile=open(folderListFile, 'r')
    folderList=folderFile.readlines()
    data_path=folderList[0][:-1]
    folderList=folderList[1:]
    folderNames=[]
    IndF=[]
    for i,f in enumerate(folderList):
        stringline = f[:].split()
        expFolderName=data_path + stringline[0]
        folderNames.append(expFolderName) # create list of folders from folderlistfile
        IndF.append(f[0:-1])
    return data_path,folderNames,IndF

def tryMkDir(path,report=0):
## Creates a new folder at the given path
## returns -1 and prints a warning if fails
## returns 1 if passed
    
    try:
        os.mkdir(path)
    except OSError:
        if(report):
            print ("Creation of the directory %s failed" % path + ", it might already exist!")
        return -1
    else:
        if(report):
            print ("Successfully created the directory %s " % path)
        return 1
        
def cycleMkDir(path,report=0):
## Creates folders and subfolder along defined path
## returns -1 and prints a warning if fails
## returns 1 if passed
    splitPath=path.split(sep="\\")
    for i,name in enumerate(splitPath):
        if(i==0):
            s=name+"\\"
        else:
            s=s+name+"\\"
        if(i!=len(splitPath)-1):    
            tryMkDir(s,report=report)
        else:
            tryMkDir(s,report=report)
            
def cycleSeg(inFiles,destFolder,copy=True,verbose=True):
    
        
    destS=[]
    for i,source in enumerate(inFiles):
        name=source.split(sep='\\',maxsplit=5)[-1]
        dest=destFolder+'/'+name
        destS.append(dest)
        if os.path.exists(dest): 
            print (dest + ' already exists... skipping') 
        else:
            cycleMkDir(dest.rsplit(sep='\\',maxsplit=1)[0])
            if copy:
                if verbose: print('Copying file ' + str(i) + ' of ' + str(len(inFiles)))
                shutil.copy(source,dest)
                if verbose: print('Done')
            else:
                if verbose:print('Moving file ' + str(i) + ' of ' + str(len(inFiles)))
                shutil.move(source,dest)
                if verbose: print('Done')
                
    return destS

def segregateSingleAnalysisFiles(inFiles,destFolder=[],copy=True,verbose=True):
    
    if len(destFolder)==0:
        inDir=inFiles[0].rsplit(sep='\\',maxsplit=7)[0]
        destFolder=inDir+'SegregatedAnalysisFiles\\'
    segFiles=cycleSeg(inFiles,destFolder,copy=copy,verbose=verbose)
    
    return segFiles
        
def segregateAnalysisFiles(folderListFile,destFolder=[],copy=True):
    wDir,folders,IndF=readIndFolderList(folderListFile)
    if len(destFolder)==0:destFolder=wDir+'SegregatedAnalysisFiles\\'
    
    for i,folder in enumerate(folders):
        xFiles=glob.glob(folder+'/Analysis/*centre_x.csv')
        yFiles=glob.glob(folder+'/Analysis/*centre_y.csv')
        stimFiles=glob.glob(folder+'\\*.csv')
        destFolderInd=destFolder+'\\'+IndF[i]+'\\Analysis\\'
        destFolderIndStim=destFolder+'\\'+IndF[i]+'\\'
        def cycleSeg(inFiles,destFolder):
            
            for i,source in enumerate(inFiles):
                _,name=source.rsplit(sep='\\',maxsplit=1)
                cycleMkDir(destFolder)
                dest=destFolder+'/'+name
                if copy:
                    shutil.copy(source,dest)
                else:
                    shutil.move(source,dest)
                
        cycleSeg(xFiles,destFolderInd)
        cycleSeg(yFiles,destFolderInd)
        cycleSeg(stimFiles,destFolderIndStim)    
#        
def cycleFixSeg(inFiles,destFolder,copy=True):
        
        for i,source in enumerate(inFiles):
            _,name=source.rsplit(sep='\\',maxsplit=1)
            cycleMkDir(destFolder)
            dest=destFolder+'\\'+name
            if copy:
                shutil.copy(source,dest)
            else:
                shutil.move(source,dest)

def fixSeg(folderListFile,copy=True):
    wDir,folders,IndF=readIndFolderList(folderListFile)
    for i,folder in enumerate(folders):
        stimFiles=[]
        chckFiles=glob.glob(folder+'\\Analysis\\*.csv')
        for chck in chckFiles:
            if chck[-5]!='x' and chck[-5]!='y' and chck[-5]!='v' and chck[-5]!='a':
                stimFiles.append(chck)
        destFolder=folder
        
        cycleFixSeg(stimFiles,destFolder)  

def grabFrame(avi,frame):
# grab frame and return the image from loaded cv2 movie
    vid=cv2.VideoCapture(avi)
    vid.set(cv2.CAP_PROP_POS_FRAMES, frame)
    ret, im = vid.read()
    vid.set(cv2.CAP_PROP_POS_FRAMES, 0)
    vid.release()
    im = np.uint8(cv2.cvtColor(im, cv2.COLOR_BGR2GRAY))
    return im
#segregateAnalysisFiles(folderListFile,destFolder='F:\ElisaDataSeg')
#segregateAnalysisFiles(folderListFile,destFolder=[],copy=False)