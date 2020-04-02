# -*- coding: utf-8 -*-
"""
Measure tail parameters

@author: Adam Kampff & Tom Ryan & Elisa C
"""
#import os
#import cv2
import numpy as np
import pandas as pd 
#import matplotlib.pyplot as plt
import glob
#import timeit
#import csv as csv

### Helper Functions

def findStartEndAndCutEC(volt_data,uncut_data):
    diffvolt = np.diff(volt_data)
    starts = np.where(diffvolt > 0)[0] #because that yields a tuple
    starts = starts + 1
    volts = volt_data[starts]
    num_trials = len(starts)
    startFrame = starts - 400
    endFrame = starts + 4200
    
    cut_data = pd.DataFrame()
    for t in range(num_trials):
        cut_data[t]=uncut_data[startFrame[t]:endFrame[t]]
    
    cut_data = cut_data.T
    cut_data['Voltage']=volts
    cut_data = cut_data.sort_values(by=['Voltage'])
    return cut_data


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


### Algorithm #################################################################

# Specify data file root
    
folderListFile='C:\TestData\FolderListTest.txt'
data_path,folderNames=readFolderList(folderListFile)

# Loop through folders
for idF, folder in enumerate(folderNames):
    csvXFiles=glob.glob(folder+'/Analysis/*centre_x.csv')
    csvYFiles=glob.glob(folder+'/Analysis/*centre_y.csv')
    stimFiles=glob.glob(folder+'/*.csv')
    
    # Loop through files (rounds)
    for x,csvXFile in enumerate(csvXFiles):
        csvYFile=csvYFiles[x]
        stimFile=stimFiles[x]
        
        volt_data = np.genfromtxt(stimFile, delimiter=',')
        x_data = np.genfromtxt(csvXFile, delimiter=',')
        y_data = np.genfromtxt(csvYFile, delimiter=',')
        volt_data=volt_data[:,0]
        
  
        # Determine number of frames, segments and angles
        num_frames = np.shape(x_data)[0]
        num_segments = np.shape(x_data)[1]
        num_angles = num_segments - 1
        
        
        # Allocate space for measurements - creates empty arrays to store data
        cumulAngles = np.zeros(num_frames)
        curvatures = np.zeros(num_frames)
        motion = np.zeros(num_frames)
        
        # Measure tail motion, angles, and curvature
        
        # Take x and y values of first frame for each segment
        prev_xs = x_data[0, :]
        prev_ys = y_data[0, :]
  
        for f in range(num_frames):
            delta_thetas = np.zeros(num_angles) # Make an array of zeros with the same size as num angles (num seg-1)
            prev_theta = 0.0 # set first theta to zero
            for a in range(num_angles):
                dx = x_data[f, a+1] - x_data[f, a] # dx between each segment for the same frame
                dy = y_data[f, a+1] - y_data[f, a] # dy between each segment for the same frame
                theta = np.arctan2(dx, dy) * 360.0 / (2.0*np.pi) # calc arctangent bt dx and dy, convert to deg
                delta_thetas[a] = theta - prev_theta
                prev_theta = theta # prev theta is set to current theta
            cumulAngles[f] = np.sum(delta_thetas) # sum all angles for each frame
            curvatures[f] = np.mean(np.abs(delta_thetas)) # mean of abs vale of angles
        
            # Measure motion
            diff_xs = x_data[f,:] - prev_xs # difference between current x and prev x from each segment and each frame
            diff_ys = y_data[f,:] - prev_ys
            motion[f] = np.sum(np.sqrt(diff_xs*diff_xs + diff_ys*diff_ys)) # motion as sqrt of sq diff of x & y
            
            # Store previous tail
            prev_xs = x_data[f, :] # So that we don't always take the 1st frame of the movie to calculate motion
            prev_ys = y_data[f, :]
            
        # Report
        #plt.plot(cumulAngles, 'r')
        #plt.plot(curvatures, 'b')
        #plt.plot(motion, 'g')
        #plt.show()
        ####### END OF FRAME LOOP ###############  
            
        # Create csv paths
        round_root = csvXFile[0:-12]
        cumulAngles_path = round_root +'cumulative_angles.csv'
        curvatures_path = round_root +'curvatures.csv'
        motion_path = round_root +'motion.csv'
        cumulAnglesCut_path=round_root +'CumulAngCut.csv'
        curvaturesCut_path=round_root +'CurvCut.csv'
        motionCut_path=round_root +'MotCut.csv'
        
        # Cut data by volt epochs, add volt to dataframe, sort according to volt
        cumulAnglesCut = findStartEndAndCutEC(volt_data,cumulAngles)
        curvaturesCut = findStartEndAndCutEC(volt_data,curvatures)
        motionCut = findStartEndAndCutEC(volt_data,motion)
        
        # Save csv files
        pd.DataFrame(cumulAngles).to_csv(cumulAngles_path, header=None, index=None)
        pd.DataFrame(curvatures).to_csv(curvatures_path, header=None, index=None)
        pd.DataFrame(motion).to_csv(motion_path, header=None, index=None)
        pd.DataFrame(cumulAnglesCut).to_csv(cumulAnglesCut_path, header=None, index=None)
        pd.DataFrame(curvaturesCut).to_csv(curvaturesCut_path, header=None, index=None)
        pd.DataFrame(motionCut).to_csv(motionCut_path, header=None, index=None)
    ####### END OF FILE LOOP ##############
######## END OF FOLDER LOOP ###########

# FIN