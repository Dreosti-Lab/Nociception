# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:15:11 2021

@author: thoma
"""

import os
# Set Library Paths
lib_path=r'C:\Users\thoma\OneDrive\Documents\GitHub\Nociceptive_Zebrafish\libs'
import sys
sys.path.append(lib_path)
import NZ_utilities as NZU
import NZ_dictionary as NZD
import numpy as np
import matplotlib.pyplot as plt
import NZ_utilities as NZU
import NZ_ROITools as NZR
import cv2

vidPath='S:/WIBR_Dreosti_Lab/Elisa/Channelopathies/CRISPR/Prelim_results/NGF_KO/6dpf_KO_heat_A1-E6-hot.avi'
ROI_masks=np.load('C:/Users/thoma/OneDrive/Documents/GitHub/Nociceptive_Zebrafish/ROIMasks.npy')

vid = cv2.VideoCapture(vidPath)
width=int(vid.get(cv2.CAP_PROP_FRAME_WIDTH))
height=int(vid.get(cv2.CAP_PROP_FRAME_HEIGHT))
#numFrames = int(vid.get(cv2.CAP_PROP_FRAME_COUNT))-10 # skip last 10 (possibly corrupt) frames
numFrames=10000
Rows=['A','B','C','D','E']
Columns=['1','2','3','4','5','6']

ROINames=[]
for i in Rows:
    for k in Columns:
        ROINames.append(i+k)

# this module displays a sample frame of the selected movie to define new ROIs.
numROIs=len(ROINames)
#ROI_masks=[]

if len(ROI_masks)!=numROIs:
        print('Define your ROIS...')   
            
        # loop through ROIs, and define each one by polygon selection performed on a pop-up of video
        ROI_poly=[]
        ROI_masks=[]
        while len(ROI_poly)!=numROIs:
            for ROIName in ROINames:
                print('ROI : ' + ROIName)   
                polydata,mask=NZR.drawPoly(vidPath,ROIName) # ROI GUI
                ROI_poly.append(polydata.points)
                ROI_masks.append(mask)

# Now we compute a difference movie for the whole thing
vid.set(cv2.CAP_PROP_POS_FRAMES, 0)
#diffMovie=[]
ROIMotionTraces=[]
for thisROI in ROI_masks:
    temp=[]
    for f in range(numFrames):
#        print(f)
        ret,im=vid.read()
        im = np.uint8(cv2.cvtColor(im, cv2.COLOR_BGR2GRAY))
        im=im*thisROI
        im=im[im>0]
        if f==0:
            prev_frame=im
#            diffMovie.append(np.abs(np.subtract(prev_frame,im)))
            
        temp.append(np.sum(np.abs(np.subtract(prev_frame,im))))
        prev_frame=im
        
    ROIMotionTraces.append(temp-np.min(temp))
    
plt.figure()
for i in range (len(ROIMotionTraces)):
    plt.plot(ROIMotionTraces[i],color='black',alpha=0.2)
plt.plot(np.mean(ROIMotionTraces,axis=0),color='Black',alpha=1,linewidth=3)

        