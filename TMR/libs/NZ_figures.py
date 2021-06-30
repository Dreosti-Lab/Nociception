# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:50:56 2021

@author: Elisa C & Tom Ryan (Dreosti-Lab, UCL)
"""
# created as library for figures on 01/02/21 TR
import matplotlib.pyplot as plt
import numpy as np

def PlotAllParamAvgMaxPerInt(intens,avgMaxRd,SdMax,fignumb,ylabel):
    # Plot max vs intensity for each parameter, save as 1 tiff file
    plt.subplot(3,2,fignumb)
    plt.errorbar(intens,avgMaxRd,SdMax,None,c='k',linestyle='None',marker='.',markersize=9,capsize=3)
    # plt.scatter(intens,avgMaxRd,c='k')
    # plt.scatter(intens,NegSdMax,s=11,c='k',marker="_")
    # plt.scatter(intens,PosSdMax,s=11,c='k',marker="_")
    plt.xlabel('Intensity (mA)')
    plt.ylabel(ylabel)
    bottom, top = plt.ylim()
    new_bottom = min(bottom,0)
    new_top = top + 20
    plt.ylim(new_bottom,new_top)
    
def GetSdTraces(avData,sdData):
    all_neg_SD = []
    all_pos_SD = []
    
    for j in range(len(avData)):
        neg_SD = np.subtract(avData[j],sdData[j])
        pos_SD = np.add(avData[j],sdData[j])
        all_neg_SD.append(neg_SD)
        all_pos_SD.append(pos_SD)
    
    return all_neg_SD, all_pos_SD