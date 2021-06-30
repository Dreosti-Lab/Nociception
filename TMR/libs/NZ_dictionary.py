# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:55:51 2021

@author: Tom Ryan (Dreosti-Lab,UCL)
"""

## Library for dictionary handling for Nociceptive_Zebrafish repo
lib_path='C:/Users/thoma/OneDrive/Documents/GitHub/Elisa/scripts_toRun/libs/'  ####### CHANGE THIS
import sys
sys.path.append(lib_path)
import NZ_utilities as NZU
import os
import shutil
import numpy as np


def createSingleFishDict(labels,values,savePath,metaDict=[],overwrite=True,report=True,save=True):
    fail,savePath=createDictionaryLocs(savePath,overwrite=overwrite)
    if fail:return -1,-1
    Dict=makeDictionary(labels=labels,values=values)
    if len(metaDict)!=0:
        Dict['MetaData']=metaDict
    else: 
        print('No metadata provided!')
    if save:
        np.save(savePath,Dict)
    if report:
        print('Dictionary saved at ' + savePath)
    return Dict,savePath

def mergeDictionariesAfter(newKeyList,dictList=[]):
    Dict={}
    if len(dictList)==0:
        dictList=[{}] # if dictionary not provided, then we are making a new dictionary, otherwise add to a new dictionary
    for i, dic in enumerate(dictList):
        if len(dic.keys())==0:
            afterKeyList=['NULL']
            print('Warning, creating new dictionary with no larger field names... You might not have wanted to do that')
        else: 
            afterKeyList=list(dic.keys())
                
        for afterKey in afterKeyList:
            try: 
                test=Dict[afterKey]
                del test
                Dict[afterKey].update({newKeyList[i]:dic[afterKey]})
            except:
                Dict.update({afterKey:{newKeyList[i]:dic[afterKey]}})
            
    return Dict

def mergeDictionariesBefore(newKeyList,dictList,Dict=[]):
    if len(Dict)==0:Dict={} # if dictionary not provided, then we are making a new dictionary, otherwise add to a new dictionary
    for i,newKey in enumerate(newKeyList):
        Dict[newKey]=dictList[i]
    return Dict

def makeDictionary(newKey=[],labels=[],values=[],Dict=[]):
    if len(Dict)==0:Dict={} # if dictionary not provided, then we are making a new dictionary, otherwise add to a new dictionary
    for thisValI,thisLab in enumerate(labels): # cycle through values and labels
        thisVal=values[thisValI]
        if len(newKey)!=0: # if a newKey has been given, all values and labels will be stored in a dictionary within the given or new dictionary under the newKey
            subDict=makeDictionary(labels=labels,values=values)
            Dict[newKey]=subDict
            break # break the for loop if we have added a dictionary to the dictionary
        else:
            Dict.update( {thisLab : thisVal} )
    return Dict
    
def createDictionaryLocs(savePath,overwrite=True): # fix this so it makes the old dict when overwrite set to false, and removes then creates new, overwriting existing dict if set to true
    wDir,name=savePath.rsplit(sep='\\',maxsplit=1)
    NZU.cycleMkDir(wDir)
    oldDictDir=wDir+'\\OldDictionaries\\'
    
    if os.path.exists(savePath) and overwrite==False:
        print('Aborting... Dictionary already exists at ' + savePath + '. Set overwite = True or change savePath')
        return True,-1
    if os.path.exists(savePath) and overwrite:
        print('Overwriting dictionary at ' + savePath + '. The old dictionary will be stored in seperate folder as *_OLD#.npy')
        if os.path.exists(oldDictDir)==False:
            NZU.cycleMkDir(oldDictDir)
            
        dest=oldDictDir+'/'+name[0:-4]+'_OLD'
        i=0
        while os.path.exists(dest): 
            i+1
            dest=dest[0:-4]+'_'+str(i)
        shutil.move(savePath,dest)
    elif  os.path.exists(savePath)==False:print('Creating new dictionary at ' + savePath)
    return False,savePath
    
    