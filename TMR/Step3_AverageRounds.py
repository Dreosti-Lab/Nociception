# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:49:28 2021

@author: Tom Ryan, Adam Kampff, Elisa C  
Compiled and edited, Tom Ryan
"""

############ Set default inputs and options DO NOT CHANGE ##########
folderListFile=[]
Dose=''
#####################################################
# Set "Library Path" - Nociceptive Zebrafish Repo
lib_path='C:/Users/thoma/OneDrive/Documents/GitHub/Nociception/TMR/libs/'
import sys
sys.path.append(lib_path)

# Generic Libraries
import glob
import numpy as np
import os
import pandas as pd
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
paired='Yes'
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
if paired == 'No':
    timing = '_' + Dose
    whichRounds=['round_1','round_2','round_3']
    paired=False
elif paired == 'Yes':
    timing='pre_' + Dose
    preDrug = ['round_1','round_2','round_3']
    postDrug = ['round_4','round_5','round_6']
    paired=True
    
if grabDict: 
    print('grabbing dictionaries')
else: print('grabbing csvs')
if saveCSV and saveDict: print('saving as both csv files and dictionaries')
if saveCSV and not saveDict: print('saving as csv files')
if saveDict and not saveCSV: print('saving as dictionaries')

data_path,folderNames=NZU.readFolderList(folderListFile)

## start of folder Loop ##
for idF, folder in enumerate(folderNames):
    if ServerExchange:
        folder_serv=folder
        fpath=folder.split(sep='\\',maxsplit=5)[-1]
        folder=root+'Data/'+fpath
    ## grab corresponding metadata
    if grabDict:    
        DictFolder=folder+'\\Dictionaries\\'
        
        if paired:
            preDrugDictFiles=[]
            postDrugDictFiles=[]
            for i,k in enumerate(preDrug):
                ff= DictFolder + k + '*.npy'
                find_it = "\\"
                repl_it = "/"
                if ff.find(find_it):
                    ff = ff.replace(find_it, repl_it)
                preDrugDictFiles.append(glob.glob(ff)[0])
                ff= DictFolder + postDrug[i] + '*.npy'
                find_it = "\\"
                repl_it = "/"
                if ff.find(find_it):
                    ff = ff.replace(find_it, repl_it)
                postDrugDictFiles.append(glob.glob(ff)[0])
            
            
            preDrugkeys,preDrugaverages,preDrugSDs=NZR.AverageRoundsDict(preDrugDictFiles)
            postDrugkeys,postDrugaverages,postDrugSDs=NZR.AverageRoundsDict(postDrugDictFiles)
            
            preAVG_Dict=NZD.makeDictionary(labels=preDrugkeys,values=preDrugaverages)
            preSD_Dict=NZD.makeDictionary(labels=preDrugkeys,values=preDrugSDs)
            
            postAVG_Dict=NZD.makeDictionary(labels=postDrugkeys,values=postDrugaverages)
            postSD_Dict=NZD.makeDictionary(labels=postDrugkeys,values=postDrugSDs)
            
            preAVGDictList=[preAVG_Dict,preSD_Dict]
            preComb_Dict=NZD.mergeDictionariesAfter(['Mean','SD'],preAVGDictList)
            
            postAVGDictList=[postAVG_Dict,postSD_Dict]
            postComb_Dict=NZD.mergeDictionariesAfter(['Mean','SD'],postAVGDictList)
            del preAVGDictList, postAVGDictList
            
            if group and saveDict:
                print('Collecting rounds into paired dictionaries')
                if deleteIndDicts: print(' will delete individual files if saving')
                predictList,prenewKeyList=[],[]
                postdictList,postnewKeyList=[],[]
                for i in range(len(preDrugDictFiles)):
                    predictList.append(np.load(preDrugDictFiles[i],allow_pickle=True).item())
                    postdictList.append(np.load(postDrugDictFiles[i],allow_pickle=True).item())
                    prenewKeyList.append('round_'+str(i+1))
                    postnewKeyList.append('round_'+str(i+4))
                    
            preAVG_Dict=NZD.mergeDictionariesBefore(prenewKeyList,predictList,Dict=preComb_Dict)
            postAVG_Dict=NZD.mergeDictionariesBefore(postnewKeyList,postdictList,Dict=postComb_Dict)
            pairedAVGDict=NZD.mergeDictionariesBefore(['preDrug','postDrug'],[preAVG_Dict,postAVG_Dict])
            
            if group and saveCSV:
                print('Not possible to save csvs from dictionary inputs yet. Overriding and saving as Dictionary')
                saveCSV=False
                saveDict=True
                
            if saveDict:
                # save paired pre vs post dictionary
                wDir,name=preDrugDictFiles[0].rsplit(sep='\\',maxsplit=1)
                names=name.rsplit(sep='_')[2:-1]
                saveName=wDir+'\\'+('_'.join(names))+'paired_AverageDict.npy'
                np.save(saveName,pairedAVGDict)
                
                print('Saved Fish dictionary at ' + saveName)
                if group and deleteIndDicts:
                    print('Deleting individual round files')
                    for thisOne in preDrugDictFiles:
                        os.remove(thisOne)
                    for thisOne in postDrugDictFiles:
                        os.remove(thisOne)
                        
        else: # not sure this is functional
            
            dictFiles=glob.glob(DictFolder+'*.npy')
            keys,averages,SDs=NZR.AverageRoundsDict(dictFiles)
            AVG_Dict=NZD.makeDictionary(labels=keys,values=averages)
            SD_Dict=NZD.makeDictionary(labels=keys,values=SDs)
            AVGDictList=[AVG_Dict,SD_Dict]
            AVG_Dict=NZD.mergeDictionariesAfter(['Mean','SD'],AVGDictList)
            del AVGDictList
            
            if group and saveDict:
                print('Collecting rounds into single dictionary')
                if deleteIndDicts: print(' will delete individual files if saving')
                dictList,newKeyList=[],[]
                for i in range(len(dictFiles)):
                    dictList.append(np.load(dictFiles[i],allow_pickle=True).item())
                    newKeyList.append('round_'+str(i))
                AVG_Dict=NZD.mergeDictionariesBefore(newKeyList,dictList,Dict=AVG_Dict)
            if saveDict:
                wDir,name=dictFiles[0].rsplit(sep='\\',maxsplit=1)
                names=name.rsplit(sep='_')[2:-1]
                saveName=wDir+'\\'+('_'.join(names))+'_AverageDict.npy'
                np.save(saveName,AVG_Dict)
                print('Saved Fish dictionary at ' + saveName)
                if group and deleteIndDicts:
                    print('Deleting individual round files')
                    for thisOne in dictList:
                        os.remove(thisOne)
            if saveCSV:
                print('Not possible to save csvs from dictionary inputs yet. Overriding and saving as Dictionary')
                saveCSV=False
    else:
        listofInSuffs=['CumulAngCut',
                       'CurvCut',
                       'MotCut',
                       'LatencyPerInt',
                       'RespFreqPerInt',
                       'BoutsPerInt',
                       'RdRespFreq',
        #               'VelocityPerBout',
                       'FirstVelocityPerInt',
                       'AvgVelocityPerInt',
                       'maxRdAng',
                       'maxRdCurv',
                       'maxRdMot',
                       'maxRdVelo',
                       'TimeTOmaxRdAng',
                       'TimeTOmaxRdCurv',
                       'TimeTOmaxRdMot',
                       'BoutOFmaxRdVelo',
                       'ReflexLatencyPerInt',
                       'ReflexRespFreqPerInt',
                       'ReflexFirstVelocityPerInt',
                       'ReflexmaxRdAng',
                       'ReflexmaxRdCurv',
                       'ReflexmaxRdMot',
                       'ReflexTimeTOmaxRdAng',
                       'ReflexTimeTOmaxRdCurv',
                       'ReflexTimeTOmaxRdMot']

        ##### names of output variables
        listofOutSuffs = ['avgRdAng',
                          'avgRandCurv',
                          'avgRdMot',
                          'avgRdLat',
                          'avgRdFreqPI',
        #                  'avgVelPerBout',
                          'avgRdBout',
                          'avgRdFreq',
                          'avgRdFirstVelo',
                          'avgRdAvVelo',
                          'avgMaxRdAng',
                          'avgMaxRdCurv',
                          'avgMaxRdMot',
                          'avgMaxRdVelo',
                          'avgTMaxRdAng',
                          'avgTMaxRdCurv',
                          'avgTMaxRdMot',
                          'avgBMaxRdVelo',
                          'REFLEXavgRdLat',
                          'REFLEXavgRdFreqPI',
                          'REFLEXavgRdFirstVelo',
                          'REFLEXavgMaxRdAng',
                          'REFLEXavgMaxRdCurv',
                          'REFLEXavgMaxRdMot',
                          'REFLEXavgTMaxRdAng',
                          'REFLEXavgTMaxRdCurv',
                          'REFLEXavgTMaxRdMot']

        ######### list of file strings to load/build csv files  
        csvIn=[]
        csvOut=[]
        for i,k in enumerate(listofInSuffs):
            csvIn.append('*_' + k + '.csv')
            csvOut.append('_' + listofOutSuffs[i] + '.csv')
            
        ##### Start run ########
        data_path,folderNames=NZU.readFolderList(folderListFile)
        for idF, folder in enumerate(folderNames):
            
            if ServerExchange:
                folder_serv=folder
                fpath=folder.split(sep='\\',maxsplit=5)[-1]
                folder=root+fpath
            
            pathD=folder+'\\Analysis\\'
        
            # type is 0 for rounds, 1 for single and 2 for parameters
            typeOfMetric=[0,0,0,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        
            # generate averaged rounds
            avgS_preDrug=[]
            stdevS_preDrug=[]
            avgS_postDrug=[]
            stdevS_postDrug=[]
    
            ###### Create csv paths ######
            # Create root path and root fname and new summary folder
            root_path = folder + '\\Analysis\\Summary\\'
            if not os.path.exists(root_path):
                os.mkdir(root_path)
        
            fname = folder[-18:] # works +/- w/ 2 digit dates & 1/2 digit fish #
            fname = fname.replace('_','')
            fname = fname.replace('\\','') + timing
        
            for i,thisSuff in enumerate(csvIn):
                path = pathD + thisSuff
                
                avg_preDrug, stdev_preDrug = NZR.grabCSVandAverage(path,preDrug,typeOfMetric[i])
                avgS_preDrug.append(avg_preDrug)
                stdevS_preDrug.append(stdev_preDrug)
            
                avg_postDrug, stdev_postDrug = NZR.grabCSVandAverage(path,postDrug,typeOfMetric[i])
                avgS_postDrug.append(avg_postDrug) # all averages are now stored in this single matrix: the order is the same as listofInSuffs
                stdevS_postDrug.append(stdev_postDrug)
            if saveCSV:
                for i,thisOutSuff in enumerate(listofOutSuffs):
                    path = root_path + fname + '_' + thisOutSuff + '_preDrugAVG.csv'
                    pd.DataFrame(avgS_preDrug[i]).to_csv(path, header=None, index=None) # save as csv as before
                    path = root_path + fname + '_' + thisOutSuff + '_preDrugSD.csv'
                    pd.DataFrame(stdevS_preDrug[i]).to_csv(path, header=None, index=None) # save as csv as before
                    path = root_path + fname + '_' + thisOutSuff + '_postDrugAVG.csv'
                    pd.DataFrame(avgS_postDrug[i]).to_csv(path, header=None, index=None) # save as csv as before
                    path = root_path + fname + '_' + thisOutSuff + '_postDrugSD.csv'
                    pd.DataFrame(stdevS_postDrug[i]).to_csv(path, header=None, index=None) # save as csv as before
        
        
            
if plot:
    
    figure=plt.figure(figsize=(10,8))
    NZF.PlotAllParamAvg(AVG_Dict['CumulAngCut']['Mean'],"Cumul. Angles (deg)",1,2500,"Average Cumulative Angles Per Fish")
    NZF.PlotAllParamAvg(AVG_Dict['CurvCut']['Mean'],"Curvature (deg)",2,2500,"Average Curvatures Per Fish")
    NZF.PlotAllParamAvg(AVG_Dict['MotCut']['Mean'],"Motion (a.u.)",3,2500,"Average Motion Per Fish")
    figPath=folder+outFolder
    figname = fname + '_avgs_short.tiff'
    figpath = figPath + figname
    figure.tight_layout(pad=1)
    if os.path.isfile(figpath):
        os.remove(figpath)
    plt.savefig(figpath,dpi=set_dpi)
    
    
    latency=AVG_Dict['LatencyPerInt']
    avgRdLat, stdevRdLat = motCut['Mean'],motCut['SD']

#avgRdLat, stdevRdLat = AverageRoundsParam(LatencyFiles,6)
#    _, avgRdFreqPI, stdevRdFreqPI = AverageRoundsParam(RespFreqPerIntFiles,6)
#    _, avgRdBout, stdevRdBout = AverageRoundsParam(BoutPerIntFiles,6)
#    _, avgRdFreq, stdevRdFreq = AverageRoundsSingle(RespFreqFiles)
#    _, avgRdFVelo, stdevRdFVelo = AverageRoundsParam(FirstVeloPerIntFiles,6)
#    _, avgRdAvVelo, stdevRdAvVelo = AverageRoundsParam(AvgVeloPerIFiles,6)
#    
#    _, avgRdmaxAng, stdevRdmaxAng = AverageRoundsParam(maxRdAngFiles,6)
#    _, avgRdmaxCurv, stdevRdmaxCurv = AverageRoundsParam(maxRdCurvFiles,6)
#    _, avgRdmaxMot, stdevRdmaxMot = AverageRoundsParam(maxRdMotFiles,6)
#    _, avgRdmaxVelo, stdevRdmaxVelo = AverageRoundsParam(maxRdVeloFiles,6)
#    
#    _, avgRdTmaxAng, stdevRdTmaxAng = AverageRoundsParam(TmaxRdAngFiles,6)
#    _, avgRdTmaxCurv, stdevRdTmaxCurv = AverageRoundsParam(TmaxRdCurvFiles,6)
#    _, avgRdTmaxMot, stdevRdTmaxMot = AverageRoundsParam(TmaxRdMotFiles,6)
#    _, avgRdBmaxVelo, stdevRdBmaxVelo = AverageRoundsParam(TmaxRdVeloFiles,6)
#    
#    _, RFavgRdLat, RFstdevRdLat = AverageRoundsParam(ReflexLatencyFiles,6)
#    _, RFavgRdFreqPI, RFstdevRdFreqPI = AverageRoundsParam(ReflexRespFreqPerIntFiles,6)
#    _, RFavgRdFVelo, RFstdevRdFVelo = AverageRoundsParam(ReflexFirstVeloPerIntFiles,6)
#    _, RFavgRdmaxAng, RFstdevRdmaxAng = AverageRoundsParam(ReflexmaxRdAngFiles,6)
#    _, RFavgRdmaxCurv, RFstdevRdmaxCurv = AverageRoundsParam(ReflexmaxRdCurvFiles,6)
#    _, RFavgRdmaxMot, RFstdevRdmaxMot = AverageRoundsParam(ReflexmaxRdMotFiles,6)
#    _, RFavgRdTmaxAng, RFstdevRdTmaxAng = AverageRoundsParam(ReflexTmaxRdAngFiles,6)
#    _, RFavgRdTmaxCurv, RFstdevRdTmaxCurv = AverageRoundsParam(ReflexTmaxRdCurvFiles,6)
#    _, RFavgRdTmaxMot, RFstdevRdTmaxMot = AverageRoundsParam(ReflexTmaxRdMotFiles,6)
#            
#            
##            NegSdAng,PosSdAng = GetSdTraces(avgRdAng,stdevRdAng)
##            NegSdCurv,PosSdCurv = GetSdTraces(avgRdCurv,stdevRdCurv)
##            NegSdMot,PosSdMot = GetSdTraces(avgRdMot,stdevRdMot)
            
#        csvYFiles=glob.glob(AnalysisFolder+'*centre_y.csv')
#        stimFiles=glob.glob(folder+'/*.csv')

    
