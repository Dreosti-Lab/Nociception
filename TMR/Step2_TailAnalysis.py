# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 11:45:15 2021

@author: Tom Ryan, Adam Kampff, Elisa C  
Compiled and edited, Tom Ryan
"""
############ Set default inputs and options DO NOT CHANGE ##########
folderListFile=[]
exptType=[]
dose=[]
#####################################################
# Set "Library Path" - Nociceptive Zebrafish Repo
lib_path='C:/Users/thoma/OneDrive/Documents/GitHub/Nociceptive_Zebrafish/libs'
import sys
sys.path.append(lib_path)
# Generic Libraries
import glob
import numpy as np
import timeit
# Custom libraries
import NZ_tailAnalysis as NZA
import NZ_utilities as NZU
import NZ_dictionary as NZD
################################## INPUTS #####################################
############ Directory info and options ##############
root='D:/Elisa/'
folderListFileN='PairedLido5mg.txt'
folderListFile=root+folderListFileN
exptType='lidocaine'
dose='5mgL'
paired=True
saveCSV=True
saveDict=True
ServerExchange=True
if saveCSV and saveDict: print('saving as both csv files and dictionaries')
if saveCSV and not saveDict: print('saving as csv files')
if saveDict and not saveCSV: print('saving as dic')
#################### Acquisition and evoked bout info ########################
# expt=expt+'_TR' # optional hack to add your own suffix to analysis run
FPS=400
noiseK=2 # the lower standard deviation range bound for accepted background noise detection window # I don't like this, see email comments

#################### Evoked bout extraction info ############################
preStimWindow=1 # the time (s) before each stimulus to capture for evoked bout extraction
responseWindow = 2 # the time (s) within which a response is counted as synched to laser (thus 'response')
postStimWindow=2 # the time (s) after each stimulus to capture for evoked bout extraction
reflexWindow = 100 # in ms, time before which a response is considered 'reflex' # resp < 100 ms are reflex(?) -> find 2nd resp too
boutK = 10 # accepted deviation from avg to be considered movement, thus 'response'
l = 2 # accepted deviation from avg for baseline return from movement, thus bout ended
interbout_intervT = 100 # min gap (in ms) between bouts
num_stim = 6

######### Start timer ##########################
tic=timeit.default_timer()

######## Check and define input dependents RARELY IF EVER CHANGE ##############
folderListFile,expt,numRounds,paired=NZU.checkInputs(folderListFile=folderListFile,exptType=exptType,dose=dose,paired=paired)

##### Start run ########
data_path,folderNames=NZU.readFolderList(folderListFile)
    
########### Start of Folder Loop (idF) ############
for idF, folder in enumerate(folderNames):
    # find files created by S1: tracking and stimulus
    AnalysisFolder=folder+'\\Analysis\\'
    chckFiles=glob.glob(AnalysisFolder+'\*.csv')
    stimFiles=glob.glob(folder+'\*.csv')
    csvXFiles,csvYFiles=[],[]
    for chck in chckFiles:
        if chck[-5]=='x':
            csvXFiles.append(chck)
        elif chck[-5]=='y':
            csvYFiles.append(chck)

    ########## Start of file (round) Loop (i) ###############
    for i,csvXFile in enumerate(csvXFiles):
        csvYFile=csvYFiles[i]
        stimFile=stimFiles[i]
        if ServerExchange:
            print('Server Exchange ON: Downloading data file from server...')
            [csvXFile_serv,csvYFile_serv,stimFile_serv]=[csvXFile,csvYFile,stimFile]
            [csvXFile,csvYFile,stimFile]=NZU.segregateSingleAnalysisFiles([csvXFile,csvYFile,stimFile],destFolder=root+'Data/')

        ###### Measure tail motion, angles, and curvature ######
        cumulAngles,curvatures,motion=NZA.computeTailCurvatures(csvXFile,csvYFile)
        
        # Get volt data
        volt_data = np.genfromtxt(stimFile, delimiter=',')
        
        if isinstance(volt_data[0], float) == False: #normal volts file
            volt_data=volt_data[:,0]
            
        elif isinstance(volt_data[0], float) == True: # volts file when I made a mistake on Bonsai ##
            volt_data = np.genfromtxt(stimFile, delimiter=',',dtype=str)
            myVolts = []
            for v in range(len(volt_data)):
                #string = str(volt_data[v])
                splits = volt_data[v].rsplit('.',1)
                splitx = float(splits[0])
                myVolts.append(splitx)
            volt_data = np.asarray(myVolts)

                
        avg_BL,SD_BL = NZA.findStartAndSetBL(volt_data,motion,k=noiseK)
        
        #convert a bunch of stuff from time to frames for bout extraction
        msPerFrame = 1/FPS*1000 # multiply by this to convert frame to ms
        preStimFrames=np.int(np.round((preStimWindow*1000)/msPerFrame))
        postStimFrames=np.int(np.round((postStimWindow*1000)/msPerFrame))
        responseFrames=np.int(np.round((responseWindow*1000)/msPerFrame))
        reflexFrames=np.int(np.round((reflexWindow)/msPerFrame))
        interbout_interv=np.int(np.round(interbout_intervT/msPerFrame))
    
        ### Cut data ### by volt epochs, add volt to dataframe, sort by volt
        cumulAnglesCut = NZA.findStartEndAndCutEC(volt_data,cumulAngles,preStimFrames=preStimFrames,postStimFrames=postStimFrames)
        curvaturesCut = NZA.findStartEndAndCutEC(volt_data,curvatures,preStimFrames=preStimFrames,postStimFrames=postStimFrames)
        motionCut = NZA.findStartEndAndCutEC(volt_data,motion,preStimFrames=preStimFrames,postStimFrames=postStimFrames)
        
        # Set thresholds for response detection
        threshold_down = avg_BL + l*SD_BL # What is a "return to baseline"?
        threshold_up = avg_BL + boutK*SD_BL # What is a "movement"?
        
        ## Measure bout frequency, latency and velocity ##
        mo=np.array(motionCut)
        ca=np.array(cumulAnglesCut)
        timeResp_PerI,reflex_timeResp_PerI,veryFirst_timeResp_PerI_all,freqResp_PI,reflex_freqResp_PI,NumResp_PerI,velocities,avgVelo_PerI,firstVelo_PerI,reflex_firstVelo_PerI = NZA.computeEvokedBouts(
                mo,
                ca,
                boutK,
                l,
                preStimFrames,
                responseFrames,
                reflexFrames,
                interbout_interv,
                threshold_down,
                threshold_up,
                avg_BL,
                SD_BL,
                msPerFrame,
                num_stim)
        
        ## Make response plots based on previous intensity... look for trends to see if pain has memory
        unsorted_motionCut=NZA.findStartEndAndCutEC(volt_data,motion,preStimFrames=preStimFrames,postStimFrames=postStimFrames,sort=False) # hack to grab original order of voltages
        volts=np.array(unsorted_motionCut['Voltage']) # grab voltages
#        preRespPI=NZA.currentVsPrevResponses(volts,freqResp_PI)
        
        freqResp_Rd = (num_stim-timeResp_PerI.count(None))/num_stim
        freqResp_Rd = [freqResp_Rd] # overall response frequency per round
        
        ###### Find peak cumulAng, curv, motion and velocity ######
        maxRdAng,TMaxRdAng = NZA.MaxTimeDF(cumulAnglesCut,preStimFrames,freqResp_PI)
        maxRdCurv,TMaxRdCurv = NZA.MaxTimeDF(curvaturesCut,preStimFrames,freqResp_PI)
        maxRdMot,TMaxRdMot = NZA.MaxTimeDF(motionCut,preStimFrames,freqResp_PI) 
        maxRdVelo,BoutMaxRdVelo = NZA.MaxTimeParam(velocities,6,freqResp_PI)
        
        ###### Find REFLEX peak cumulAng, curv, motion and velocity ######
        RmaxRdAng,RTMaxRdAng = NZA.MaxTimeReflexDF(cumulAnglesCut,preStimFrames,reflexFrames,reflex_freqResp_PI)
        RmaxRdCurv,RTMaxRdCurv = NZA.MaxTimeReflexDF(curvaturesCut,preStimFrames,reflexFrames,reflex_freqResp_PI)
        RmaxRdMot,RTMaxRdMot = NZA.MaxTimeReflexDF(motionCut,preStimFrames,reflexFrames,reflex_freqResp_PI) 
        
####### Here is good example of good practice in efficiency and making code alter friendly #######
# Don't do the same thing (like save csv files) over and over again.
# Instead, make a function that saves it (if you need/want), then make a list of the data and loop through it
# This is a) less typing therefore less opportunity for mistake
#         b) easier to read, test, understand and verify
#         b) if you want to add or change a name or variable, you only have to add it to the lists
    
        #List variables you want to save
        data1 =     [cumulAnglesCut,
                     curvaturesCut,
                     motionCut,
                     timeResp_PerI,
                     reflex_timeResp_PerI,
                     freqResp_PI,
                     reflex_freqResp_PI,
                     NumResp_PerI,
                     freqResp_Rd,
#                     velocities,
                     avgVelo_PerI,
                     firstVelo_PerI,
                     reflex_firstVelo_PerI]
            
        # List their names
        labels1   = ['CumulAngCut',
                     'CurvCut',
                     'MotCut',
                     'LatencyPerInt',
                     'ReflexLatencyPerInt',
                     'RespFreqPerInt',
                     'ReflexRespFreqPerInt',
                     'BoutsPerInt',
                     'RdRespFreq',
#                     'VelocityPerBout',
                     'AvgVelocityPerInt',
                     'FirstVelocityPerInt',
                     'ReflexFirstVelocityPerInt']
        
        data2 =    [ maxRdAng,
                     maxRdCurv,
                     maxRdMot,
                     maxRdVelo,
                     TMaxRdAng,
                     TMaxRdCurv,
                     TMaxRdMot,
                     BoutMaxRdVelo]
            
        labels2  = ['maxRdAng',
                    'maxRdCurv',
                    'maxRdMot',
                    'maxRdVelo',
                    'TimeTOmaxRdAng',
                    'TimeTOmaxRdCurv',
                    'TimeTOmaxRdMot',
                    'BoutOFmaxRdVelo']
        
        data3 =    [RmaxRdAng,
                    RmaxRdCurv,
                    RmaxRdMot,
                    RTMaxRdAng,
                    RTMaxRdCurv,
                    RTMaxRdMot]
        
        labels3 = ['ReflexmaxRdAng',
                   'ReflexmaxRdCurv',
                   'ReflexmaxRdMot',
                   'ReflexTimeTOmaxRdAng',
                   'ReflexTimeTOmaxRdCurv',
                   'ReflexTimeTOmaxRdMot']
        
        lastlabels=['VeryFirstLatencyPerInt'] # only has one value so does not average like the others
        
        lastData=[veryFirst_timeResp_PerI_all]
        ## Save parameters used this round
        metaLabels = ['expType',
                      'dose',
                      'numStim',
                      'voltOrder',
                      'paired?',
                      'FPS',
                      'noiseK',
                      'preStimWindow',
                      'responseWindow',
                      'postStimWindow',
                      'reflexWindow',
                      'responseK',
                      'returnK',
                      'interBoutInt']
        
        metaValues = [exptType,
                      dose,
                      num_stim,
                      volts,
                      paired,
                      FPS,
                      noiseK,
                      preStimWindow,
                      responseWindow,
                      postStimWindow,
                      reflexWindow,
                      boutK,
                      l,
                      interbout_intervT]
        
        allLabels=labels1+labels2+labels3+lastlabels
        allData=data1+data2+data3+lastData
        
        ## Create a metaData analysis dictionary for this round
        metaDict=NZD.makeDictionary(labels=metaLabels,values=metaValues)
        
        if saveDict: # if saving in a dictionary, create a new directory on the level of fish, and save dictionaries for each round
            wDir,_,name=csvXFile.rsplit(sep='\\',maxsplit=2)
            dictPath=wDir+'\\Dictionaries\\'
            NZU.cycleMkDir(dictPath)
            dictPath=dictPath+name[0:-13]+'_Analysis.npy'
            thisFishDict,dictPath=NZD.createSingleFishDict(allLabels,allData,dictPath,metaDict=metaDict)
            if thisFishDict==-1: 
                print('Creating Dictionary failed... saving csvs instead')
                saveCSV=True
                
        if saveCSV:
            # Create csv paths and save csv files
            round_root = csvXFile[0:-13]
            # save metaData as csv
            pathMeta=round_root+'_meta.csv'
            print('Saving metaData for this round at ' + pathMeta )
            with open(pathMeta, 'w') as f:
        
                for key in metaDict.keys():
                    f.write("%s,%s\n"%(key,metaDict[key]))   
                
            # loop through data and labels, build a path to each in the Analysis folder and save
            for i,label in enumerate(allLabels):
                dat=allData[i]
                thisSuff='_'+label+'.csv'
                savePath=round_root+thisSuff
                NZU.saveCSV(dat,savePath)
                
        
        ### End of Round Loop
    ### end of Folder (Fish) loop

toc = timeit.default_timer()
time_elapsed_sec = toc-tic
time_elapsed_min = time_elapsed_sec/60
time_elapsed_hr = time_elapsed_min/60

time_elapsed_sec_r = round(time_elapsed_sec)
time_elapsed_min_r = round(time_elapsed_min)
time_elapsed_hr_r = round(time_elapsed_hr,1)

print('Finished. Time elapsed:' + str(time_elapsed_sec_r) + ' sec, i.e. ~' + str(time_elapsed_min_r) + ' mins, i.e. ~' + str(time_elapsed_hr_r) + 'h')

# FIN    
