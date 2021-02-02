# -*- coding: utf-8 -*-
"""
Fit fixed-length segments to zebrafish tail

@author: Adam Kampff
"""
import os
import cv2
import numpy as np
import pandas as pd 
# import matplotlib.pyplot as plt
import timeit
import glob

### Helper Functions

# Batch processing
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


# Get intensity with subpixel accuracy using bilinear interpolation
def subpixel_intensity(image, x, y):
    xi = np.int(np.round(x))
    yi = np.int(np.round(y))
    dx = x - xi
    dy = y - yi
    weight_tl = (1.0 - dx) * (1.0 - dy)
    weight_tr = (dx)       * (1.0 - dy)
    weight_bl = (1.0 - dx) * (dy)
    weight_br = (dx)       * (dy)
    accum = 0.0
    accum += weight_tl * image[yi, xi]
    accum += weight_tr * image[yi, xi+1]
    accum += weight_bl * image[yi+1, xi]
    accum += weight_br * image[yi+1, xi+1]
    return accum

# Get intensity profile along arc
def arc_profile(image, x, y, radius, start_theta, stop_theta, step_theta):
    profile = []
    xs = []
    ys = []
    for t in np.arange(start_theta, stop_theta, step_theta):
        nx = x + (radius * np.cos(t))
        ny = y + (radius * np.sin(t))
        value = subpixel_intensity(image, nx, ny)
        xs.append(nx)
        ys.append(ny)
        profile.append(value)    
    return np.array(profile), np.array(xs), np.array(ys)


### Algorithm

tic=timeit.default_timer()


# Specify data file root
if not 'folderListFile' in locals():
    folderListFile = input('Path for folder list txt file:').strip("'")

data_path,folderNames=readFolderList(folderListFile)

# Get csv file with pixel coordinates of tail base
if not 'tailBase_coords_file' in locals():
    tailBase_coords_file = input('Path for tail base coords csv file:').strip("'")
    
tailBase_coords = np.genfromtxt(tailBase_coords_file, delimiter=',')

# Define expt type: 'WT' or 'control' (laser offset) or 'lidocaine'
if not 'exptType' in locals():
    exptType = input('Experiment type: WT, control, lidocaine or AITC? ')

if exptType == 'WT':
    expt = ''
elif exptType == 'control':
    expt = '_control'
elif exptType == 'lidocaine':
    expt = '_lidocaine'
elif exptType == 'AITC':
    expt = '_AITC'
    
# Set whether experiment was paired
if not 'pairing' in locals():
    pairing = input('Was the experiment paired?: Yes / No ')

if pairing == 'No':
    numRounds = 3
elif pairing == 'Yes':
    numRounds = 6

# Loop through folders
for idF, folder in enumerate(folderNames):

    # Create analysis folder within each fish 
    fish_analysis_path = folder + '/Analysis/'
    if not os.path.exists(fish_analysis_path):
        os.mkdir(fish_analysis_path)
    
    for r in range(numRounds):
        video_path = glob.glob(folder + '/round_' + str(r+1) + '*.avi')
        head,tail = os.path.split(video_path[0])
        output_path = head + '/Analysis/' + tail[:-4] + expt + '_tracking.avi'
        
    
        # Specify tracking parameters
        display = True
        save = True
        if r < 3:
            tail_base_x = int(tailBase_coords[idF,1])
            tail_base_y = int(tailBase_coords[idF,2])
        else:
            tail_base_x = int(tailBase_coords[idF,3])
            tail_base_y = int(tailBase_coords[idF,4])
        num_segments = 12
        initial_search_angle = (np.pi / 2.0) # Down
        segment_size = 10.0
        
        # Open video
        video = cv2.VideoCapture(video_path[0])
        frame_width = np.int(video.get(cv2.CAP_PROP_FRAME_WIDTH))
        frame_height = np.int(video.get(cv2.CAP_PROP_FRAME_HEIGHT))
        frame_count = np.int(video.get(cv2.CAP_PROP_FRAME_COUNT))
        
        # Open output video (?)
        if save:
                output = cv2.VideoWriter(output_path, cv2.VideoWriter_fourcc('F','M','P','4'), 60, (frame_width,frame_height))
        
        # Open window (?)
        if display:
            cv2.namedWindow("Display")
        
        # Process video
        #video.set(cv2.CAP_PROP_POS_FRAMES, 1)
        #frame_count = np.int(video.get(cv2.CAP_PROP_FRAME_COUNT))
        
        # Create a matrix with num_frames in x, and num_segments in y
        centre_x_Alltemp=[]
        centre_y_Alltemp=[] 
        
        ####
        for f in range(frame_count):
            ret, frame = video.read()
            grayscale = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            smoothed = cv2.GaussianBlur(grayscale, (7,7), 2)
        
            # Displaying?
            if display:
                ret = cv2.drawMarker(frame, (tail_base_x, tail_base_y), (0,255,255), cv2.MARKER_CROSS, 7, thickness=1, line_type=8)
        
            # Start segment search along arc
            centre_x = tail_base_x
            centre_y = tail_base_y
            search_angle = initial_search_angle
            for s in range(num_segments):
                start_theta = search_angle - (np.pi/2.0)
                stop_theta = search_angle + (np.pi/2.0)
                step_theta = (2.0*np.pi)/360.0
                profile, xs, ys = arc_profile(smoothed, centre_x, centre_y, segment_size, start_theta, stop_theta, step_theta)
                print(s)
                # Search arc profile for valley (next search angle)
                search_angle = start_theta + (np.argmin(profile) * step_theta)
        
                # Update centre point
                centre_x = centre_x + (segment_size * np.cos(search_angle))
                centre_y = centre_y + (segment_size * np.sin(search_angle))
                
                # Put centre_x and centre_y into their corresponding lists
                centre_x_Alltemp.append(centre_x)
                centre_y_Alltemp.append(centre_y)
                
                # Displaying?
                if display:
                    ret = cv2.drawMarker(frame, (np.int(centre_x), np.int(centre_y)), (0,0,255), cv2.MARKER_DIAMOND, 7, thickness=1, line_type=8)
                    for i in range(len(profile)):
                        ret = cv2.drawMarker(frame, (np.int(xs[i]), np.int(ys[i])), (0,255,0), cv2.MARKER_SQUARE, 1, thickness=1, line_type=8)
            
            # Saving?
            if save:
                output.write(frame)
        
            # Displaying?
            if display:
                cv2.imshow("Display", frame)
                retval = cv2.waitKey(1)
                if retval == ord('q'):
                    break
        
        # Convert list into array and reshape it ('break it down')
        centre_x_All=np.asarray(centre_x_Alltemp).reshape((frame_count,num_segments))
        centre_y_All=np.asarray(centre_y_Alltemp).reshape((frame_count,num_segments))
        
        # Save a csv file
        pd.DataFrame(centre_x_All).to_csv(head + '/Analysis/' + tail[:-4] + expt + '_centre_x.csv', header=None, index=None)
        pd.DataFrame(centre_y_All).to_csv(head + '/Analysis/' + tail[:-4] + expt + '_centre_y.csv', header=None, index=None)
        
        
        # Saving?
        if save:
            output.release()
        
        # Displaying?
        if display:
            cv2.destroyAllWindows()
        
        # Close video
        video.release()

# End timer
toc = timeit.default_timer()
time_elapsed_sec = toc-tic
time_elapsed_min = time_elapsed_sec/60
time_elapsed_hr = time_elapsed_min/60

time_elapsed_sec_r = round(time_elapsed_sec)
time_elapsed_min_r = round(time_elapsed_min)
time_elapsed_hr_r = round(time_elapsed_hr,1)

print('Finished. Time elapsed:' + str(time_elapsed_sec_r) + ' sec, i.e. ~' + str(time_elapsed_min_r) + ' mins, i.e. ~' + str(time_elapsed_hr_r) + 'h')

#FIN