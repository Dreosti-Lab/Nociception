# Nociception
This repository contains all the Bonsai and Python scriots that Elisa Clemente is using to analyse her exciting nociception experiments.
The experiments consists in zebrafish embedded in agarose and stimulated with an infra-red laser. 

# Bonsai Scripts
Name: ...
Aim: Script used to collect zebrafish nociception data. 
What does it do: the Bonsai script does the following
What are the files generated:


# Python Analysis

Script 1 : S1_tail_fit_curvature_ARK_EC_extracting_data_final Aim: Track the tail of hed-fixed zebrafish What does it do: Fits fixed-length segments to zebrafish tail and extracts coordinates of edges of segments What are the files generated: 2 csv files with the x/y coordinates of the edges of the tail segments for each frame.

Script 2: S2_tail_parameters_with_intensities_ARK_EC_TR_new Aim: Measure different tail parameters What does it do: It uses the x/y coordinates from script 1 to calculate tail motion, curvature and cumulative angles What are the files generated: 3 csv files with each of the above paramters and 3 csv files with the raw data "cut" along specific time windows (around the laser stimulations).

Script 3: S3_average_rounds_plot_max_latest Aim: Average data from csv files across trials for each fish, determine peak values of each parameter and plot these What does it do: It averages data from the csv files generated with script 2, across trials for each fish, and plots it. It also determines the peak values of each parameter (tail motion, curvature and cumulative angles) and plots them against stimulus intensities. What are the files generated: csv files with average values for each paramater (and standard deviations); csv files with maximum/peak values for each paramater (and standard deviations); png files with plots of the csv files.

Script 4: S4_average_fish_plot_max_latest Aim: As for script 3, but across fish What does it do: As for script 3, but across fish What are the files generated: As for script 3, but across fish.
