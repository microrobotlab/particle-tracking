## Purpose
Tracking videos with circular moving objects, in order to acquire their trajectory as a .csv file containing  X and Y coordinates for each video frame (file: track_particles.jl).
The .csv can then be imported through a second script that calculate the mean square displacement (MSD) and plots a series of useful graph to visualize the singles MSD, the trend of the trajectories and the MSD with its standard deviation (file: MSD_singVid.jl).
A third script is then used to perform either a ballistic or a linear fit contemporary on a data file containing active and inactive particles (file: fitBothMSDs.jl)

## Authors
Gaia Petrucci, Yasmine Abedour, Stefano Palagi

## Affiliation
The BioRobotics Institute, Scuola Superiore Sant’Anna, Pontedera (Pisa), Italy

## Release date
31-July-2023

## Prerequisites
Julia

##  Instructions to run the code
To track a video, open the file track_particles.jl.
Set the file name and the path of the folder containing the video to be tracked, and folder and destination path for the outputs.
In case you have to mask any part of the video, you can select the pixel you want to exclude and the activate the mask inside the function "BlobTracker". Same for the preprocessor. 
As an output you will have a .png file with the trajectories plotted over the first frame of the video, plus a .csv  file containing  X and Y coordinates for each video frame. 

The .csv can then be imported through a second script, MSD_singVid.jl, that calculate the mean square displacement (MSD) and plots a series of graphs useful to visualize the singles MSD, the trend of the trajectories and the MSD with its standard deviation. Also here, set the proper folder and filename, plus the experimental parameters, as written in the comments. 

A third script, fitBothMSDs.jl, is then used to perform either a ballistic or a linear fit contemporary on data files containing active and inactive particles.

## References
A Practical Guide to Analyzing and Reporting the Movement of Nanoscale Swimmers
Wei Wang and Thomas E. Mallouk
ACS Nano 2021 15 (10), 15446-15460
DOI: 10.1021/acsnano.1c07503