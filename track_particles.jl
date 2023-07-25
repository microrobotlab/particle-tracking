using BlobTracking, Images, VideoIO, ImageView, FileIO
using CSV, DataFrames
include("save_data.jl")

##INSERT --- opens the video, creates a iterable stack of frames stored in "vid"
filename="20230523_NIK_P04_E017_08_GP_06-VID00057" # NOTES: poi 2, poi capire preprocessor..

pathORIG="C:\\Users\\g.petrucci\\Scuola Superiore Sant'Anna\\Microscale Robotics Laboratory - RESEARCH - Research\\Data\\NIK_Nikon-phase-contrast\\P04\\2023\\20230523_NIK_P04_E017_GP\\"
folderDEST="20230523_NIK_P04_E017_GP - J28\\" 
pathDEST="C:\\Users\\g.petrucci\\OneDrive - Scuola Superiore Sant'Anna\\tracking_code\\Results\\3_um\\Pt\\GUVs"*folderDEST

#-------------------------------------------------------------------------------

pathTOT=pathORIG*filename*".avi"
io   = VideoIO.open(pathTOT)
vid  = VideoIO.openvideo(io)
img  = first(vid)

#creates a blob tracker with the desired parameters.

#----- For Hirox ---- WITH NIKON MASK IS NOT REQUIRED; COMMENT IT INSIDE bt!!
#--- AND change blobtracker from 5:11 for Hirox 800x & 1000x to
#mask=trues(1530,2040)
#mask[1300:1530,1700:2040].=false

#---- For NIKON :-----
function preprocessor(storage, img)
    storage .= Float32.(img)
    #update!(medbg, storage) # update the background model
    storage .= abs.(1 .- img)  # You can save some computation by not calculating a new background image every sample
end

bt = BlobTracker(6:8, #array of blob sizes we want to detect --> era 5 e 11 per Hirox, tr2 H2O2: 3:1, tr3 cambia noise2 to 15 tr4 noise2 to 20. Per Nikon alla fine 5:6
                3.0, # σw Dynamics noise std. (kalman filter param)  --> era 3.0
                10.0,  # σe Measurement noise std. (pixels) (kalman filter param) --> Per Hirox era 10.0, ALZA: Portato a 15.0 per Nikon.
#                mask=mask, #image processing before the detection, not implemented here because unecessary
                preprocessor = preprocessor, #image processing before the detection, not implemented here because unecessary
                amplitude_th = 0.008, ## with less, like 0.007, it detects false positives (in the Hirox videos) --> era 0.008. Mantenuto per Nikon
                correspondence = HungarianCorrespondence(p=0.5, dist_th=4), # dist_th is the number of sigmas away from a predicted location a measurement is accepted.--> era p=0.5, dist_th=4
)

#tune_size can be used to automatically tune the size array in bt based on img (the first img of vid). not mandatory.
#tune_sizes(bt, img)


result = track_blobs(bt, vid,
                        display = Base.display, # use nothing to omit displaying.
                        recorder = Recorder(),) # records result to video on disk


#plots trajectories and start-end points for each blob

traces = trace(result, minlife=15) # Filter minimum lifetime of 5
measurement_traces = tracem(result, minlife=5)
drawimg = RGB.(img)
draw!(drawimg, traces, c=RGB(0,0,0.5))
draw!(drawimg, measurement_traces, c=RGB(0.5,0,0))

save(pathDEST*"tracking_"*filename*".png", drawimg)
#save(path*"tracking_"*filename*".svg", drawimg)
#-----> if we just need the coordinates whitout tracking, use this
#coords = get_coordinates(bt, vid)


#saves data in a dataframe in .csv file. 4 columns: blob ID, time, x and y for each frame.
#framerate is the frame rate of the video
#-----> WRITE the ACTUAL framerate as second entry
resultfilename=pathDEST*"coordinates_"*filename*".csv"
save_data(result,20,resultfilename) #the secon entry is the framerate, change it iif you want to have the proper time in the excel file

