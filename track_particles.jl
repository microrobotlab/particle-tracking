using BlobTracking, Images, VideoIO, ImageView, FileIO
using CSV, DataFrames
include("save_data.jl")


##INSERT --- opens the video, creates a iterable stack of frames stored in "vid"
folder="J1\\"
filename="J1brown_27_3march"
#folder="J10+H2O2\\"
#filename="J10H2O2_m58"
#-------------------------------------------------------------------------------

pathV = "..\\tracking_videos\\Results\\" *folder #Path dal PC nuovo, per vecchio sostituisci g.petrucci con petru
pathTOT=pathV*filename*".avi"
io   = VideoIO.open(pathTOT)
vid  = VideoIO.openvideo(io)
img  = first(vid)

#creates a blob tracker with the desired parameters.
mask=trues(1530,2040)
mask[1300:1530,1700:2040].=false
bt = BlobTracker(5:11, #array of blob sizes we want to detect --> era 5 e 11
                3.0, # σw Dynamics noise std. (kalman filter param)  --> era 3.0
                10.0,  # σe Measurement noise std. (pixels) (kalman filter param) --> era 10.0
                mask=mask, #image processing before the detection, not implemented here because unecessary
                #preprocessor = preprocessor, #image processing before the detection, not implemented here because unecessary
                amplitude_th = 0.008, ## with less, like 0.007, it detects false positives (in the Hirox videos) --> era 0.008
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

path="C:\\Users\\g.petrucci\\OneDrive - Scuola Superiore Sant'Anna\\tracking_code\\Results\\"*folder
save(path*"tracking_"*filename*".png", drawimg)
#-----> if we just need the coordinates whitout tracking, use this
#coords = get_coordinates(bt, vid)


#saves data in a dataframe in .csv file. 4 columns: blob ID, time, x and y for each frame.
#framerate is the frame rate of the video
#-----> WRITE the ACTUAL framerate as second entry
resultfilename=path*"coordinates_"*filename*".csv"
save_data(result,20,resultfilename)

