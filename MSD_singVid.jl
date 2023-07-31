using CSV, Pkg, StatsPlots, DataFrames, CategoricalArrays, Plots, NaNStatistics, LsqFit, CurveFit, Statistics, JSON3, FileIO
gr()    #plot's backend
include("drift_corr.jl")
include("MSDfun.jl")

##-- -FOLDER IN WHICH THE .csv ARE STORED -------------
folder="..."   # name of the folder in wich the CSV file with the results of the tracking is stored
filename="..."   # name of the raw video

##--- Made-up INFOs -----------------------------------
diamPart=3  # mean diameter of the particles to be tracked, in microns
a=1.5                   # part of the microns to pixel convertion, as our Nikon microscope can work with an additional 1.5x lens
um_px = (10/82.978)/a   # micron to pixel convertion for Nikon microscope using the 20x objective: COMMENT ALL BUT THE ONE IN USE!
#um_px = 50/318          # micron to pixel convertion for Hirox microscope using the mid objective with 1000x lens (1/6.32)
#um_px = 50/255          # micron to pixel convertion for Hirox microscope using the mid objective with 1000x lens (1/6.32)# HRX mid 800x
#um_px = 100/382         # micron to pixel convertion for Hirox microscope using the mid objective with 1000x lens (1/6.32)# HRX mid 600x #
framerate = 25          # fps of the video in analysis


##--- Read the data file and save it to a dataframe ---
path="Results\\"*folder
df = CSV.read(path*"coordinates_"*filename*".csv", DataFrame)
## --- Operation on the DataFrame: --------------------
df[!,:BlobID] = convert.(Int64,df[!,:BlobID]);
df[!,:Frame] = round.(Int64,df[:,:Time]./df[1,:Time])
df[!,:x] = df[!,:x]*um_px
df[!,:y] = df[!,:y]*um_px


##--- Rendering of the plots -------------------------
boxtrack=20    # Max X & Y in the box plots 
YlimMSD=10.    # Y Max in MSD plot
lYlimMSD=-0.1  # Y min in MSD plot

##---- On the basis of the entry, calculate: ---------
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     # Diffusive coefficient (um^2/s)
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)       # Rotational diffusive coefficient
tr=(Dr)^(-1)                                              # Rotational time scale

##--- Apply the drift correction through the function "drift_corr.jl" ---
##--- Return gdf_clean_corrected, immobile_tracks, jump_tracks, short_tracks, discard_tracks ---
gdf_clean_corrected, immobile_tracks, jump_tracks, short_tracks, discard_tracks = drift_corr(df,um_px,framerate,filename)

##--- Calculate number of detected traks -------------
nTraks=length(gdf_clean_corrected)
##--- Find Max length of time vector -----------------
lt = maximum([size(g)[1] for g in gdf_clean_corrected])
tauMax=ceil(Int,lt/10) # Decide the time interval (ceil is the excess approximation, floor for defect one)


##--- Plots a restricted number of track, traslated to zero ---> you may have to change axes limits!!!

##--- Initialize 0TRAJECTORY Plot -------------------
graphSDtrck =plot(); 

idx = []
for i in 1:length(gdf_clean_corrected)
    push!(idx,i)
end

for i in rand(idx,min(length(idx),10))

    ## Define the zero position
    x0i= gdf_clean_corrected[i][1,:x]
    y0i = gdf_clean_corrected[i][1,:y]
    ## Add column and fill it with data
    plot!(graphSDtrck, gdf_clean_corrected[i][!,:x].-x0i,gdf_clean_corrected[i][!,:y].-y0i, xlims=(-boxtrack,boxtrack), ylims=(-boxtrack,boxtrack),legend=false,aspect_ratio = 1,framestyle = :box)         
end

display(graphSDtrck)


##--- Calculates & Plots the MSD   ---> you may have to change axes limits!!!

##--- Initialize Plot single MSDs -----------------
graphsingMSD=plot()
matrMSD=fill(NaN, tauMax+1, length(idx))
xMSD=Array(0:1/framerate:tauMax/framerate)


##--- Calculate the MSD through the function "MSDfun.jl" ---
##--- Return MSD -----------------------------------

for i in 1:length(idx)-7
    matrMSD[1:tauMax+1, i] = MSDfun(gdf_clean_corrected[idx[i]],tauMax)
end

MSD=vec(nanmean(matrMSD, dims=2))
dsMSD=vec(nanstd(matrMSD; dims=2))
plot!(xMSD,matrMSD, ylims=(lYlimMSD,YlimMSD), legend=false)    # plots of the single MSDs of each track
plot!(xMSD,MSD, yerror=dsMSD, ylims=(lYlimMSD,YlimMSD), marker=true,legend=false);    # plots the MSD over the single traks
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
display(graphsingMSD)



##--- Initialize Plot MEDIA MSD ------------------
graphMSD=plot();
plot!(xMSD,MSD, yerror=dsMSD, ylims=(lYlimMSD,YlimMSD), marker=true,legend=false);  # plots the MSD alone
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
display(graphMSD)


##--- SAVE WITHOUT the fit --> this is done in a different script ---
png(graphsingMSD, path*"singMSD_PROVA"*filename)
png(graphMSD, path*"MSD_PROVA"*filename)
png(graphSDtrck, path*"tracks_PROVA"*filename)

##--- Save a .csv with the MSD to overlay plots in a second moment ---
MSD_df=DataFrame(xMSD=xMSD, MSD=MSD, yerror=dsMSD)
CSV.write(path*"MSD_PROVA"*filename*".csv", MSD_df)

##--- Save variables --------------------------------
d=Dict("length_idx"=>length(idx), "tauMax"=>tauMax,"nTracks"=>nTraks,"um_px"=>um_px, "framerate"=>framerate, "diamPart"=>diamPart,"idx"=>idx,"D"=>D,"Dr"=>Dr,"tr"=>tr)
JSON3.write(path*"var_PROVA"*filename*".json", d)
#--to read the JSON3 file and get back the variables--
#d2= JSON3.read(read("file.json", String))
#for (key,value) in d2
#        @eval $key = $value
#end