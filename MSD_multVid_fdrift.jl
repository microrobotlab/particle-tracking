using CSV, Pkg, StatsPlots, DataFrames, CategoricalArrays, Plots, NaNStatistics, LsqFit, CurveFit, Statistics, JSON3, FileIO
gr()    #backend dei plot, cerca figure interattive
include("drift_corr.jl")

##---IT NEEDS TO HAVE ONLY THE VIDEOS IN THE ORIGIN FOLDER; AND ALL OF THEM SHOULD BE AT THE SAME MAGNIFICATION

pathORIG="C:\\Users\\g.petrucci\\Scuola Superiore Sant'Anna\\Microscale Robotics Laboratory - RESEARCH - Research\\Data\\NIK_Nikon-phase-contrast\\P01\\20230320_NIK_P01_E009-J27_Pt@PS_GP\\"

#---FOLDER IN WHICH THE .csv ARE STORED -------------
folder="PS_1.1_um\\Pt\\20230320_P01_E3009 - J27\\"
#----------------------------------------------------

#--- made-up INFO -----------------------------------
diamPart=1.1  # in microns
a=1.5                   #for lens 1x or 1.5x , Nikon
um_px = (10/82.978)/a   # micron / pixel : NIK 20x
#um_px = 50/318 # HRX mid 1000x #1/6.32
#um_px = 50/255 # HRX mid 800x
#um_px = 100/382 # HRX mid 600x #
framerate = 25          # fps
#----------------------------------------------------

#definition to be added in the name
boxtrack=20
YlimMSD=10.1
lYlimMSD=-0.1 

path="Results\\"*folder

list=readdir(pathORIG)
for i in list[1:end]
    filename=i
    
    ## Read the data file and save it to a dataframe
    df = CSV.read(path*"coordinates_"*filename*".csv", DataFrame)
    #--- operation on the DataFrame:
    df[!,:BlobID] = convert.(Int64,df[!,:BlobID]);
    df[!,:Frame] = round.(Int64,df[:,:Time]./df[1,:Time])
    df[!,:x] = df[!,:x]*um_px
    df[!,:y] = df[!,:y]*um_px

    #---on the basis of the entry, calculate:
    D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
    Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
    tr=(Dr)^(-1)

    #---Apply the drift correction through the function,
    #---Return gdf_clean_corrected, immobile_tracks, jump_tracks, short_tracks, discard_tracks ##REGISTERED IN OUTPUTS???
    gdf_clean_corrected, immobile_tracks, jump_tracks, short_tracks, discard_tracks = drift_corr(df,um_px,framerate,filename)

    # Calculate number of detected traks
    nTraks=length(gdf_clean_corrected)
    # Find Max length of time vector
    lt = maximum([size(g)[1] for g in gdf_clean_corrected])
    tauMax=ceil(Int,lt/10) #ceil approssima per eccesso, floor per difetto 


    # Plots a restricted number of track, scaled to zero ---> you may have to change axes limits!!!

    #---> Initialize 0TRAJECTORY Plot
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

    function MSDfun(track,tauMax)
        ltrack=length(track[!,:Time])
        tMax=ceil(Int, ltrack/10)    #5 that was 10, is it rigorous or related to something..? #taken from the article Wei Wang(bibbia MSD)
        msd=fill(NaN, tauMax+1)
        msd[1:tMax+1].=0
        for tau in 1:tMax 
            for i in tau+1:ltrack
                msd[tau+1]+=((track[i,:x]-track[i-tau,:x])^2+(track[i,:y]-track[i-tau,:y])^2)/(ltrack-tau)
            end
        end
    #    println(length(msd))
        return msd

    end


        # Plots the MSD   ---> you may have to change axes limits!!!

        #---> Initialize Plot singoli MSD
        graphsingMSD=plot();

        matrMSD=fill(NaN, tauMax+1, length(idx))

        xMSD=Array(0:1/framerate:tauMax/framerate)

        for i in 1:length(idx)
            matrMSD[1:tauMax+1, i] = MSDfun(gdf_clean_corrected[idx[i]],tauMax)
            #    plot!(graphMSD, xMSD, matrMSD[1:tauMax+1, i],label = string("MSD ",i),legend=false)         
        end

        MSD=vec(nanmean(matrMSD, dims=2))
        dsMSD=vec(nanstd(matrMSD; dims=2))
        plot!(xMSD,matrMSD, ylims=(lYlimMSD,YlimMSD), legend=false)  #plotta i singoli MSD di ogni traccia
        plot!(xMSD,MSD, yerror=dsMSD, ylims=(lYlimMSD,YlimMSD), marker=true,legend=false);  #plotta la media
        #Layout(xaxis_range=[0, 1], yaxis_range=[0,2])
        xlabel!("Δt [s]");
        ylabel!("MSD [μm²]")

        display(graphsingMSD)



        #---> Initialize Plot MEDIA MSD
        graphMSD=plot();
        plot!(xMSD,MSD, yerror=dsMSD, ylims=(lYlimMSD,YlimMSD), marker=true,legend=false);  #plotta la media
        #Layout(xaxis_range=[0, 1], yaxis_range=[0,2])
        xlabel!("Δt [s]");
        ylabel!("MSD [μm²]")

        display(graphMSD)


        #---SAVE WITHOUT the fit-------------------------------
        #tauM="tauM"*string(tauMax, base = 10, pad = 2)  #per inserire tauMax nel titolo, sostituito con datetime
        #Dates.format(DateTime, "e, dd u yyyy HH:MM:SS")

        png(graphsingMSD, path*"singMSD_"*filename)
        png(graphMSD, path*"MSD_"*filename)
        png(graphSDtrck, path*"tracks_"*filename)

        #---Save a .csv with the MSD to overlay plots in a second moment
        MSD_df=DataFrame(xMSD=xMSD, MSD=MSD, yerror=dsMSD)
        CSV.write(path*"MSD_"*filename*".csv", MSD_df)

        #---Save variables------------------------------------
        d=Dict("length_idx"=>length(idx), "tauMax"=>tauMax,"nTracks"=>nTraks, "um_px"=>um_px, "framerate"=>framerate, "diamPart"=>diamPart,"idx"=>idx,"D"=>D,"Dr"=>Dr,"tr"=>tr)
        JSON3.write(path*"var_"*filename*".json", d)
        #--to read the JSON3 file and get back the variables--
        #d2= JSON3.read(read("file.json", String))
        #for (key,value) in d2
        #        @eval $key = $value
        #end
end