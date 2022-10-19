using CSV, Pkg, StatsPlots, DataFrames, CategoricalArrays, Plots, NaNStatistics, LsqFit, CurveFit, Statistics, JSON3, FileIO, Dates
gr()    #backend dei plot, cerca figure interattive

##---INNSERT---same as track_particles------------------
#folder="J11milliQlong\\"
#filename="movie057"
folder="J9H2O2long\\Diluted\\"
filename="J9_H2O2_tot_"

filename1="m48_dd"  #48_dil non fa #traiettorie brutte, ultima g ha un solo punto, viene NaN meandr, mettendo meandr[1:end-1] al calcolo dei quartili va ma non filtra bene
filename2="m49_dd"
filename3="m50_dd"
filename4="m52_dd"
filename5="m53_dd"
#----- INSERT FROM FIJI -------
convFact=50/318 #mid 1000x #1/6.32
#convFact= 50/255 # mid 800x

diamPart=1  # in microns

#---parameters for the filtering---
framerate=12
ltrackmin=framerate*10 #tauMax> 1 sec 
jump=2 #max jump allowed between 2 frames
#------------------------------
YlimMSD=20.1

## Read the data file and save it to a dataframe
path="Results\\"*folder
df1 = CSV.read(path*"coordinates_"*filename1*".csv", DataFrame)
df1[!,:BlobID]=df1[!,:BlobID].+1000
df2 = CSV.read(path*"coordinates_"*filename2*".csv", DataFrame)
df2[!,:BlobID]=df2[!,:BlobID].+2000
df3 = CSV.read(path*"coordinates_"*filename3*".csv", DataFrame)
df3[!,:BlobID]=df3[!,:BlobID].+3000
df4 = CSV.read(path*"coordinates_"*filename4*".csv", DataFrame)
df4[!,:BlobID]=df4[!,:BlobID].+4000
df5 = CSV.read(path*"coordinates_"*filename5*".csv", DataFrame)
df5[!,:BlobID]=df5[!,:BlobID].+5000

df=vcat(df1,df2,df3,df4,df5,cols=:union)
# Make x and y columns float and convert pixels to microns

#---on the basis of the entry, calculate tauMax for the parabolic fitting, for other fits tauMax< 1/5 (o 1/10) *length video
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)
tm=tr*framerate  #tauMax sensato per fare un fitting parabolico, in base al frame rate

#----------------------------------

df[!,:x] = df[!,:x]*convFact
df[!,:y] = df[!,:y]*convFact

## Make one column categorical
df[!,:BlobID] = categorical(df[!,:BlobID],compress=true)
## Group dataframe by values in categorical column
gdf = groupby(df,:BlobID,sort=true)


# filter out trajectories shorter than 5 sec (12 frames), those with jumps and the outliers

filter=histogram();
idx=[]   # save IDnumber of good traks
meandr = [mean(sqrt.((diff(g[!,:x])).^2+(diff(g[!,:y]).^2))) for g in gdf]
display(histogram!(meandr, bins = 0:0.01:0.3))
meandr[meandr.<convFact*(sqrt(2)/2)].=NaN
quartiles = [nanquantile(meandr,0.25), nanquantile(meandr,0.75)]
IQR=quartiles[2]-quartiles[1] #diff(y)
lowFen=quartiles[1]-1.5*IQR
upFen=quartiles[2]+1.5*IQR
#boxplot(meandr)
display(histogram!(meandr, bins = 0:0.01:0.3))
for i in 1:length(gdf)
    flag=false
    ltrack=length(gdf[i][!,:x])
    dr=sqrt.((diff(gdf[i][!,:x])).^2+(diff(gdf[i][!,:y]).^2))  #vector with the instant dr of each track
    if (ltrack>ltrackmin) && (meandr[i]>lowFen && meandr[i]<upFen) #track length> 1 sec + medium displacement between two adiacents point > 1 pixel - filters out still particles (rumore su part che stanno ferme)
        for n in 1:(ltrack-1)
            if (dr[n]>jump*diamPart) #max jump allowed between 2 frames x diamPart
                flag=true
                break
            end
        end  
        if flag==false
            push!(idx, i)
        end
    end
end

display(histogram!(meandr[idx], bins = 0:0.01:0.3))

print("n. of tracks = ")
println(length(gdf))
print("n. of filtered tracks = ")
println(length(idx))

# Calculate number of detected traks
nTraks=length(gdf)
# Find Max length of time vector
lt = maximum([size(g)[1] for g in gdf])
tauMax=ceil(Int,lt/10) #ceil approssima per eccesso, floor per difetto 

# Filters out the drift due to the addition of liquid and save results on a new column of gdf
matrVx=fill(NaN, lt, nTraks)
matrVx[1,:].=0.0
matrVy=fill(NaN, lt, nTraks)
matrVy[1,:].=0.0

for i in idx
    dx=diff(gdf[i][!,:x])
    dy=diff(gdf[i][!,:y])
    matrVx[2:length(dx)+1,i]=dx*framerate
    matrVy[2:length(dy)+1,i]=dy*framerate
end
Vx=vec(nanmean(matrVx, dims=2))
Vy=vec(nanmean(matrVy, dims=2))
Δx=cumsum(Vx/framerate)
Δy=cumsum(Vy/framerate)
plot(Δx)
display(plot!(Δy))

for i in 1:nTraks
    h=length(gdf[i][!,:Time])
    gdf[i][!,:xdc]=gdf[i][!,:x].-Δx[1:h]
    gdf[i][!,:ydc]=gdf[i][!,:y].-Δy[1:h]
end


# Plots a restricted number of track, scaled to zero ---> you may have to change axes limits!!!

#---> Initialize 0TRAJECTORY Plot
graphSDtrck =plot();

for i in rand(idx,min(length(idx),10))
     
    ## Define the zero position
    x0i= gdf[i][1,:xdc]
    y0i = gdf[i][1,:ydc]
    ## Add column and fill it with data
    plot!(graphSDtrck, gdf[i][!,:xdc].-x0i,gdf[i][!,:ydc].-y0i, xlims=(-10.0,10.0), ylims=(-10.0,10.0),legend=false,aspect_ratio = 1,framestyle = :box)         
end

display(graphSDtrck)


#---> Initialize Plot track and DRIFT corrected
graphSDtrck_dc =plot(yflip=true);

for i in idx
    ## Add column and fill it with data
    plot!(graphSDtrck_dc, gdf[i][!,:xdc],gdf[i][!,:ydc], color=:red, legend=false, aspect_ratio = 1,framestyle = :box)         
    plot!(graphSDtrck_dc, gdf[i][!,:x],gdf[i][!,:y], color=:green, legend=false, aspect_ratio = 1,framestyle = :box)         
end

display(graphSDtrck_dc)


function MSDfun(track,tauMax)
    ltrack=length(track[!,:Time])
    tMax=min(tauMax, ceil(Int, ltrack/5))
    msd=fill(NaN, tauMax+1)
    msd[1:tMax+1].=0
    for tau in 1:tMax 
        for i in tau+1:ltrack
            msd[tau+1]+=((track[i,:xdc]-track[i-tau,:xdc])^2+(track[i,:ydc]-track[i-tau,:ydc])^2)/(ltrack-tau)
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
    matrMSD[1:tauMax+1, i] = MSDfun(gdf[idx[i]],tauMax)
    #    plot!(graphMSD, xMSD, matrMSD[1:tauMax+1, i],label = string("MSD ",i),legend=false)         
end

MSD=vec(nanmean(matrMSD, dims=2))
dsMSD=vec(nanstd(matrMSD; dims=2))
plot!(xMSD,matrMSD, ylims=(-0.10,YlimMSD), legend=false)  #plotta i singoli MSD di ogni traccia
plot!(xMSD,MSD, yerror=dsMSD, ylims=(-0.10,YlimMSD), marker=true,legend=false);  #plotta la media
#Layout(xaxis_range=[0, 1], yaxis_range=[0,2])
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

display(graphsingMSD)



#---> Initialize Plot MEDIA MSD
graphMSD=plot();
plot!(xMSD,MSD, yerror=dsMSD, ylims=(-0.10,YlimMSD), marker=true,legend=false);  #plotta la media
#Layout(xaxis_range=[0, 1], yaxis_range=[0,2])
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

display(graphMSD)


#---SAVE WITHOUT the fit-------------------------------
#tauM="tauM"*string(tauMax, base = 10, pad = 2)  #per inserire tauMax nel titolo, sostituito con datetime
DateTime= Dates.format(now(), "dduyy_HHMM") #Dates.now() #Dates.format(now(), "HH:MM")
#Dates.format(DateTime, "e, dd u yyyy HH:MM:SS")

png(graphsingMSD, path*"singMSD_"*filename*DateTime)
png(graphMSD, path*"MSD_"*filename*DateTime)
png(graphSDtrck, path*"tracks_"*filename*DateTime)
png(graphSDtrck_dc, path*"tracks_dc_"*filename*DateTime)
png(filter, path*"filter_"*filename*DateTime)


#---Save a .csv with the MSD to overlay plots in a second moment
MSD_df=DataFrame(xMSD=xMSD, MSD=MSD, yerror=dsMSD)
CSV.write(path*"MSD_"*filename*DateTime*".csv", MSD_df)

#---Save variables------------------------------------
d=Dict("length_idx"=>length(idx), "tauMax"=>tauMax,"nTracks"=>nTraks, "ltrackmin"=>ltrackmin, "jump"=>jump, "convFact"=>convFact, "framerate"=>framerate, "diamPart"=>diamPart,"idx"=>idx)
JSON3.write(path*"var_"*filename*DateTime*".json", d)
#--to read the JSON3 file and get back the variables--
#d2= JSON3.read(read("file.json", String))
#for (key,value) in d2
#        @eval $key = $value
#end
