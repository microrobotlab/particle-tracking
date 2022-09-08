using CSV, Pkg, StatsPlots, DataFrames, CategoricalArrays, Plots, NaNStatistics, LsqFit, CurveFit, Statistics, JSON3, FileIO, Dates
gr()    #backend dei plot, cerca figure interattive

##---INNSERT---same as track_particles------------------
#folder="J11_PS_H2O2esaurita\\"
#filename="J11Brown_m44"

#folder="J11_PS+H2O2\\"
#filename="J11H2O2_m42"

#folder="J10_in milliQ\\"
#filename="movie052_J10_milliQ"

#folder="J10+H2O2\\"
#filename="J10H2O2_m55"

#folder="J1\\"
#filename="J1H2O2_32_3march"
#filename="J1+H2O2_m40_2march"

folder="Sio2inMilliQ\\"
filename="SiO2millliQ_frame_m43"

#filename="movie057_J10+H2O2"


#----- INSERT FROM FIJI -------
#convFact=50/318 #mid 1000x #1/6.32
convFact= 50/255 # mid 800x
diamPart=1  # in microns

#---parameters for the filtering---
framerate=12
ltrackmin=framerate*5 #tauMax> 1 sec 
displminpx=0.5 #medium displacement between two adiacents point > 1 pixel - filters out still particles (rumore su part che stanno ferme)
jump=2 #max jump allowed between 2 frames
#------------------------------
YlimMSD=8.1

## Read the data file and save it to a dataframe
path="Results\\"*folder
df = CSV.read(path*"coordinates_"*filename*".csv", DataFrame)

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
idx=[]   # save IDnumber of good traks
meandr = [mean(sqrt.((diff(g[!,:x])).^2+(diff(g[!,:y]).^2))) for g in gdf]
quartiles = quantile(meandr,[0.25, 0.75])
IQR=quartiles[2]-quartiles[1] #diff(y)
lowFen=quartiles[1]-1.5*IQR
upFen=quartiles[2]+1.5*IQR
boxplot(meandr)
for i in 1:length(gdf)
    flag=false
    ltrack=length(gdf[i][!,:x])
    dr=sqrt.((diff(gdf[i][!,:x])).^2+(diff(gdf[i][!,:y]).^2))  #vector with the instant dr of each track
    if (ltrack>ltrackmin) && (meandr[i]>lowFen && meandr[i]<upFen) #track length> 1 sec + medium displacement between two adiacents point > 1 pixel - filters out still particles (rumore su part che stanno ferme)
        for n in 1:(ltrack-1)
            if (dr[n]>jump*diamPart) # || (gdf[i][n,:x]>1700*convFact && gdf[i][n,:y]>1420*convFact) #max jump allowed between 2 frames x diamPart + cut scalebar (cutted a priori)
                flag=true
                break
            end
        end  
        if flag==false
            push!(idx, i)
        end
    end
end

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
    tMax=min(tauMax, ceil(Int, ltrack/10))
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


## ---Fit----------------------------------------------

if  (0.5*tr*framerate)>5  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi
    # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
    model(t,p)=4*p[1].*t.+(p[2]).*(t).^2
    p0=[D,0.1] #first guess
    fit2=LsqFit.curve_fit(model,xMSD[1:5],MSD[1:5],p0,lower=[0.2*D,0.0],upper=[5*D,10])
    p=fit2.param
    plot!(xMSD,model(xMSD,p), label="Fit")
    #confidence_inter = confidence_interval(fit2, 0.05)
    velox=string(round(sqrt(p[2]),digits=2))
    print("v bal = ")
    println(velox)
    title!("v= "*velox)

elseif tauMax>(2*tr*framerate)+5  #tr<<tauMax, DIFFUSIVE regime, linear fit, ma abbiamo video troppo corti, se riesci a farli di circa 1 min moltiplica per 5 invece
    pfit=linear_fit(xMSD[end-4:end],MSD[end-4:end])
    yfit=pfit[2].*xMSD.+pfit[1]
    plot!(xMSD,yfit)
    Deff=pfit[2]/4
    velox=sqrt((pfit[2]-4D)*tr)
    velox=sqrt(-2pfit[1]/tr^2)

    #    fit=curve_fit(LinearFit, xMSD[end-4:end],MSD[end-4:end])
#    model2(t,p)=4*p[1]*D.*t.-p[2]
#    model(t,p)=4*(p[1]*D+(1/4)*(p[2])*tr)*t.-(p[2])*(((tr)^2)/2)
#    p0=[1.0,0.5] #first guess
#    fit2=LsqFit.curve_fit(model,xMSD[floor(Int,(2*tr*framerate)):end],MSD[floor(Int,(2*tr*framerate)):end],p0)
#    fit2=LsqFit.curve_fit(model2,xMSD[end-4:end],MSD[end-4:end],p0, lower=[1.0,0.0],upper=[10.0,Inf])
#    p=fit2.param
#    plot!(xMSD,model2(xMSD,p), label="Fit")
#    confidence_inter = confidence_interval(fit2, 0.05)
#    velox=string(round(sqrt(p[2]),digits=2))
#    velox=string(round(sqrt(4(p[1]*D-D)/tr),digits=2))
    print("v diff = ")
    println(velox)
    title!("v= "*velox)
end 

display(graphMSD)

png(graphSDtrck,  path*"tracks_"*filename*DateTime)
png(graphSDtrck_dc,  path*"tracks_"*filename*"_dc"*DateTime)
png(graphMSD, path*"MSDfit_"*filename*DateTime)

#---Save a .csv with the MSD to overlay plots in a second moment
MSD_df=DataFrame(xMSD=xMSD, MSD=MSD, yerror=dsMSD)
CSV.write(path*"MSD_"*filename*DateTime*".csv", MSD_df)

#---Save variables------------------------------------
d=Dict("length_idx"=>length(idx), "velox"=>velox,"tauMax"=>tauMax,"nTracks"=>nTraks, "ltrackmin"=>ltrackmin," displminpx"=>displminpx, "jump"=>jump, "convFact"=>convFact, "framerate"=>framerate, "diamPart"=>diamPart,"idx"=>idx)
JSON3.write(path*"var_"*filename*DateTime*".json", d)
#--to read the JSON3 file and get back the variables--
#d2= JSON3.read(read("file.json", String))
#for (key,value) in d2
#        @eval $key = $value
#end
