using CSV, Pkg, DataFrames, CategoricalArrays, Plots, NaNStatistics, CurveFit, LsqFit, Statistics, JSON3, FileIO
gr()    #backend dei plot, cerca figure interattive

## Read the data file and save it to a dataframe: choose between one of the two
folder="J10+H2O2\\"
filename="movie057_J10+H2O2"

path="Results\\"*folder
#df = CSV.read("brownian_coordniates.csv", DataFrame)
#df = CSV.read("h2o2_1_coordinates.csv", DataFrame)
df = CSV.read(path*"coordinates_"*filename*".csv", DataFrame)

# Make x and y columns float and convert pixels to microns

#----- INSERT FROM FIJI ------- REMEMBER to change filename and the savings at the end, if running all
#convFact=50/318 #mid 1000x #1/6.32
convFact= 50/255 # mid 800x
diamPart=1  # in microns
#------------------------------
#---on the basis of the entry, calculate tauMax for the parabolic fitting, for other fits tauMax< 1/5 length video
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)
tm=tr*12
tauMax=40 #floor(Int, tm)
#---parameters for the filtering---
ltrackmin=12 #track length> 1 sec 
displminpx=0.5 #medium displacement between two adiacents point > 1 pixel - filters out still particles (rumore su part che stanno ferme)
jump=2 #max jump allowed between 2 frames
#----------------------------------

df[!,:x] = df[!,:x]*convFact
df[!,:y] = df[!,:y]*convFact

## Make one column categorical
df[!,:BlobID] = categorical(df[!,:BlobID],compress=true)
## Group dataframe by values in categorical column
gdf = groupby(df,:BlobID,sort=true)


# filter out trajectories shorter than 1 sec (12 frames), the still ones, those with jumps and cut the scalebar

#---> Initialize Plot

idx=[]   # save IDnumber of good traks

for i in 1:length(gdf)
    flag=false
    ltrack=length(gdf[i][!,:x])
    dr=sqrt.((diff(gdf[i][!,:x])).^2+(diff(gdf[i][!,:y]).^2))  #vector with the instant dr of each track
    if (ltrack>ltrackmin) && (mean(dr)>(displminpx*convFact)) #track length> 1 sec + medium displacement between two adiacents point > 1 pixel - filters out still particles (rumore su part che stanno ferme)
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

print("n. of filtered traks = ")
println(length(idx))

# Calculate number of detected traks
nTraks=length(gdf)
# Find Max length of time vector
lt = maximum([size(g)[1] for g in gdf])

# Filters out the drift due to the addition of liquid and save results on a new column of gdf
matrVx=fill(NaN, lt, nTraks)
matrVx[1,:].=0.0
matrVy=fill(NaN, lt, nTraks)
matrVy[1,:].=0.0

for i in idx
    dx=diff(gdf[i][!,:x])
    dy=diff(gdf[i][!,:y])
    matrVx[2:length(dx)+1,i]=dx*12
    matrVy[2:length(dy)+1,i]=dy*12
end
Vx=vec(nanmean(matrVx, dims=2))
Vy=vec(nanmean(matrVy, dims=2))
Δx=cumsum(Vx/12)
Δy=cumsum(Vy/12)
plot(Δx)
plot!(Δy)

for i in 1:nTraks
    h=length(gdf[i][!,:Time])
    gdf[i][!,:xdc]=gdf[i][!,:x].-Δx[1:h]
    gdf[i][!,:ydc]=gdf[i][!,:y].-Δy[1:h]
end


# Plots a restricted number of track, scaled to zero ---> you may have to change axes limits!!!

#---> Initialize Plot
graphSDtrck =plot();

for i in rand(idx,min(length(idx),10))
     
    ## Define the zero position
    x0i= gdf[i][1,:xdc]
    y0i = gdf[i][1,:ydc]
    ## Add column and fill it with data
    plot!(graphSDtrck, gdf[i][!,:xdc].-x0i,gdf[i][!,:ydc].-y0i, xlims=(-10.0,10.0), ylims=(-10.0,10.0),legend=false,aspect_ratio = 1,framestyle = :box)         
end

display(graphSDtrck)


#---> Initialize Plot
graphSDtrck =plot();

for i in idx
    ## Add column and fill it with data
    plot!(graphSDtrck, gdf[i][!,:xdc],gdf[i][!,:ydc], color=:red, legend=false, aspect_ratio = 1,framestyle = :box)         
    plot!(graphSDtrck, gdf[i][!,:x],gdf[i][!,:y], color=:green, legend=false, aspect_ratio = 1,framestyle = :box)         
end

display(graphSDtrck)


function MSDfun(track,tauMax)
    ltrack=length(track[!,:Time])
    msd=zeros(tauMax+1)
    for tau in 1:tauMax 
        for i in tau+1:ltrack
            msd[tau+1]+=((track[i,:xdc]-track[i-tau,:xdc])^2+(track[i,:ydc]-track[i-tau,:ydc])^2)/(ltrack-tau)
        end
    end
#    println(length(msd))
    return msd

end


# Plots the MSD   ---> you may have to change axes limits!!!

#---> Initialize Plot
graphMSD=plot();

matrMSD=fill(NaN, tauMax+1, length(idx))

xMSD=Array(0:1/12:tauMax/12)

for i in 1:length(idx)   # -1 perché ho scoperto che la traiettoria che andava fuori scala erano l'ultimma, in entrambi i casi (-4 x Brownian)
    matrMSD[1:tauMax+1, i] = MSDfun(gdf[idx[i]],tauMax)
#    plot!(graphMSD, xMSD, matrMSD[1:tauMax+1, i],label = string("MSD ",i),legend=false)         
end

MSD=vec(nanmean(matrMSD, dims=2))
dsMSD=vec(nanstd(matrMSD; dims=2))
plot!(xMSD,matrMSD, ylims=(-0.10,4.10), legend=false)
plot!(xMSD,MSD, yerror=dsMSD, ylims=(-0.10,4.10), marker=true,legend=false);
#Layout(xaxis_range=[0, 1], yaxis_range=[0,2])
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

display(graphMSD)




#---SAVE WITHOUT the fit-------------------------------
tauM="tauM"*string(tauMax, base = 10, pad = 2)
png(graphMSD, path*"MSD_"*filename*tauM)

## ---Fit----------------------------------------------
# model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
model(t,p)=4*p[1].*t.+(p[2].^2).*(t).^2
p0=[D,0.1] #first guess
fit2=LsqFit.curve_fit(model,xMSD,MSD,p0,lower=[0.2*D,0.0],upper=[5*D,10])
p=fit2.param
plot!(xMSD,model(xMSD,p), label="Fit")
#confidence_inter = confidence_interval(fit2, 0.05)
velox=string(round(p[2],digits=2))
print("v = ")
println(velox)
#title!("+ H2O2: v= "*velox)

png(graphSDtrck,  path*"tracks_"*filename*tauM)
png(graphMSD, path*"MSDfit_"*filename*tauM)

#---Save a .csv with the MSD to overlay plots in a second moment


#---Save variables------------------------------------
d=Dict("length_idx"=>length(idx), "velox"=>velox,"tauMax"=>tauMax,"nTracks"=>nTraks, "ltrackmin"=>ltrackmin,"displminpx"=>displminpx, "jump"=>jump, "idx"=>idx)
JSON3.write(path*"var_"*filename*tauM*".json", d)
#--to read the JSON3 file and get back the variables--
#d2= JSON3.read(read("file.json", String))
#for (key,value) in d2
#        @eval $key = $value
#end
