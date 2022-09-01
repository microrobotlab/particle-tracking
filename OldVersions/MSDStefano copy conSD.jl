using CSV, Pkg, DataFrames, CategoricalArrays, Plots, NaNStatistics, CurveFit, LsqFit 
gr()    #backend dei plot, cerca figure interattive

## Read the data file and save it to a dataframe
df = CSV.read("h2o2_1_coordinates.csv", DataFrame)

# Make x and y columns float and convert pixels to microns

#----- INSERT FROM FIJI -------
convFact=1/6.32
diamPart=1  # in microns
#------------------------------

df[!,:x] = df[!,:x]*convFact
df[!,:y] = df[!,:y]*convFact

## Make one column categorical
df[!,:BlobID] = categorical(df[!,:BlobID],compress=true)
## Group dataframe by values in categorical column
gdf = groupby(df,:BlobID,sort=true)

# Calculate number of detected traks
nTraks=length(gdf)

#---> Initialize Plot
graphSDtrck =plot();

for i in 1:length(gdf)
     
    ## Define the zero position
    x0i= gdf[i][1,:x]
    y0i = gdf[i][1,:y]
    ## Add column and fill it with data
    gdf[i][!,:SD] = ((gdf[i][!,:x]).-x0i).^2. +((gdf[i][!,:y]).-y0i).^2
    plot!(graphSDtrck, gdf[i][!,:Time],gdf[i][!,:SD],label = string("SD ",i),legend=:bottomright)         
end

display(graphSDtrck)


# Find Max length of time vector
lt=0
for g in gdf
    if size(g)[1]>lt 
        lt= size(g)[1]
    end
end

# filter out trajectories shorter than 1 sec (12 frames), still ones and those with jumps

#---> Initialize Plot
graphSD=plot();

idx=[]   # IDnumber of good traks
matr=fill(NaN, lt, nTraks)

for i in 1:length(gdf)
    flag=false
    ltrack=length(gdf[i][!,:x])
    if (ltrack>12) && (maximum(gdf[i][!,:SD])>(diamPart))
        for n in 1:(ltrack-1)
            if (abs((gdf[i][n,:SD]-gdf[i][n+1,:SD]))>20) || (gdf[i][n,:x]>1700*convFact && gdf[i][n,:y]>1420*convFact) #max jump allowed between 2 frames = 1um ?
#            if (abs((gdf[i][n,:SD]-gdf[i][n+1,:SD]))>20) #max jump allowed between 2 frames = 1um ?

                flag=true
                break
            end
        end  
        if flag==false
            
            push!(idx, i)
            matr[1:ltrack, i]=gdf[i][!,:SD]
            plot!(gdf[i][!,:Time],gdf[i][!,:SD],legend=false)
        end
    end
end

display(graphSD)


for i in 1:idx
   
    ## Add column and fill it with data
    gdf[i][!,:MSD] = MSDfun(gdf[i])
    plot!(graphSDtrck, gdf[i][!,:Time],gdf[i][!,:SD],label = string("SD ",i),legend=:bottomright)         
end

function MSDfun(track;tauMax=9)
    ltrack=length(track[!,:Time])
    msd=zeros(tauMax+1)
    for tau in 1:tauMax 
        #MSD[tau+1]=0
        for i in tau+1:ltrack
            msd[tau+1]+=((track[i,:x]-track[i-tau,:x])^2+(track[i,:y]-track[i-tau,:y])^2)/(ltrack-tau)
        end
    end
    return msd

end


#display(graphSDtrck)



MSD=vec(nanmean(matr, dims=2))
t=Array(0:(1/12):(lt-1)/12)
plot!(t,MSD, marker=true);
xlabel!("Time [s]");
ylabel!("MSD [μ²]")

## ---Fit----------------------------------------------
# (dL)^2 = 4D*(dt) + V^2*(dt)^2
# model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
model(t,p)=4*p[1].*t.+(p[2].^2).*(t).^2
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)
p0=[D,0.04] #first guess
#fit2=LsqFit.curve_fit(model,t[1:12],MSD[1:12],p0,lower=[0.1*D,0.0],upper=[10*D,Inf])
fit2=LsqFit.curve_fit(model,t,MSD,p0,lower=[0.1*D,0.0],upper=[10*D,Inf])
p=fit2.param
plot!(t,model(t,p), label="Fit")
confidence_inter = confidence_interval(fit2, 0.05)
velox=string(round(p[2],digits=2))
title!("+ H2O2: v= "*velox)
#fit=poly_fit(t,MSD,2)
#plot!(t,fit[1].*(t.^0)+fit[2].*t+fit[3].*(t.^2))

#png(graphSD, "BrownianMotion")
png(graphSD, "+H2O2")