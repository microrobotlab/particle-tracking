using CSV, DataFrames, Plots, JSON3, LsqFit, Statistics, Dates
gr()    # backend dei plot
include("model_par.jl")
include("model_lin.jl")

## --- Made-up INFOs ----------------------------------
diamPart=3   # mean diameter of the particles to be tracked, in microns
framerate=25    # fps of the video in analysis

##--- Brownian MSD ------------------------------------
folder1="..."   # name of the folder in wich the CSV file with the results of the tracking of the Brownian inactive particles is stored
filename1="..."   # name of the CSV file 
path1=folder1*filename1
##--- Active MSD --------------------------------------
folder2="..."   # name of the folder in wich the CSV file with the results of the tracking of the active particles is stored
filename2="..."   # name of the CSV file
path2=folder2*filename2

##--- Read the data file and save it to a dataframe ----
dfMSDp = CSV.read(path1*".csv", DataFrame)
dfMSDa = CSV.read(path2*".csv", DataFrame)

##--- Needed parameters --------------------------------
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12       # diffusive coefficient (um^2/s)
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)         # rotational diffusive coefficient
tr=(Dr)^(-1)                                                # rotational time scale
tauMax=min(length(dfMSDp[!,:xMSD]),length(dfMSDa[!,:xMSD])) # used to choose what kind of fit to perform 

cut=0   # cut the final points (when there are broken trajectories leading to stairs in the MSD plot - useful only in the linear regime)
ylimMSD=12.1
xlimMSD=tauMax/framerate-cut/framerate
lfit=50 # on how many poits the fitting is performed


##--- Initialize Plot with the 2 MSD ---------------------
graphMSD=plot()
plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], yerror=dfMSDp[!,:yerror], xlims=(-0.10,xlimMSD), ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label ="Brownian")
plot!(dfMSDa[!,:xMSD], dfMSDa[!,:MSD], yerror=dfMSDa[!,:yerror], ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label="Active")


##--- TWO MODELS: Ballistic & Diffusive regime -----------
y_par(t,D,V2)=4D*t+V2*(t^2)
y_lin(t,D,V,q)=4D*t+V^2*tr*t+q


##--- Inside the IF you define which model to use and so wheter to take the start or the end points ---

if  (0.1*tr*framerate)>lfit  # tr>>tauMax, BALLISTIC regime, parabolic fitting (0.1 is one o.d.g. less than tr, same for the  coefficient 10 for the elseif)
    x1=fill(1, lfit) # id number-like
    x2=fill(2, lfit)
    prex=[dfMSDp[1:lfit,:xMSD] ; dfMSDa[1:lfit,:xMSD]]
    preid=[x1; x2] # choosing between 1 or 2 depending on which fit is to be done
    x=[prex preid]
    y=[dfMSDp[1:lfit,:MSD]; dfMSDa[1:lfit,:MSD]]
    model=model_par
    yfun = y_par
    p0=[D,0.1] # first guess


elseif tauMax>(10*tr*framerate)+lfit  #tr<<tauMax, DIFFUSIVE regime, linear fitting
    x1=fill(1, lfit) # id number-like
    x2=fill(2, lfit)
    prex=[dfMSDp[end-lfit+1-cut:end-cut,:xMSD] ; dfMSDa[end-lfit+1-cut:end-cut,:xMSD]]
    preid=[x1; x2] # choosing between 1 or 2 depending on which fit is to be done
    x=[prex preid]
    y=[dfMSDp[end-lfit+1-cut:end-cut,:MSD]; dfMSDa[end-lfit+1-cut:end-cut,:MSD]]
    model=model_lin
    yfun = y_lin
    p0=[D,0.1, 0.1] # first guess

end

##--- Then you perform the fitting with those indications ---
fit=LsqFit.curve_fit(model,x,y,p0,lower=[0.9*D,0.0],upper=[1.1*D,Inf]) 
p=fit.param
XMSD=dfMSDa[1:tauMax,:xMSD]
yfit=yfun.(XMSD,p...)
plot!(XMSD,yfit, legend=:topleft, label="Fit")
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
velox=string(round(p[2],digits=2))
print("v diff = ")
println(velox)
title!("v diff= "*velox*" μm/s")

##--- Displaying and saving the result ---------------------
display(graphMSD)
png(graphMSD, path2*"MSDap_"*filename2)