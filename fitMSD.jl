using CSV, DataFrames, Plots, JSON3, LsqFit, Statistics, Dates #CurveFit
gr()    #backend dei plot, cerca figure interattive

##--- Brownian MSD-------------------------------------
folder="Sio2inMilliQ\\"
filename="MSD_SiO2millliQ_frame_m4305Sep22_1550"
##--- Active MSD---------------------------------------
folder="Sio2inMilliQ\\"
filename="MSD_SiO2millliQ_frame_m4305Sep22_1550"

## Read the data file and save it to a dataframe
path="Results\\"*folder*filename
dfMSD = CSV.read(path*".csv", DataFrame)

##--- Needed parameters---------------------------------
framerate=12
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)





## ---Fit----------------------------------------------
if  (0.5*tr*framerate)>10  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi, già coi .vmw che sono a 20fps prende lui.
    # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
    model(t,p)=4*p[1].*t.+(p[2]).*(t).^2
    p0=[D,0.1] #first guess
    fit2=LsqFit.curve_fit(model,dfMSD[1:5,:xMSD],dfMSD[1:5,:MSD],p0,lower=[0.2*D,0.0],upper=[5*D,10])
    p=fit2.param
    plot!(xMSD,model(xMSD,p), label="Fit")
    #confidence_inter = confidence_interval(fit2, 0.05)
    velox=string(round(sqrt(p[2]),digits=2))
    print("v = ")
    println(velox)
    title!("v= "*velox)

elseif tauMax>(2*tr*framerate)+10  #tr<<tauMax, DIFFUSIVE regime, linear fit, ma abbiamo video troppo corti, se riesci a farli di circa 1 min moltiplica per 5 invece
    model(t,p)=4*p[1]*D*t.-p[2]
#    model(t,p)=4*(p[1]*D+(1/4)*(p[2])*tr)*t.-(p[2])*(((tr)^2)/2)
    p0=[1.0,1.0] #first guess
#    fit2=LsqFit.curve_fit(model,xMSD[floor(Int,(2*tr*framerate)):end],MSD[floor(Int,(2*tr*framerate)):end],p0)
    fit2=LsqFit.curve_fit(model,xMSD[end-9:end],MSD[end-9:end],p0, lower=[0.2,0.0],upper=[5.0,Inf])
    p=fit2.param
    plot!(xMSD,model(xMSD,p), label="Fit")
#    confidence_inter = confidence_interval(fit2, 0.05)
#    velox=string(round(sqrt(p[2]),digits=2))
    velox=string(round(sqrt(4(p[1]-D)/tr),digits=2))
    print("v = ")
    println(velox)
    title!("v= "*velox)
end 
function model(t,p)
    y= 4*p[1]*D*t.-p[2]
    return y
end

display(graphMSD)


#png(graphMSD, path*"MSDfit_"*filename*DateTime)