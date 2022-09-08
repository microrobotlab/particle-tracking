using CSV, DataFrames, Plots, JSON3, CurveFit, Statistics, Dates #CurveFit
gr()    #backend dei plot, cerca figure interattive

##--- Brownian MSD-------------------------------------
folder="J1\\"
filename="MSD_J1brown_27_3march05Sep22_1019"
##--- Active MSD---------------------------------------
#folder="J1\\"
filename2="MSD_J1+H2O2_m40_2march05Sep22_1153"

## Read the data file and save it to a dataframe
path="Results\\"*folder*filename
dfMSDp = CSV.read(path*".csv", DataFrame)
path="Results\\"*folder*filename2
dfMSDa = CSV.read(path*".csv", DataFrame)

##--- Needed parameters---------------------------------
framerate=12
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)

#ylims=6.1
ylimMSD=6.1
xlimMSD=2

plot()
## ---Fit Brown----------------------------------------------
plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], yerror=dfMSDp[!,:yerror], ylims=(-0.10,ylimMSD), marker=true, legend=false)
#if  (0.5*tr*framerate)>10  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi, già coi .vmw che sono a 20fps prende lui.
    # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
    pfit=poly_fit(dfMSDp[1:5,:xMSD],dfMSDp[1:5,:MSD],2) #dovresti dirgli che pfit[1]=0
    yfit=pfit[1].+pfit[2].*dfMSDp[!,:xMSD].+pfit[3].*dfMSDp[!,:xMSD].^2
    plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], label="Fit")
    #confidence_inter = confidence_interval(fit2, 0.05)
    velox=string(round(sqrt(p[2]),digits=2))
    print("v = ")
    println(velox)
    title!("v bal = "*velox)

#elseif tauMax>(2*tr*framerate)+10  #tr<<tauMax, DIFFUSIVE regime, linear fit, ma abbiamo video troppo corti, se riesci a farli di circa 1 min moltiplica per 5 invece
    pfit2=linear_fit(dfMSDp[!,:xMSD],MSD[end-9:end])
    yfit2=pfit[2].*dfMSDp[!,:xMSD].+pfit2[1]
    plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], label="Fit2")
    Deff=pfit2[2]/4
    velox=sqrt((pfit2[2]-4D)*tr)
    velox=sqrt(-2pfit2[1]/tr^2)
#    confidence_inter = confidence_interval(fit2, 0.05)
#    velox=string(round(sqrt(p[2]),digits=2))
    velox2=string(round(sqrt(4(p[1]-D)/tr),digits=2))
    print("v = ")
    println(velox2)
    title!("v diff = "*velox)
#end 

## ---Fit Active----------------------------------------------
plot!(dfMSDa[!,:xMSD], dfMSDa[!,:MSD], yerror=dfMSDa[!,:yerror],xlims=(-0.10,xlimMSD) ylims=(-0.10,ylimMSD), marker=true, legend=false)
#if  (0.5*tr*framerate)>10  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi, già coi .vmw che sono a 20fps prende lui.
    # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
    pfit3=poly_fit(dfMSDp[1:10,:xMSD],dfMSDp[1:10,:MSD],2) #dovresti dirgli che pfit[1]=0
    yfit3=pfit3[1].+pfit3[2].*xMSD.+pfit3[3].*xMSD.^2
    plot(dfMSD[:,:xMSD],dfMSD[:,:MSD], yerror=dsMSD, ylims=(-0.10,YlimMSD), marker=true,legend=false)
    plot!(xMSD,model(xMSD,p), label="Fit")
    #confidence_inter = confidence_interval(fit2, 0.05)
    velox3=string(round(sqrt(p[2]),digits=2))
    print("v = ")
    println(velox3)
    title!("v bal = "*velox)

#elseif tauMax>(2*tr*framerate)+10  #tr<<tauMax, DIFFUSIVE regime, linear fit, ma abbiamo video troppo corti, se riesci a farli di circa 1 min moltiplica per 5 invece
    pfit4=linear_fit(xMSD[end-9:end],MSD[end-9:end])
    yfit4=pfit[2].*xMSD.+pfit[1]
    plot!(xMSD,yfit4)
    Deff=pfit4[2]/4
    velox3=sqrt((pfit4[2]-4D)*tr)
    velox3=sqrt(-2pfit4[1]/tr^2)
#    confidence_inter = confidence_interval(fit2, 0.05)
#    velox=string(round(sqrt(p[2]),digits=2))
    velox4=string(round(sqrt(4(p[1]-D)/tr),digits=2))
    print("v = ")
    println(velox2)
    title!("v diff = "*velox)
#end 

display(graphMSD)


#png(graphMSD, path*"MSDfit_"*filename*DateTime)