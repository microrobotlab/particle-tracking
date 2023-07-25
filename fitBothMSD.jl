using CSV, DataFrames, Plots, JSON3, LsqFit, Statistics, Dates #, CurveFit
gr()    #backend dei plot, cerca figure interattive

diamPart=3
framerate=25

##--- Brownian MSD-------------------------------------
folder1="Results\\3_um\\Pt\\20230307_P01_E3007 - J26\\1x\\"
filename1="MSD_VID02924.AVI" #eventualmente 70 e 71
##--- Active MSD---------------------------------------
folder2="Results\\3_um\\Pt\\20230307_P01_E3007 - J26\\1x\\"
filename2="MSD_VID02922.AVI" #
path=""#*folder

## Read the data file and save it to a dataframe
path1=path*folder1*filename1
dfMSDp = CSV.read(path1*".csv", DataFrame)
path2=path*folder2*filename2
dfMSDa = CSV.read(path2*".csv", DataFrame)

##--- Needed parameters---------------------------------
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3) #only if this is highet than the max sampled rate it make sense to perform the linear fit..
tr=(Dr)^(-1)
tauMax=min(length(dfMSDp[!,:xMSD]),length(dfMSDa[!,:xMSD]))

#cut the final points (whe there are broken trajectories leading to stairs in the MSD plot)
cut=0 # useful only in the linear regime

ylimMSD=12.1
xlimMSD=tauMax/framerate-cut/framerate
lfit=50 # on how many poits it performs the ?lin? fitting


##---Initialize Plot with the 2 MSD-------------------------
graphMSD=plot()
plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], yerror=dfMSDp[!,:yerror], xlims=(-0.10,xlimMSD), ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label ="Brownian")
plot!(dfMSDa[!,:xMSD], dfMSDa[!,:MSD], yerror=dfMSDa[!,:yerror], ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label="Active")


#---TWO MODELS: Ballistic & Diffusive regime-------------
y_par(t,D,V2)=4D*t+V2*(t^2)
y_lin(t,D,V,q)=4D*t+V^2*tr*t+q

function model_par(x,p)
    t=x[:,1]
    id=x[:,2]
    D=p[1]
    V2=zeros(length(id))
    V2[id.==2.0].=p[2]
    yf= y_par.(t,D,V2)
    return yf
end

function model_lin(x,p)
    t=x[:,1]
    id=x[:,2]
    D0=p[1]
    V=zeros(length(id))
    V[id.==2.0].=p[2]
    q=zeros(length(id))
    q[id.==2.0].=p[3]
    yf= y_lin.(t,D0,V,q)
    return yf
end


#---Inside the IF you define which model to use and so wheter to thake the start or the end points--------------

if  (0.1*tr*framerate)>lfit  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio # metto 0.1, cioè un odg meno di tr, e uguale sotto x10
    x1=fill(1, lfit) #id number-like
    x2=fill(2, lfit)
    prex=[dfMSDp[1:lfit,:xMSD] ; dfMSDa[1:lfit,:xMSD]]
    preid=[x1; x2] # choosing between 1 or 2 depending on which fit is to be done
    x=[prex preid]
    #x=[dfMSDp[!,:xMSD] x1 ; dfMSDa[!,:xMSD] x2]
    y=[dfMSDp[1:lfit,:MSD]; dfMSDa[1:lfit,:MSD]]
    #p=zeros(Float64, 3, 1) 
    model=model_par
    yfun = y_par
    p0=[D,0.1] #first guess



elseif tauMax>(10*tr*framerate)+lfit  #tr<<tauMax, DIFFUSIVE regime, linear fit #tauMax>(10*tr*framerate)+10

    x1=fill(1, lfit) #id number-like
    x2=fill(2, lfit)
    prex=[dfMSDp[end-lfit+1-cut:end-cut,:xMSD] ; dfMSDa[end-lfit+1-cut:end-cut,:xMSD]]
    preid=[x1; x2] # choosing between 1 or 2 depending on which fit is to be done
    x=[prex preid]
    #x=[dfMSDp[!,:xMSD] x1 ; dfMSDa[!,:xMSD] x2]
    y=[dfMSDp[end-lfit+1-cut:end-cut,:MSD]; dfMSDa[end-lfit+1-cut:end-cut,:MSD]]
    #p=zeros(Float64, 3, 1) 
    yfun = y_lin
    p0=[D,0.1, 0.1] #first guess
    model=model_lin


end

#---Then you perform the fitting with those indications---
fit=LsqFit.curve_fit(model,x,y,p0,lower=[0.9*D,0.0],upper=[1.1*D,Inf]) #,p0,lower=[0.2*D,0.0],upper=[5*D,10])
p=fit.param
#yfit(t,D0,V,q)= 4D0.*t.+V.^2*tr.*t.+q
XMSD=dfMSDa[1:tauMax,:xMSD]
#y0b=fit.(XMSD)
yfit=yfun.(XMSD,p...)#dfMSDa[!,:xMSD],p)
#plot!(dfMSDp[!,:xMSD],yfit(dfMSDp[!,:xMSD],p[1],p[2],p[3]), legend=:topleft, label="Fit")
plot!(XMSD,yfit, legend=:topleft, label="Fit")
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
velox=string(round(p[2],digits=2))
print("v diff = ")
println(velox)
title!("v diff= "*velox*" μm/s")

display(graphMSD)
#Date_Time= Dates.format(now(), "dduyy_HHMM") 
png(graphMSD, path2*"MSDap_"*filename2)




## ---OLD Fit------------------------------------------

# if  (0.1*tr*framerate)>10  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio # metto 0.1, cioè un odg meno di tr, e uguale sotto x10
#     model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
#     model(t,p)=4*p[1].*t.+(p[2]).*(t).^2
#     p0=[D,0.1] #first guess
#     fit2=LsqFit.curve_fit(model,xMSD[1:5],MSD[1:5],p0,lower=[0.2*D,0.0],upper=[5*D,10])
#     p=fit2.param
#     plot!(xMSD,model(xMSD,p), label="Fit")
#     confidence_inter = confidence_interval(fit2, 0.05)
#     velox=string(round(sqrt(p[2]),digits=2))
#     print("v bal = ")
#     println(velox)
#     title!("v= "*velox)

# elseif tauMax>(10*tr*framerate)+10  #tr<<tauMax, DIFFUSIVE regime, linear fit
#     pfit=linear_fit(xMSD[end-4:end],MSD[end-4:end])
#     yfit=pfit[2].*xMSD.+pfit[1]
#     plot!(xMSD,yfit)
#     Deff=pfit[2]/4
#     velox=sqrt((pfit[2]-4D)*tr)
#     velox=sqrt(-2pfit[1]/tr^2)

#        fit=curve_fit(LinearFit, xMSD[end-4:end],MSD[end-4:end])
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
# end 

