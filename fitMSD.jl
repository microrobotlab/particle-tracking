using CSV, DataFrames, Plots, JSON3, LsqFit, Statistics, Dates #CurveFit
gr()    #backend dei plot, cerca figure interattive

##--- Brownian MSD-------------------------------------
filename1="MSD_movie051.avi"
##--- Active MSD---------------------------------------
filename2="MSD_movie054.avi" 

path="Results\\20221020\\J9\\"
folderDEST="vs_m51\\"

## Read the data file and save it to a dataframe
path1=path*filename1
dfMSDp = CSV.read(path1*".csv", DataFrame)
path2=path*filename2
dfMSDa = CSV.read(path2*".csv", DataFrame)
diamPart=1

##--- Needed parameters---------------------------------
framerate=12
D=(1.380649e-23*298)/(6π*1e-3*(diamPart*1e-6/2))*1e12     #um^2/s
Dr=(1.380649e-23*298)/(8π*1e-3*(diamPart*1e-6/2)^3)
tr=(Dr)^(-1)

#cut the final points (whe there are broken trajectories leading to stairs in the MSD plot)
cut=0

#ylims=6.1
ylimMSD=20.1
xlimMSD=6-cut/12
lfit=10

##---Initialize Plot------------------------------------
graphMSD=plot()
plot!(dfMSDp[!,:xMSD], dfMSDp[!,:MSD], yerror=dfMSDp[!,:yerror], xlims=(-0.10,xlimMSD), ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label ="Brownian")
plot!(dfMSDa[!,:xMSD], dfMSDa[!,:MSD], yerror=dfMSDa[!,:yerror], ylims=(-0.10,ylimMSD), marker=true, legend=:topleft, label="Active")
#tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi. nel caso implementare con  # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
#tr<<tauMax, DIFFUSIVE regime, linear fit, se riesci a farli di circa 1 min poi prendi 1/10 del video invece di 1/5 taumax
x1=fill(1, lfit)
x2=fill(2, lfit)
prex=[dfMSDp[end-9-cut:end-cut,:xMSD] ; dfMSDa[end-9-cut:end-cut,:xMSD]]
preid=[x1; x2]
x=[prex preid]
#x=[dfMSDp[!,:xMSD] x1 ; dfMSDa[!,:xMSD] x2]
y=[dfMSDp[end-9-cut:end-cut,:MSD]; dfMSDa[end-9-cut:end-cut,:MSD]]
#p=zeros(Float64, 3, 1) 

function model(x,p)
    t=x[:,1]
    id=x[:,2]
    D0=p[1]
    V=zeros(length(id))
    V[id.==2.0].=p[2]
    q=zeros(length(id))
    q[id.==2.0].=p[3]
    yf= 4D0.*t.+V.^2*tr.*t.+q
    return yf
end


p0=[D,0.1, 0.1] #first guess
fit=LsqFit.curve_fit(model,x,y,p0) #,p0,lower=[0.2*D,0.0],upper=[5*D,10])
p=fit.param
yfit(t,D0,V,q)= 4D0.*t.+V.^2*tr.*t.+q
#plot!(dfMSDp[!,:xMSD],yfit(dfMSDp[!,:xMSD],p[1],p[2],p[3]), legend=true, label="Fit")
plot!(dfMSDa[!,:xMSD],yfit(dfMSDa[!,:xMSD],p[1],p[2],p[3]), legend=:topleft, label="Fit")
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
velox=string(round(p[2],digits=2))
print("v diff = ")
println(velox)
title!("v diff= "*velox*" μm/s")

pathDEST=path*folderDEST
display(graphMSD)
#Date_Time= Dates.format(now(), "dduyy_HHMM") 
png(graphMSD, pathDEST*"MSDap_"*filename2)




## ---OLD Fit------------------------------------------

#if  (0.5*tr*framerate)>5  #tr>>tauMax, BALLISTIC regime, fit parabolico inizio, ma abbiamo frame rate troppo bassi
    # model(t,MSD,D)=4D.*t.+(V.^2).*(t)^2
    #model(t,p)=4*p[1].*t.+(p[2]).*(t).^2
    #p0=[D,0.1] #first guess
    #fit2=LsqFit.curve_fit(model,xMSD[1:5],MSD[1:5],p0,lower=[0.2*D,0.0],upper=[5*D,10])
    #p=fit2.param
    #plot!(xMSD,model(xMSD,p), label="Fit")
    #confidence_inter = confidence_interval(fit2, 0.05)
    #velox=string(round(sqrt(p[2]),digits=2))
    #print("v bal = ")
    #println(velox)
    #title!("v= "*velox)

#elseif tauMax>(2*tr*framerate)+5  #tr<<tauMax, DIFFUSIVE regime, linear fit, ma abbiamo video troppo corti, se riesci a farli di circa 1 min moltiplica per 5 invece
    #pfit=linear_fit(xMSD[end-4:end],MSD[end-4:end])
    #yfit=pfit[2].*xMSD.+pfit[1]
    #plot!(xMSD,yfit)
    #Deff=pfit[2]/4
    #velox=sqrt((pfit[2]-4D)*tr)
    #velox=sqrt(-2pfit[1]/tr^2)

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
#end 
