using DataFrames, Plots, CurveFit
df = DataFrame(X=[25,25,25,50,50,50,60,60,60,90,90,90], Y=[16.9,18.4,15.6,19.6,21.6,33.4,33.2,37.1,25.5,48,56.6,67.5])
rate=scatter(df[!,:X],df[!,:Y], legend=false)
#df[!,:X] = categorical(df[!,:X],compress=true)
#gdf = groupby(df,:X,sort=true)
#mean

fit=linear_fit(df[!,:X],df[!,:Y])
yfit=fit[2].*df[!,:X]
plot!(df[!,:X],yfit)
xlabel!("t [s]");
ylabel!("thikness [nm]")

png(rate, "Rate_Pt")


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



yfit(t,D0,V,q)= 4D0.*t.+V.^2*tr.*t.+q
#plot!(dfMSDp[!,:xMSD],yfit(dfMSDp[!,:xMSD],p[1],p[2],p[3]), legend=true, label="Fit")
plot!(dfMSDa[!,:xMSD],yfit(dfMSDa[!,:xMSD],p[1],p[2],p[3]), legend=:topleft, label="Fit")
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")
velox=string(round(p[2],digits=2))
print("v diff = ")
println(velox)
title!("v diff= "*velox*" μm/s")

display(graphMSD)