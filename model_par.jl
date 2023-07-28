using DataFrames

function model_par(x,p)
    t=x[:,1]
    id=x[:,2]
    D=p[1]
    V2=zeros(length(id))
    V2[id.==2.0].=p[2]
    yf= y_par.(t,D,V2)
    return yf
end