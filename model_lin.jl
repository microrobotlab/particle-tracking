using DataFrames

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
