using CSV, Pkg, DataFrames, CategoricalArrays, Plots, NaNStatistics 
gr()    #backend dei plot, cerca figure interattive

## Read the data file and save it to a dataframe
df = CSV.read("h2o2_1_coordinates.csv", DataFrame)
# Make x and y columns float and convert pixels to microns
df[!,:x] = df[!,:x]*1.0 #mettici fattore di conversione
df[!,:y] = df[!,:y]*1.0

## Make one column categorical
df[!,:BlobID] = categorical(df[!,:BlobID],compress=true)
## Group dataframe by values in categorical column
gdf = groupby(df,:BlobID,sort=true)

## Define the zero position for the fiirst particle
x0 = gdf[1][1,:x]
y0 = gdf[1][1,:y]
## Add column and fill it with data
gdf[1][!,:SD] = ((gdf[1][!,:x]).-x0).^2+((gdf[1][!,:y]).-y0).^2
nTraks=length(gdf)

#------> Azzera Plot
#graphSDtrck =plot(gdf[1][!,:Time],gdf[1][!,:SD],label="SD "*"1",legend=:bottomright);
graphSDtrck =plot()

for i in 1:length(gdf)
   
    ## Define the zero position
    x0i= gdf[i][1,:x]
    y0i = gdf[i][1,:y]
    ## Add column and fill it with data
    gdf[i][!,:SD] = ((gdf[i][!,:x]).-x0i).^2. +((gdf[i][!,:y]).-y0i).^2
    plot!(graphSDtrck, gdf[i][!,:Time],gdf[i][!,:SD],label = string("SD ",i),legend=:bottomright)         
end

display(graphSDtrck)

#filtra le traiettorie troppo corte o con salti troppo elevati tra 2 punti

#------> Azzera Plot
graphSD=plot();

# crea un altro oggetto con solo le traiettorie filtrate per mediare solo quelle x MSD
#dff = DataFrame(BlobID = Any[], Time = Any[], x = Any[], y = Any[], SD = Any[])

# Trova Max length vettore tempo
lt=0
for g in gdf
    if size(g)[1]>lt 
        lt= size(g)[1]
    end
end

idx=[]
matr=fill(NaN, lt, nTraks)

for i in 1:length(gdf)
    #    if (length(gdf[i][!,:x])>50) && (sum(gdf[i][!,:SD])) #sommare tutti gli el di SD per eliminare le traiettorie stuopide dove non si muove nulla
    flag=false
    ltrack=length(gdf[i][!,:x])
    if (ltrack>12) && (maximum(gdf[i][!,:SD])>10)
        for n in 1:(ltrack-1)
            if abs((gdf[i][n,:SD]-gdf[i][n+1,:SD]))>800
                flag=true
                break
            end
        end  
        if flag==false
            # crea un altro oggetto con solo le traiettorie filtrate per mediare solo quelle x MSD
 #           push!(dff, [[gdf[i][!, :BlobID]] [gdf[i][!, :Time]] [gdf[i][!, :x]] [gdf[i][!, :y]] [gdf[i][!, :SD]]])
            #push!(dff, [gdf[i]])
            push!(idx, i)
            matr[1:ltrack, i]=gdf[i][!,:SD]
            #graphSD = plot!(gdf[i][!,:Time],gdf[i][!,:SD],label = string("SD ",i),legend=:bottomright)
            plot!(gdf[i][!,:Time],gdf[i][!,:SD],legend=false)
        end
    end
end

display(graphSD)

MSD=nanmean(matr, dims=2)
t=0:(1/12):(lt-1)/12
plot!(t,MSD, marker=true)

# MSD =  PROBLEMA! Se non cancello dal df le cose che ho filtrato, qua le medio!
#           meglio rifare il ciclo o nel ciclo precedente dirgli di cancellare 
#           o di creare un'altra struttura con quelle che salvo?

## Make one column categorical
dff[!,:Time] = categorical(df[!,:Time],compress=true)
## Group dataframe by values in categorical column
tdf = groupby(dff,:Time,sort=true)
# Dimmi quale è il tdf più lungo, e prendo lui come vettore t tempo (ora è 1)
#...

t=gdf[1][!,:Time]
#l=length(tdf[1][!,:Time])
l=length(tdf)
MSD= Array{Float64}(undef, l, 1)

for i in 1:length(tdf)
   
    MSD[i]= sum(tdf[i][!,:SD])/length(tdf[i][!,:SD])
      
end

graphSD = plot!(gdf[1][!,:Time],MSD,marker=true,label = "MSD",legend=false)   


# OOOOOOor:


#df[!,:MSD]=[!,!]
#for i in 1:length(tdf)
   
#    ## Add column and fill it with data
#    tdf[i][1,:MSD]= sum(tdf[i][!,:SD])/length(tdf[i][!,:SD])
      
#end

#graphSD = plot!(gdf[1][!,:Time], df[!,:MSD],label = "MSD",legend=:bottomright)       
