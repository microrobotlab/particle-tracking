using CSV, Pkg, DataFrames, CategoricalArrays, Plots
gr()    #backend dei plot, cerca figure interattive

## Read the data file and save it to a dataframe
df = CSV.read("h2o2_1_coordinates.csv", DataFrame)
# Make x and y columns float so as to allow scientific notation in df
df[!,:x] = df[!,:x]*1.0
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
graphSD=plot(gdf[1][!,:Time],gdf[1][!,:SD],label="SD "*"1",legend=:bottomright)

#gdf qui dice sono 243 element, ma poi mi stampa 222 etichette nel ciclo..

for i in 2:length(gdf)

#    if (length(gdf[i][!,:x])>50) && (sum(gdf[i][!,:SD])) #sommare tutti gli el di SD per eliminare le traiettorie stuopide dove non si muove nulla
    if length(gdf[i][!,:x])>50
        for n in 1:length(gdf[i][!,:x])
            if skipmissings((gdf[i][n,:SD]+gdf[i][n+1,:SD])<1000)
          
                ## Define the zero position
                x0i= gdf[i][1,:x]
                y0i = gdf[i][1,:y]
                ## Add column and fill it with data
                gdf[i][!,:SD] = ((gdf[i][!,:x]).-x0i).^2. +((gdf[i][!,:y]).-y0i).^2
                plot!(gdf[i][!,:Time],gdf[i][!,:SD],label = string("SD ",i),legend=:bottomright)
                #graphSD = plot!(gdf[i][!,:Time],gdf[i][!,:SD],label="SD",legend=:bottomright)
            end
        end    
    end
end

# MSD =
#filtrare traiettorie troppo corte
display(graphSD)
