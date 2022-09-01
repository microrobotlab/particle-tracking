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

## Define the zero position
x0 = df[1,:x]
y0 = df[1,:y]
## Add column and fill it with data
df[!,:SD] = ((df[!,:x]).-x0).^2+((df[!,:y]).-y0).^2
#df[!,:MSD] = 

plot(gdf[1][!,:Time],gdf[1][!,:SD],marker=true,label="SD",legend=:bottomright)

x0 = gdf[1][1,:x]
y0 = gdf[1][1,:y]
gdf[1][!,:SD] = ((gdf[1][!,:x]).-x0).^2+((gdf[1][!,:y]).-y0).^2

plot!(gdf[1][!,:Time],gdf[1][!,:SD],label="SD",legend=:bottomright)

x0 = gdf[22][1,:x]
y0 = gdf[22][1,:y]
gdf[22][!,:SD] = ((gdf[22][!,:x]).-x0).^2+((gdf[22][!,:y]).-y0).^2

plot!(gdf[22][!,:Time],gdf[22][!,:SD],label="SD",legend=:bottomright)