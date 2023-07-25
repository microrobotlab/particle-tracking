# Importing packages 
using DataFrames, CSV

## Read the data file and save it to a dataframe
df = CSV.File("Results\\test\\data_example.csv") |> DataFrame

# copia del dataframe originale
df0 = df
df0_copy = copy(df)

println("\ndf")
display(df)

gdf = groupby(df,:title2)

for g in gdf
    g[:,:title4] = ones(nrow(g))
end

println("\ndf after changes in gdf")
display(df)

gdf2 = gdf
for g in gdf2
    g[:,:title4] = zeros(nrow(g))
end

println("\ndf after changes in second gdf")
display(df)

println("\ndf 'saved' at the beginning")
display(df0)
println("\ndf COPIED at the beginning")
display(df0_copy)

CSV.write("data_example_modified.csv",df)

# gdf3 = copy(gdf)