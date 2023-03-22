using CSV
using DataFrames

#saves data in a dataframe in .csv file. 4 columns: blob ID, time, x and y for each frame.
#framerate is the frame of the video
function save_data(result,vid_framerate,resultfilename)
    
    blobid = []
    time = []
    coord_x = []
    coord_y = []
    for i in eachindex(result.blobs)
        for j in eachindex(result.blobs[i].trace)
            coords = result.blobs[i].trace[j]
            push!(blobid,i)
            push!(time,j/vid_framerate)
            push!(coord_x,coords[2]) #invertiti??
            push!(coord_y,coords[1])

        end
    end

    #creates new file
    touch(resultfilename)

    efg = open(resultfilename, "w")
    
    #creating DataFrame
    data = DataFrame(BlobID = blobid,
    Time = time,
    x = coord_x,
    y= coord_y) 

    CSV.write(resultfilename, data)
    #return blobid, time, coord_x, coord_y
end

#brouillon 

# function save_data(result)
#     blobid = []
#     time = []
#     coord_x = []
#     coord_y = []
#     for i in eachindex(result.blobs)
#         for j in eachindex(result.blobs[i].trace)
#             coords = result.blobs[i].trace[j]
#             push!(blobid,i)
#             push!(time,j/12)
#             push!(coord_x,coords[1])
#             push!(coord_y,coords[2])

#         end
#     end
#     return blobid, time, coord_x, coord_y
# end

# blobid, time, coord_x, coord_y = save_data(result)


# #Creating DataFrame
# data = DataFrame(BlobID = blobid,
#                Time = time,
#                x = coord_x,
#                y= coord_y) 

# # modifying the content of myfile.csv using write method
# CSV.write("C:\\Users\\Yasmine\\SANDBOX\\git files\\pontedera\\results\\gaia\\coordinates\\coordinates.csv", data)

#xml