using DataFrames
using NaNStatistics
using Plots
using Statistics
using CSV

function drift_corr(df, um_px, framerate, filename)
    
    # derived information:
    px_um = 1/um_px     # pixel / micron
    max_act_vel = 100   # micron / s
    max_disp_frame = max_act_vel / framerate
    min_duration_s = 1.0 # minimum track length, s
    min_duration_frame = floor(Int,min_duration_s*framerate)
    Δt_drift = 10 #s

    # group and plot tracks
    gdf = groupby(df,:BlobID)

    plt_tracks=plot()
    for g in gdf
        plot!(g[:,:x],g[:,:y],legend=:none,color=:black)
    end
    xlabel!("x (μm)"); ylabel!("y (μm)")

    # find tracks of immobile (stuck) particles and camera artefacts
    immobile_tracks = Int[]
    for g in gdf
        rms_x_px = sqrt(sum((g[:,:x].-g[1,:x]).^2)/nrow(g)) * px_um
        rms_y_px = sqrt(sum((g[:,:y].-g[1,:y]).^2)/nrow(g)) * px_um
        if rms_x_px < 1.0 && rms_y_px < 1.0
            push!(immobile_tracks,g[1,:BlobID])
        end
    end

    # find tracks with unphysical jumps
    jump_tracks = Int[]
    for g in gdf
        dr = sqrt.(diff(g[:,:x]).^2 + diff(g[:,:y]).^2)
        max_dr = maximum(dr; init=0)
        if max_dr > max_disp_frame
            push!(jump_tracks,g[1,:BlobID])
        end
    end

    # find too-short tracks
    short_tracks = Int[]
    for g in gdf
        if nrow(g) < min_duration_frame
            push!(short_tracks,g[1,:BlobID])
        end
    end

    # tracks to be discarded
    discard_tracks = unique([immobile_tracks; jump_tracks;short_tracks])

    # clean grouped DataFrame without discarded tracks and plot tracks
    gdf_clean = gdf[[g[1,:BlobID] ∉ discard_tracks for g in gdf]]

    print("n. of tracks = ")
    println(length(gdf))
    print("n. of filtered tracks = ")
    println(length(gdf_clean))

    for g in gdf_clean
        plot!(plt_tracks,g[:,:x],g[:,:y],legend=:none,color=:green)
    end
    display(plt_tracks)

    # calculate drift
    # ===

    n_tracks = length(gdf)
    n_frames = maximum(nrow(g) for g in gdf) 

    # put good tracks in matrix
    x_m, y_m = fill(NaN,n_frames,n_tracks), fill(NaN,n_frames,n_tracks)
    p = [x_m,y_m]
    for g in gdf
        track = g[1,:BlobID]
        if track ∉ discard_tracks
            x_m[(g[:,:Frame]),track]=g[:,:x] #CAMBIATO FRAME A TIME
            y_m[(g[:,:Frame]),track]=g[:,:y] #CAMBIATO FRAME A TIME
        end
    end

    # calculate global displacement between frames
    dp = [nanmean(diff(pi,dims=1),dims=2) for pi in p]
    # plot(dp[1]); plot!(dp[2])

    # calculate drift
    drift=[[0;cumsum(dpi,dims=1)] for dpi in dp]
    plot(drift[1],label="drift x"); plot!(drift[2],label="drift y")

    # smoothing of the drift: is that really needed
    function mymovmean(s,w)
        ssm = fill(NaN,size(s))
        for i in eachindex(s)
            wi = minimum([w,2i-1,2(length(s)-i)+1])
            ssm[i] = mean(s[i-(wi÷2):i+(wi÷2)])
        end
        return ssm
    end

    w = Δt_drift*framerate +1 -mod(Δt_drift*framerate,2) #how much the drift is changing? Calculate over Δt to smooth
    # w = n_frames÷4+1-mod(n_frames÷4,2) #n.frames/4, but it needs to be odd
    # w = 1
    drift_sm = mymovmean.(drift,w)
    plot!(drift_sm[1],label="drift x sm"); plot!(drift_sm[2],label="drift y sm")

    # correct drift
    gdf_clean_corrected = gdf_clean
    for g in gdf_clean_corrected
        g[!,:x] -= drift_sm[1][g[:,:Frame]]
        g[!,:y] -= drift_sm[2][g[:,:Frame]]
    end

    for g in gdf_clean_corrected
        plot!(plt_tracks,g[:,:x],g[:,:y],legend=:none,color=:red)
    end
    display(plt_tracks)
    
    #---SAVE dataframe with drift correction & Plot-----------
    CSV.write(path*"MSDdriftCorr_"*filename*".csv", df)
    png(plt_tracks, path*"tracks_dc_"*filename)

    return gdf_clean_corrected, immobile_tracks, jump_tracks, short_tracks, discard_tracks
end