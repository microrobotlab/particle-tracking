using DataFrames
using NaNStatistics

function MSDfun(track,tauMax)   # tracks are yours gdf_clean_corrected, tauMax=ceil(Int,lt/10) with lt = max size gdf 
    ltrack=length(track[!,:Time])   # here uses :Time, but could have been anything else, it's just to get the length of each blob in the cycle
    tMax=ceil(Int, ltrack/10)    # divided by 10, value taken from the article Wang2021 (Practical guide MSD) --> SAME AS TAUMAX BUT on the single track => at the greatest it can be same, otherwise smaller
    msd=fill(NaN, tauMax+1)
    msd[1:tMax+1].=0
    for tau in 1:tMax 
        for i in tau+1:ltrack
            msd[tau+1]+=((track[i,:x]-track[i-tau,:x])^2+(track[i,:y]-track[i-tau,:y])^2)/(ltrack-tau)
        end
    end
#    println(length(msd))
    return msd

end
