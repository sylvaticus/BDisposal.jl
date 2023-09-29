###############################################################
# Part of BDisposal.jl package
###############################################################

"""
    prodIndexFB(gI,gO,bO,bI;remarcable_obs_dmu,remarcable_obs_period)

Compute productivity indexes with a fixed base.


Given a set of measures of inputs, "good" ("desiderable") and "bad" ("undesiderable") outputs for different decision making units, and a specific dmu/time observation to consider as base ("remarcable"), provides the productivity indexes at various time periods.

## Parameters:
- Positional
  *  `gI`: "Good" inputs (3D array by DMUs, input items and periods)
  *  `gO`: "Good" outputs (3D array by DMUs, input items and periods)
  *  `bO`: "Bad" outputs (3D array by DMUs, input items and periods)
  *  `bI`: "Bad" inputs (optional 3D array by DMUs, input items and periods) [default: empty array]
- Keyword
  *  `remarcable_obs_dmu`: Which observation to consider the "remarcable" one [def: `1`]
  *  `remarcable_obs_period`: Which tipe period to consider the "remarcable" one [def: `1``]

## Returns:
* A named touple where each element is a matrix of base-related production indexes by DMUs and time periods.
* Currently the value reported are:
    - `prodIndexes`:          Overall production indexes

## Description of the function

`prodIndexFB` is a fixed base version of the multiplicative `prodIndex` function. As a result, `prodIndexFB` is based upon the identification of a fixed base observation allowing to compare spatial and/or temporal observations. Specifically, the choice of the base observation for comparisons is affected by the identification of remarkable spatial and/or temporal units.
    
## Interpretation of the results
    
`prodIndexFB` > 1 indicates environmentally-adjusted performance improvement. In such case, the remarkable observation produces more good outputs and less bad goods than the compared observation for a given level of good and bad inputs. Alongside, less good and bad inputs are employed by the remarkable observation relatively to the compared observation for given good and bad outputs. If `prodIndexFB < 1`, the reverse reasonning holds.

## Notes:
* The function assumes convex, multiplicative production functions with variable returns to scale

"""
function prodIndexFB(gI::Array{Float64,3},gO::Array{Float64,3},bO::Array{Float64,3},bI::Array{Float64,3}=Array{Float64}(missing, (size(gI,1),0,size(gI,3)));remarcable_obs_dmu=1,remarcable_obs_period=1)

    # Working on the logs of the data
    gI = log.(gI); gO = log.(gO); bO = log.(bO); bI = log.(bI)

    nDMUs = size(gI,1)
    nPer  = size(gI,3)
    ngI, nbI, ngO, nbO = size(gI,2), size(bI,2), size(gO,2), size(bO,2)

    # Output containers

    # - Global
    prodIndexes = Array{Union{Float64,Missing}}(missing, (nDMUs,nPer))

    noBadIndexDefault = 1.0

    #=
    # Observation averaged across time..
    gIₐ = reshape(mean(gI,dims=(1,3)), ngI) 
    bIₐ = reshape(mean(bI,dims=(1,3)), nbI)
    gOₐ = reshape(mean(gO,dims=(1,3)), ngO)
    bOₐ = reshape(mean(bO,dims=(1,3)), nbO)
    =#

    # Min of the various observations across time and obs....
    #gIₘᵢₙ = reshape(minimum(gI,dims=(1,3)), ngI)
    bIₘᵢₙ = reshape(minimum(bI,dims=(1,3)), nbI)
    gOₘᵢₙ = reshape(minimum(gO,dims=(1,3)), ngO)
    #bOₘᵢₙ = reshape(minimum(bO,dims=(1,3)), nbO)

    # Mam of the various observations across time and obs....
    gIₘₐₓ = reshape(maximum(gI,dims=(1,3)), ngI)
    #bIₘₐₓ = reshape(maximum(bI,dims=(1,3)), nbI)
    #gOₘₐₓ = reshape(maximum(gO,dims=(1,3)), ngO)
    bOₘₐₓ = reshape(maximum(bO,dims=(1,3)), nbO)

    # Remarcable obs
    gIᵣ = gI[remarcable_obs_dmu,:,remarcable_obs_period]
    bIᵣ = bI[remarcable_obs_dmu,:,remarcable_obs_period]
    gOᵣ = gO[remarcable_obs_dmu,:,remarcable_obs_period]
    bOᵣ = bO[remarcable_obs_dmu,:,remarcable_obs_period]

    # The observations with all the time observations merged together for a given DMU/observation item
    gIₘ = vcat([s for s in eachslice(gI,dims=3)]...)
    bIₘ = vcat([s for s in eachslice(bI,dims=3)]...)
    gOₘ = vcat([s for s in eachslice(gO,dims=3)]...)
    bOₘ = vcat([s for s in eachslice(bO,dims=3)]...)

    dirGI = (1,0,0,0) 
    dirBI = (0,-1,0,0)
    dirGO = (0,0,1,0) 
    dirBO = (0,0,0,1) 

    idx_gi_r = problemFB(gIᵣ,bIₘᵢₙ,gOₘᵢₙ,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirGI) |> exp
    idx_bi_r = problemFB(gIₘₐₓ,bIᵣ,gOₘᵢₙ,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirBI) |> exp
    idx_go_r = problemFB(gIₘₐₓ,bIₘᵢₙ,gOᵣ,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirGO) |> exp
    idx_bo_r = problemFB(gIₘₐₓ,bIₘᵢₙ,gOₘᵢₙ,bOᵣ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirBO) |> exp

    for t in 1:nPer

        gIₜ = gI[:,:,t]
        bIₜ = bI[:,:,t]
        gOₜ = gO[:,:,t]
        bOₜ = bO[:,:,t]

        for z in 1:nDMUs
            gIₜ₀ = gIₜ[z,:]
            bIₜ₀ = bIₜ[z,:]
            gOₜ₀ = gOₜ[z,:]
            bOₜ₀ = bOₜ[z,:]

            idx_gi_t = problemFB(gIₜ₀,bIₘᵢₙ,gOₘᵢₙ,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirGI) |> exp
            idx_bi_t = (nbI == 0) ? noBadIndexDefault : problemFB(gIₘₐₓ,bIₜ₀,gOₘᵢₙ,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ, directions=dirBI) |> exp
            idx_go_t = problemFB(gIₘₐₓ,bIₘᵢₙ,gOₜ₀,bOₘₐₓ,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirGO) |> exp
            idx_bo_t = problemFB(gIₘₐₓ,bIₘᵢₙ,gOₘᵢₙ,bOₜ₀,gIₘ,bIₘ,gOₘ,bOₘ,directions=dirBO) |> exp

            # Computing aggregated for the full (gI,Bi,bO,gO) partition...
            #if prodStructure == "multiplicative"
                    idx_i_t = (idx_gi_t/idx_gi_r) * (idx_bi_t/idx_bi_r)
                    idx_o_t = (idx_go_r/idx_go_t) * (idx_bo_t/idx_bo_r)
                    idx_t   = idx_o_t/idx_i_t
                    idx = idx_t
            #else # additive
            #        idx_i_t = (idx_gi_r - idx_gi_t) + (idx_bi_r - idx_bi_t)
            #        idx_o_t = (idx_go_t - idx_go_r) - (idx_bo_r - idx_bo_t)
            #        idx_t   = idx_o_t - idx_i_t
            #        idx     = idx_t
            #end
            prodIndexes[z,t] = idx
        end # end for each DMU
    end # end for each period
    return (prodIndexes        =  prodIndexes,)
end
