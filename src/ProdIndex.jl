###############################################################
# Part of BDisposal.jl package
###############################################################

"""
    prodIndex(gI,gO,bO,bI;retToScale,prodStructure,convexAssumption,startθ,startμ,startλ)

Compute productivity indexes


Given a set of measures of inputs, "good" ("desiderable") and "bad" ("undesiderable") outputs for different decision making units, compute their
productivity indexes improvements (or declines) between consecutive time periods.

## Parameters:
- Positional
  *  `gI`: "Good" inputs (3D matrix by DMUs, input items and periods)
  *  `gO`: "Good" outputs (3D matrix by DMUs, input items and periods)
  *  `bO`: "Bad" outputs (3D matrix by DMUs, input items and periods)
  *  `bI`: "Bad" inputs (optional 3D matrix by DMUs, input items and periods) [default: empty array]
- Keyword
  *  `retToScale`: Wheter the return to scales should be assumed "constant" (default) or "variable"
  *  `prodStructure`: Wheter the production structure should be assumed "additive" (default) or "multiplicative"
  *  `convexAssumption`: Wheter a convex production frontier should be assumed [default: `true`]

## Returns:
* A matrix of production indexes by DMUs and period passages (e.g. "year2 on year1" and "year3 on year2")

## Description of the function

`prodIndex()` displays environmental productivity indices. These
productivity measures are implemented for different time periods (eg., years, months
etc.) or spatial units (eg., countries, cities etc.), based on the  environmental efficiency
indicators described in `efficiencyScores`. Hence, the environmental productivity indices inherit the
structure of additive and multiplicative productivity measures.

## Interpretation of the results

The additive productivity indicator shows combined desirable and undesirable outputs
productivity improvement (respectively decline) when it takes positive (respectively negative) values.
In the multiplicative context, if the productivity measure is greater (respectively lesser) than
1 then, desirable and undesirable outputs productivity increase (respectively decrease)
arises. The BDisposal package underscores the prominent sources of environmental productivity
change. The main drivers of productivity variation are technological change,
technical efficiency variation and scale efficiency change. For the additive background,
when the technological change is greater (respectively, lesser) than 0 then, technolog-
ical improvement (respectively, deterioration) occurs. A similar reasoning applies for
the additive technical and scale efficiency components. In the multiplicative context, if
the technological change is greater (respectively, lesser) than 1 then, technological in-
crease (respectively, decrrease) arises. A similar reasonnng applies for the multiplicative
technical and scale efficiency components.

## Notes:
* The non-convex problem still need to consider bad inputs, consider directions, consider multiplicative/additive structures

"""
function prodIndex(gI::Array{Float64,3},gO::Array{Float64,3},bO::Array{Float64,3},bI::Array{Float64,3}=Array{Float64}(undef, (size(gI,1),0,size(gI,3)));
                   retToScale="constant",prodStructure="multiplicative",convexAssumption=true,)

    nDMUs = size(gI,1)
    nPer = size(gI,3)
    nBI  = size(bI,2)
    prodIndexes = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))

    noBadIndexDefault = prodStructure == "multiplicative" ? 1.0 : 0.0

    for t in 1:nPer-1
        gIₜ = gI[:,:,t]
        gIᵤ = gI[:,:,t+1]
        bIₜ = bI[:,:,t]
        bIᵤ = bI[:,:,t+1]
        gOₜ = gO[:,:,t]
        gOᵤ = gO[:,:,t+1]
        bOₜ = bO[:,:,t]
        bOᵤ = bO[:,:,t+1]
        dirGI = prodStructure == "multiplicative" ? (-1,0,0,0) : (1,0,0,0)
        dirBI = prodStructure == "multiplicative" ? (0,-1,0,0) : (0,1,0,0)
        dirGO = prodStructure == "multiplicative" ? (0,0,1,0)  : (0,0,1,0)
        dirBO = prodStructure == "multiplicative" ? (0,0,0,-1) : (0,0,0,-1)
        for z in 1:nDMUs
            gIₜ₀ = gIₜ[z,:]
            bIₜ₀ = bIₜ[z,:]
            gOₜ₀ = gOₜ[z,:]
            bOₜ₀ = bOₜ[z,:]
            gIᵤ₀ = gIᵤ[z,:]
            bIᵤ₀ = bIᵤ[z,:]
            gOᵤ₀ = gOᵤ[z,:]
            bOᵤ₀ = bOᵤ[z,:]

            forceLinearModel = (convexAssumption == true) ?  true : false

            # t...
            idx_gi_t̃ = problem(gIᵤ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI,crossTime=true)
            idx_gi_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI)
            idx_bi_t̃ = (nBI == 0) ? noBadIndexDefault : problem(gIₜ₀,bIᵤ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI,crossTime=true)
            idx_bi_t = (nBI == 0) ? noBadIndexDefault : problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI)
            idx_go_t̃ = problem(gIₜ₀,bIₜ₀,gOᵤ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO,crossTime=true)
            idx_go_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO)
            idx_bo_t̃ = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO,crossTime=true)
            idx_bo_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO)

            # u...
            idx_gi_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI)
            idx_gi_ũ = problem(gIₜ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI,crossTime=true)
            idx_bi_u = (nBI == 0) ? noBadIndexDefault : problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI)
            idx_bi_ũ = (nBI == 0) ? noBadIndexDefault : problem(gIᵤ₀,bIₜ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI,crossTime=true)
            idx_go_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO)
            idx_go_ũ = problem(gIᵤ₀,bIᵤ₀,gOₜ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO,crossTime=true)
            idx_bo_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO)
            idx_bo_ũ = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO,crossTime=true)

            if prodStructure == "multiplicative"
               idx_i_t = (idx_gi_t̃/idx_gi_t) * (idx_bi_t̃/idx_bi_t)
               idx_o_t = (idx_go_t̃/idx_go_t) * (idx_bo_t̃/idx_bo_t)
               idx_t   = idx_o_t/idx_i_t
               idx_i_u = (idx_gi_u/idx_gi_ũ) * (idx_bi_u/idx_bi_ũ)
               idx_o_u = (idx_go_u/idx_go_ũ) * (idx_bo_u/idx_bo_ũ)
               idx_u   = idx_o_u/idx_i_u
               idx     = (idx_t * idx_u)^(1/2)

            else
                idx_i_t = (idx_gi_t̃ - idx_gi_t) + (idx_bi_t̃ - idx_bi_t)
                idx_o_t = (idx_go_t - idx_go_t̃) + (idx_bo_t - idx_bo_t̃)
                idx_t   = idx_o_t - idx_i_t
                idx_i_u = (idx_gi_u - idx_gi_ũ) + (idx_bi_u - idx_bi_ũ)
                idx_o_u = (idx_go_ũ - idx_go_u) + (idx_bo_ũ - idx_bo_u)
                idx_u   = idx_o_u - idx_i_u
                idx     = (idx_t + idx_u) / 2
            end

            prodIndexes[z,t] = idx
        end
    end
    return prodIndexes
end
