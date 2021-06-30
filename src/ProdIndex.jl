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
* A named touple where each element is a matrix of production indexes by DMUs and period passages (e.g. "year2 on year1" and "year3 on year2") or one of their decompositions.
* Currently the value reported are:
    - `prodIndexes`:          Overall production indexes
    - `prodIndexes_G`:        Decomposition for "good" inputs and outputs
    - `prodIndexes_B`:        Decomposition for "bad" inputs and outputs
    - `prodIndexes_T`:        Decomposition for the technological component, overall
    - `prodIndexes_T_O`:      Decomposition for the technological component, outputs
    - `prodIndexes_T_G_O`:    Decomposition for the technological component, good outputs
    - `prodIndexes_T_B_O`:    Decomposition for the technological component, bad outputs
    - `prodIndexes_T_I`:      Decomposition for the technological component, inputs
    - `prodIndexes_T_G_I`:    Decomposition for the technological component, good inputs
    - `prodIndexes_T_B_I`:    Decomposition for the technological component, bad inputs
    - `prodIndexes_E`:        Decomposition for the efficiency component, overall
    - `prodIndexes_E_O`:      Decomposition for the efficiency component, outputs
    - `prodIndexes_E_G_O`:    Decomposition for the efficiency component, good output
    - `prodIndexes_E_B_O`:    Decomposition for the efficiency component, bad outputs
    - `prodIndexes_E_I`:      Decomposition for the efficiency component, inputs
    - `prodIndexes_E_G_I`:    Decomposition for the efficiency component, good inputs
    - `prodIndexes_E_B_I`:    Decomposition for the efficiency component, bad inputs
    - `prodIndexes_S`:        Decomposition for the scale (residual) component, overall
    - `prodIndexes_S_O`:      Decomposition for the scale (residual) component, outputs
    - `prodIndexes_S_G_O`:    Decomposition for the scale (residual) component, good output
    - `prodIndexes_S_B_O`:    Decomposition for the scale (residual) component, bad outputs
    - `prodIndexes_S_I`:      Decomposition for the scale (residual) component, inputs
    - `prodIndexes_S_G_I`:    Decomposition for the scale (residual) component, good inputs
    - `prodIndexes_S_B_I`:    Decomposition for the scale (residual) component, bad inputs


 The second and third element of the tuple are respectively the "good inputs/outputs" and "bad/inputs/outputs" components. The first matrix can be retrieved from the two components by multiplying them (for `multiplicative` production structure) or summing them (for `additive` production strucure).


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
                   retToScale="constant",prodStructure="multiplicative",convexAssumption=true)

    nDMUs = size(gI,1)
    nPer = size(gI,3)
    nBI  = size(bI,2)

    # Output containers

    # - Global
    prodIndexes = Array{Union{Float64,Missing}}(undef, (nDMUs,nPer-1))

    # - decomposition by type of the inputs or outputs (Good vs Bad)
    prodIndexes_G = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_B = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))

    # - decomposition by tecnological/efficiency/scale (residual)
    prodIndexes_T     = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1)) # technological
    prodIndexes_T_O   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_T_G_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_T_B_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_T_I   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_T_G_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_T_B_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))

    prodIndexes_E     = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1)) # Efficiency
    prodIndexes_E_O   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_E_G_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_E_B_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_E_I   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_E_G_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_E_B_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))

    #prodIndexes_S     = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1)) # Scale (residual)
    prodIndexes_S_O   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S_G_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S_B_O = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S_I   = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S_G_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S_B_I = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))
    prodIndexes_S     = Array{Union{Float64,Missing}}(undef, (nDMUs,nPer-1))



    noBadIndexDefault = prodStructure == "multiplicative" ? 1.0 : 0.0

    for t in 1:nPer-1
        #println(t)
        gIₜ = gI[:,:,t]
        gIᵤ = gI[:,:,t+1]
        bIₜ = bI[:,:,t]
        bIᵤ = bI[:,:,t+1]
        gOₜ = gO[:,:,t]
        gOᵤ = gO[:,:,t+1]
        bOₜ = bO[:,:,t]
        bOᵤ = bO[:,:,t+1]
        #if convexAssumption == true
            dirGI = prodStructure == "multiplicative" ? (-1,0,0,0) : (1,0,0,0)
            dirBI = prodStructure == "multiplicative" ? (0,-1,0,0) : (0,1,0,0)
            dirGO = prodStructure == "multiplicative" ? (0,0,1,0)  : (0,0,1,0)
            dirBO = prodStructure == "multiplicative" ? (0,0,0,-1) : (0,0,0,1)
        #else
        #    dirGI = prodStructure == "multiplicative" ? (1,0,0,0)  : (1,0,0,0)
        #    dirBI = prodStructure == "multiplicative" ? (0,1,0,0)  : (0,1,0,0)
        #    dirGO = prodStructure == "multiplicative" ? (0,0,-1,0) : (0,0,1,0)
        #    dirBO = prodStructure == "multiplicative" ? (0,0,0,1)  : (0,0,0,-1)
        #end
        for z in 1:nDMUs
            gIₜ₀ = gIₜ[z,:]
            bIₜ₀ = bIₜ[z,:]
            gOₜ₀ = gOₜ[z,:]
            bOₜ₀ = bOₜ[z,:]
            gIᵤ₀ = gIᵤ[z,:]
            bIᵤ₀ = bIᵤ[z,:]
            gOᵤ₀ = gOᵤ[z,:]
            bOᵤ₀ = bOᵤ[z,:]

           #if z == 1 && t == nPer-1
            #   println("here we stop and check")
           #end



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

            # DMU at t, frontier at t+1...
            idx_gi_tu = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI,crossTime=true)
            idx_bi_tu = (nBI == 0) ? noBadIndexDefault : problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI,crossTime=true)
            idx_go_tu = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO,crossTime=true)
            idx_bo_tu = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO,crossTime=true)

            # DMU at t+1, frontier at t...
            idx_gi_ut = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGI,crossTime=true)
            idx_bi_ut = (nBI == 0) ? noBadIndexDefault : problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBI,crossTime=true)
            idx_go_ut = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirGO,crossTime=true)
            idx_bo_ut = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,forceLinearModel=forceLinearModel,
                        directions=dirBO,crossTime=true)

            # Computing aggregated for the full (gI,Bi,bO,gO) partition...
            if prodStructure == "multiplicative"
                if convexAssumption == true
                    idx_i_t = (idx_gi_t̃/idx_gi_t) * (idx_bi_t̃/idx_bi_t)
                    idx_o_t = (idx_go_t/idx_go_t̃) * (idx_bo_t/idx_bo_t̃)
                    idx_t   = idx_o_t/idx_i_t
                    idx_i_u = (idx_gi_u/idx_gi_ũ) * (idx_bi_u/idx_bi_ũ)
                    idx_o_u = (idx_go_ũ/idx_go_u) * (idx_bo_ũ/idx_bo_u)
                    idx_u   = idx_o_u/idx_i_u
                    idx     = (idx_t * idx_u)^(1/2)
                else
                    idx_i_t = (idx_gi_t/idx_gi_t̃) * (idx_bi_t/idx_bi_t̃)
                    idx_o_t = (idx_go_t/idx_go_t̃) / (idx_bo_t̃/idx_bo_t)
                    idx_t   = idx_o_t/idx_i_t
                    idx_i_u = (idx_gi_ũ/idx_gi_u) * (idx_bi_ũ/idx_bi_u)
                    idx_o_u = (idx_go_ũ/idx_go_u) / (idx_bo_u/idx_bo_ũ)
                    idx_u   = idx_o_u/idx_i_u
                    idx     = (idx_t * idx_u)^(1/2)
                end


            else
                if convexAssumption == true
                    idx_i_t = (idx_gi_t̃ - idx_gi_t) + (idx_bi_t̃ - idx_bi_t)
                    idx_o_t = (idx_go_t - idx_go_t̃) - (idx_bo_t̃ - idx_bo_t)
                    idx_t   = idx_o_t - idx_i_t
                    idx_i_u = (idx_gi_u - idx_gi_ũ) + (idx_bi_u - idx_bi_ũ)
                    idx_o_u = (idx_go_ũ - idx_go_u) - (idx_bo_u - idx_bo_ũ)
                    idx_u   = idx_o_u - idx_i_u
                    idx     = (idx_t + idx_u) / 2
                else

                end
            end
            prodIndexes[z,t] = idx

            # Computing aggregated for the good inputs/outputs only...
            if prodStructure == "multiplicative"
               idx_Gi_t = (idx_gi_t̃/idx_gi_t)
               idx_Go_t = (idx_go_t/idx_go_t̃)
               idx_Gt   = idx_Go_t/idx_Gi_t
               idx_Gi_u = (idx_gi_u/idx_gi_ũ)
               idx_Go_u = (idx_go_ũ/idx_go_u)
               idx_Gu   = idx_Go_u/idx_Gi_u
               idx_G    = (idx_Gt * idx_Gu)^(1/2)
            else
                idx_Gi_t = (idx_gi_t̃ - idx_gi_t)
                idx_Go_t = (idx_go_t - idx_go_t̃)
                idx_Gt   = idx_Go_t - idx_Gi_t
                idx_Gi_u = (idx_gi_u - idx_gi_ũ)
                idx_Go_u = (idx_go_ũ - idx_go_u)
                idx_Gu   = idx_Go_u - idx_Gi_u
                idx_G     = (idx_Gt + idx_Gu) / 2
            end
            prodIndexes_G[z,t] = idx_G

            # Computing aggregated for the bad inputs/outputs only...
            if prodStructure == "multiplicative"
               idx_Bi_t = (idx_bi_t̃/idx_bi_t)
               idx_Bo_t = (idx_bo_t̃/idx_bo_t)
               idx_Bt   = 1/(idx_Bo_t * idx_Bi_t)
               idx_Bi_u = (idx_bi_u/idx_bi_ũ)
               idx_Bo_u = (idx_bo_u/idx_bo_ũ)
               idx_Bu   = 1/(idx_Bo_u*idx_Bi_u)
               idx_B    = (idx_Bt * idx_Bu)^(1/2)
            else
                idx_Bi_t = (idx_bi_t̃ - idx_bi_t)
                idx_Bo_t = (idx_bo_t̃ - idx_bo_t)
                idx_Bt   = - idx_Bo_t - idx_Bi_t
                idx_Bi_u = (idx_bi_u - idx_bi_ũ)
                idx_Bo_u = (idx_bo_u - idx_bo_ũ)
                idx_Bu   = - idx_Bo_u - idx_Bi_u
                idx_B    = (idx_Bt + idx_Bu) / 2
            end
            prodIndexes_B[z,t] = idx_B

            # Computing tecnological/efficiency/residual decomposition
            if prodStructure == "multiplicative"
                # Disaggregation T/E/S...
                idx_T_G_O = ((idx_go_t/idx_go_tu) * (idx_go_ut/idx_go_u) )^-(1/2)
                idx_T_B_O = ((idx_bo_t/idx_bo_tu) * (idx_bo_ut/idx_bo_u) )^-(1/2)
                idx_T_O   = idx_T_G_O * idx_T_B_O
                idx_T_G_I = ((idx_gi_t/idx_gi_tu) * (idx_gi_ut/idx_gi_u) )^-(1/2)
                idx_T_B_I = ((idx_bi_t/idx_bi_tu) * (idx_bi_ut/idx_bi_u) )^-(1/2)
                idx_T_I   = idx_T_G_I * idx_T_B_I
                idx_T     = (idx_T_O * idx_T_I)

                idx_E_G_O = (idx_go_u/idx_go_t)^-1
                idx_E_B_O = (idx_bo_u/idx_bo_t)^-1
                idx_E_O   = idx_E_G_O * idx_E_B_O
                idx_E_G_I = (idx_gi_u/idx_gi_t)^-1
                idx_E_B_I = (idx_bi_u/idx_bi_t)^-1
                idx_E_I   = idx_E_G_I * idx_E_B_I
                idx_E     = (idx_E_O * idx_E_I)

                idx_S_I   = idx / (idx_T_I * idx_E_I)
                idx_S_O   = idx / (idx_T_O * idx_E_O)
                idx_S_G_O = ((idx_go_t̃ / idx_go_ut) * (idx_gi_t̃ / idx_gi_t) * (idx_go_tu / idx_go_ũ) * (idx_gi_u/idx_gi_ũ)  )^-(1/2)
                idx_S_B_O = ((idx_bo_t̃ / idx_bo_ut) * (idx_bi_t̃ / idx_bi_t) * (idx_bo_tu / idx_bo_ũ) * (idx_bi_u/idx_bi_ũ)  )^-(1/2)
                idx_S_G_I = ((idx_gi_t̃ / idx_gi_ut) * (idx_go_t̃ / idx_go_t) * (idx_gi_tu / idx_gi_ũ) * (idx_go_u/idx_go_ũ)  )^-(1/2)
                idx_S_B_I = ((idx_bi_t̃ / idx_bi_ut) * (idx_bo_t̃ / idx_bo_t) * (idx_bi_tu / idx_bi_ũ) * (idx_bo_u/idx_bo_ũ)  )^-(1/2)

                idx_S     = idx / (idx_T * idx_E)


            else # additive production structure
                # Disaggregation T/E/S...
                idx_T_G_O = ((-idx_go_t+idx_go_tu) + (-idx_go_ut+idx_go_u) )/2
                idx_T_B_O = ((-idx_bo_t-idx_bo_tu) + (+idx_bo_ut+idx_bo_u) )/2
                idx_T_O   = idx_T_G_O + idx_T_B_O
                idx_T_G_I = ((-idx_gi_t+idx_gi_tu) + (-idx_gi_ut+idx_gi_u) )/2
                idx_T_B_I = ((-idx_bi_t+idx_bi_tu) + (-idx_bi_ut+idx_bi_u) )/2
                idx_T_I   = idx_T_G_I + idx_T_B_I
                idx_T     = (idx_T_O + idx_T_I)

                idx_E_G_O = -(idx_go_u-idx_go_t)
                idx_E_B_O = -(idx_bo_u-idx_bo_t)
                idx_E_O   = idx_E_G_O + idx_E_B_O
                idx_E_G_I = -(idx_gi_u-idx_gi_t)
                idx_E_B_I = -(idx_bi_u-idx_bi_t)
                idx_E_I   = idx_E_G_I + idx_E_B_I
                idx_E     = (idx_E_O + idx_E_I)

                idx_S_I   = idx - (idx_T_I + idx_E_I)
                idx_S_O   = idx - (idx_T_O + idx_E_O)
                idx_S_G_O = -((idx_go_t̃ - idx_go_ut) + (idx_gi_t̃ - idx_gi_t) + (idx_go_tu - idx_go_ũ) + (idx_gi_u-idx_gi_ũ)  )/2
                idx_S_B_O = -((idx_bo_t̃ + idx_bo_ut) + (idx_bi_t̃ - idx_bi_t) + (-idx_bo_tu - idx_bo_ũ) + (idx_bi_u-idx_bi_ũ)  )/2
                idx_S_G_I = -((idx_gi_t̃ - idx_gi_ut) + (idx_go_t̃ - idx_go_t) + (idx_gi_tu - idx_gi_ũ) + (idx_go_u-idx_go_ũ)  )/2
                idx_S_B_I = -((idx_bi_t̃ - idx_bi_ut) + (idx_bo_t̃ - idx_bo_t) + (idx_bi_tu - idx_bi_ũ) + (idx_bo_u-idx_bo_ũ)  )/2
                idx_S     = idx - (idx_T + idx_E)

                #=
                idx_T_G_O = ((idx_go_t-idx_go_ũ) + (idx_go_t̃-idx_go_u) )/2
                idx_T_B_O = ((idx_bo_t-idx_bo_ũ) + (idx_bo_t̃-idx_bo_u) )/2
                idx_T_O   = idx_T_G_O + idx_T_B_O
                idx_T_G_I = ((idx_gi_t-idx_gi_ũ) + (idx_gi_t̃-idx_gi_u) )/2
                idx_T_B_I = ((idx_bi_t-idx_bi_ũ) + (idx_bi_t̃-idx_bi_u) )/2
                idx_T_I   = idx_T_G_I + idx_T_B_I
                idx_T     = (idx_T_O + idx_T_I)/2

                idx_E_G_O = idx_go_u-idx_go_t
                idx_E_B_O = idx_bo_u-idx_bo_t
                idx_E_O   = idx_E_G_O + idx_E_B_O
                idx_E_G_I = idx_gi_u-idx_gi_t
                idx_E_B_I = idx_bi_u-idx_bi_t
                idx_E_I   = idx_E_G_I + idx_E_B_I
                idx_E     = (idx_E_O + idx_E_I)/2

                idx_S_G_O = idx_G -  (idx_T_G_O + idx_E_G_O)
                idx_S_B_O = idx_B -  (idx_T_B_O + idx_E_B_O)
                idx_S_O   = idx_S_G_O + idx_S_B_O
                idx_S_G_I = idx_G -  (idx_T_G_I + idx_E_G_I)
                idx_S_B_I = idx_B -  (idx_T_B_I + idx_E_B_I)
                idx_S_I   = idx_S_G_I + idx_S_B_I
                #idx_S     = (idx_S_O + idx_S_I)/2
                =#
            end

            prodIndexes_T[z,t]      =  idx_T
            prodIndexes_T_O[z,t]    =  idx_T_O
            prodIndexes_T_G_O[z,t]  =  idx_T_G_O
            prodIndexes_T_B_O[z,t]  =  idx_T_B_O
            prodIndexes_T_I[z,t]    =  idx_T_I
            prodIndexes_T_G_I[z,t]  =  idx_T_G_I
            prodIndexes_T_B_I[z,t]  =  idx_T_B_I
            prodIndexes_E[z,t]      =  idx_E
            prodIndexes_E_O[z,t]    =  idx_E_O
            prodIndexes_E_G_O[z,t]  =  idx_E_G_O
            prodIndexes_E_B_O[z,t]  =  idx_E_B_O
            prodIndexes_E_I[z,t]    =  idx_E_I
            prodIndexes_E_G_I[z,t]  =  idx_E_G_I
            prodIndexes_E_B_I[z,t]  =  idx_E_B_I
            prodIndexes_S[z,t]      =  idx_S
            prodIndexes_S_O[z,t]    =  idx_S_O
            prodIndexes_S_G_O[z,t]  =  idx_S_G_O
            prodIndexes_S_B_O[z,t]  =  idx_S_B_O
            prodIndexes_S_I[z,t]    =  idx_S_I
            prodIndexes_S_G_I[z,t]  =  idx_S_G_I
            prodIndexes_S_B_I[z,t]  =  idx_S_B_I

        end # end for each DMU
    end # end for each period

    return (prodIndexes        =  prodIndexes       ,
            prodIndexes_G      =  prodIndexes_G     ,
            prodIndexes_B      =  prodIndexes_B     ,
            prodIndexes_T      =  prodIndexes_T     ,
            prodIndexes_T_O    =  prodIndexes_T_O   ,
            prodIndexes_T_G_O  =  prodIndexes_T_G_O ,
            prodIndexes_T_B_O  =  prodIndexes_T_B_O ,
            prodIndexes_T_I    =  prodIndexes_T_I   ,
            prodIndexes_T_G_I  =  prodIndexes_T_G_I ,
            prodIndexes_T_B_I  =  prodIndexes_T_B_I ,
            prodIndexes_E      =  prodIndexes_E     ,
            prodIndexes_E_O    =  prodIndexes_E_O   ,
            prodIndexes_E_G_O  =  prodIndexes_E_G_O ,
            prodIndexes_E_B_O  =  prodIndexes_E_B_O ,
            prodIndexes_E_I    =  prodIndexes_E_I   ,
            prodIndexes_E_G_I  =  prodIndexes_E_G_I ,
            prodIndexes_E_B_I  =  prodIndexes_E_B_I ,
            prodIndexes_S      =  prodIndexes_S     ,
            prodIndexes_S_O    =  prodIndexes_S_O   ,
            prodIndexes_S_G_O  =  prodIndexes_S_G_O ,
            prodIndexes_S_B_O  =  prodIndexes_S_B_O ,
            prodIndexes_S_I    =  prodIndexes_S_I   ,
            prodIndexes_S_G_I  =  prodIndexes_S_G_I ,
            prodIndexes_S_B_I  =  prodIndexes_S_B_I )
end
