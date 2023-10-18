###############################################################
# Part of BDisposal.jl package
###############################################################


"""
    efficiencyScores(gI,gO,bO,bI;retToScale,prodStructure,dirGI,dirBI,dirGO,dirBO,startθ,startμ,startλ)

Compute efficiency indicators and convexity test results considering the joined production of good and bad outputs.

Given a set of measures of inputs, "good" ("desiderable") and "bad" ("undesiderable") outputs for different decision making units, compute their distance to the production frontier, i.e. their degree of efficiency.

## Parameters:
- Positional
  *  `gI`: "Good" inputs (3D matrix by DMUs, input items and periods)
  *  `gO`: "Good" outputs (3D matrix by DMUs, input items and periods)
  *  `bO`: "Bad" outputs (3D matrix by DMUs, input items and periods)
  *  `bI`: "Bad" inputs (optional 3D matrix by DMUs, input items and periods) [default: empty array]
- Keyword
  *  `retToScale`: Wheter the returns to scale should be assumed "constant" or "variable" (default). Non-convex distance and test is computed only under the "variable" assumption, except for the `(0,0,1,-1)` distance under multiplicative production function
  *  `prodStructure`: Wheter the production structure should be assumed "additive" (default) or "multiplicative"
  *  `dirGI`,`dirBI`,`dirGO`,`dirBO`: The directions toward where to measure the efficiency (see notes) [default: `(0,0,1,0)`]
  *  `startθ`,`startμ`,`startλ`: Initial values in the convex optimisation problem [default: `(0,0,1.1)`]

## Returns:
- A named tuple with the following items:
  * `λs`: The efficiency indicators coherent with the convexity test (2D matrix by DMUs and periods)
  * `λs_convex`: The efficiency indicators assuming a convex production frontier (2D matrix by DMUs and periods)
  * `λs_nonconvex`:  The efficiency indicators assuming a non-convex production frontier (2D matrix by DMUs and periods)
  * `nonConvTest`: The result of the non-convexity test (2D matrix of boolean by DMUs and periods). Note that "true" here refers to non-convex.
  * `nonConvTest_value`: The value of the non-convexity test (2D matrix by DMUs and periods)

## Description of the function

`efficiencyScores` defines environmental efficiency indicators for a set of Decision Making Units (ie., observations). The nature of the DMUs (eg. firms, countries etc.) depends on the aim of the environmental efficiency analysis.
    
All the DMUs are characterized by a set of inputs and outputs separated into undesirable (eg., polluting) and desirable (eg., no polluting) ones. Note that, inputs separation into desirable (eg., non emissions-causing) and undesirable (eg., emissions-causing) components is optional for the BDisposal package.

The BDisposal DEA modelling defines the smallest environmental production process that contains all DMUs informations (inputs/outputs) and satisfies the B-disposal axiomatic pattern.

The BDisposal algorithms for efficiency assessment compute the distance to the best
environmental production procedures; ie., to the efficient production frontier. These
efficiency indices inherit the structure of additive and multiplicative distance functions.

Notice that the BDisposal package computes general shape of additive and multiplicative distance functions. It follows that, the Bdisposal package introduces a flexible framework for environmental efficiency assessment, with the possibility to analyse different directions for the distance functions.
In addition, the BDisposal package allows to test the (economists’) properties of convexity and disposability.

## Interpretation of the results

### Efficiency indicators

When the additive shape of the distance function is greater than 0 then, the production unit doesn’t operate efficiently.
If the additive efficiency score is equal to 0 then, the DMU is efficient. It follows that, the DMU employs one of the best production technology to transform inputs into outputs; ie., the producer is a benchmark for the other observations. A similar reasoning holds for the multiplicative distance functions.

If the multiplicative efficiency index is greater than 1 then, the DMU doesn’t perform
efficiently. When the multiplicative efficiency score is equal to 1 then, the production
unit is a benchmark for the other observations.

### Convexity test

When the multiplicative test of convexity is greater than 1 (respectively, equal to 1) then, the production process frontier is non convex (respectively, convex). If the multiplicative test of disposability is greater than 1 (respectively, equal to 1) then, the frontier of the production set displays B-disposability; ie., there exists costs for undesirable component decrease. A similar reasoning applies for the additive framework.

The informations provided by these tests define the shape of the production process boundaries.

## Notes:
* Directions toward where to measure the efficiency distance can be tuned using the `dirGI`,`dirBI`,`dirGO` and `dirBO` parameters.
* The following directions are supported:
  - multiplicative pr. struct.: (-1,0,0,0), (0,-1,0,0), (0,0,1,0), (0,0,0,-1) or (0,0,1,-1)
  - addittive pr. struct.: (1,0,0,0), (0,1,0,0), (0,0,1,0) or (0,0,0,1)
  Other directions can be used to compute the direction under the convex frontier assumption, but the convexity test is not performed and, for the multiplicative case, no guarantee is given on solving the underliying (non-linear) problem (different initial value `startθ`,`startμ`,`startλ` can be attempted)
"""
function efficiencyScores(gI::AbstractArray{TgI,3},gO::AbstractArray{TgO,3},bO::AbstractArray{TbO,3},bI::AbstractArray{TbI,3}=zeros(size(gI,1),0,size(gI,3));
                          retToScale="variable",prodStructure="additive",
                          dirGI=0,dirBI=0,dirGO=1,dirBO=0,
                          startθ=0.01,startμ=0.01,startλ=1.1) where {TgI <: Number, TgO <: Number, TbO <: Number, TbI <: Number}

    # Data processing
    nDMUs, nPer, nGI, nBI, nGO, nBO = size(gI,1), size(gI,3), size(gI,2), size(bI,2), size(gO,2), size(bO,2)
    λs_convex    = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)
    λs_nonconvex = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)
    directions   = (dirGI,dirBI,dirGO,dirBO)


    # Tests if doing also convex analysis or not
    doNonConvexAnalysis = true

    if retToScale == "constant" && directions != (0,0,1,-1)
        @warn "Non-convex analysis is supported only with variable returns to scale or with direcitons (0,0,1,-1) and multiplicative production structure."
        doNonConvexAnalysis = false
    end

    if  prodStructure == "multiplicative"
        if !( directions in [(-1,0,0,0),(0,-1,0,0),(0,0,1,0),(0,0,0,-1),(0,0,1,-1)])
            @warn "Non-convex analysis not supported with these directions. For the multiplicative case,
            non-convex analysis is supported only with directions (-1,0,0,0),(0,-1,0,0),(0,0,1,0),(0,0,0,-1) or (0,0,1,-1)"
             doNonConvexAnalysis = false
        end
    else
        if !( directions in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)])
            @warn "Non-convex analysis not supported with these directions. For the addittive case,
            non-convex analysis is supported only with directions (1,0,0,0),(0,1,0,0),(0,0,1,0) or (0,0,0,1)"
             doNonConvexAnalysis = false
        end
    end
    # temp todo:
    doNonConvexAnalysis = true

    # Loping over the periods
    for t in 1:nPer
        gIₜ = gI[:,:,t]; gOₜ = gO[:,:,t]; bOₜ = bO[:,:,t]
        #bIₜ = nBI > 0 ? bI[:,:,t] : Array{Float64}(undef, 0,0)
        bIₜ = bI[:,:,t]

        #println("Period: $t")
        # Solving for each dmu
        for z in 1:nDMUs
            gIₜ₀ = gIₜ[z,:]; gOₜ₀ = gOₜ[z,:]; bOₜ₀ = bOₜ[z,:]
            #bIₜ₀ = nBI > 0 ? bIₜ[z,:] : Float64[]
            bIₜ₀ = bIₜ[z,:]

        #    println("DMU: $z")


            # CONVEX analysis....
            linearModelPossibleDirections = [(-1,0,0,0),(0,-1,0,0),(0,0,1,0),(0,0,0,-1)]
            forceLinearModel = (dirGI,dirBI,dirGO,dirBO) in linearModelPossibleDirections ?  true : false

            λs_convex[z,t] = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ;
                                     retToScale=retToScale,prodStructure=prodStructure,
                                     convexAssumption=true,
                                     directions=(dirGI,dirBI,dirGO,dirBO),
                                     startValues=(startθ,startμ,startλ),
                                     forceLinearModel=forceLinearModel)

            # Non-convex analysis (when implemented)...
            λs_nonconvex[z,t] = doNonConvexAnalysis ? problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ;
                                        retToScale=retToScale,prodStructure=prodStructure,
                                        directions=(dirGI,dirBI,dirGO,dirBO),
                                        convexAssumption=false) : missing
        end # end of each dmu
    end # endof each period

#    println("Done all DMU problems")
    nonConvTest_value = λs_nonconvex ./ λs_convex
    nonConvTest       = nonConvTest_value .< (1.0 - eps())
    λs                = doNonConvexAnalysis ? min.(λs_nonconvex,λs_convex) : λs_convex

    return λs, λs_convex, λs_nonconvex, nonConvTest, nonConvTest_value
end
