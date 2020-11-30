module BDisposal

using JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex

# Define the function that compute the efficiency measures
"""
  efficiencyScores(inputs,goodOutputs,badOutputs,data;retToScale="constant",effAnalysisType="static",formattedOutput=false)

Compute efficiency scores considering the joined production of good and bad outputs.

Given a set of measures of inputs, "good" ("desiderable") and "bad" ("undesiderable") outputs for different decision making units, compute their
distance to the production frontier, i.e. their degree of efficiency.

# Parameters:
* `inputs`: The headers of the inputs (vector of strings)
* `goodOutputs`: The headers of the good putputs (vector of strings)
* `badOutputs`: The headers of the bad putputs (vector of strings)
* `data`: A dataframe where the first column must have header `period`, the second one `dmu`
   (standing for "decision making unit") and then the measures of inputs, good and bad outputs
   with the headers as defined in the first 3 parameters
* `retToScale`: Wheter to return to scales should be assumed `constant` (default) or `variable`
* `effAnalysisType`: Wheter to perform a `static` (default) or `dynamic` efficiency analysis
* `formattedOutput`: wheater the output should be formatted or not (i.e. raw output) [def: true]
# Returns:
* A printable matrix of the computed efficiency scores
# Notes:
* Only the static efficiency analysis is implemented at the moment
"""
#function efficiencyScores(inputs,goodOutputs,badOutputs,data;
#                                   retToScale="constant",formattedOutput=false,
#                                   badInputs=[],prodStructure="additive",
#                                   dirGI=0,dirBI=0,dirGO=1,dirBO=-1,
#                                   startθ=0,startμ=0,startλ=1.1)


function efficiencyScores(gI::Array{Float64,3},gO::Array{Float64,3},bO::Array{Float64,3},bI::Array{Float64,3}=Array{Float64}(undef, (size(gI,1),0,size(gI,3)));
                          retToScale="constant",prodStructure="additive",
                          dirGI=0,dirBI=0,dirGO=1,dirBO=-1,
                          startθ=0,startμ=0,startλ=1.1)

    # Data processing
    nDMUs, nPer, nGI, nBI, nGO, nBO = size(gI,1), size(gI,3), size(gI,2), size(bI,2), size(gO,2), size(bO,2)
    λs_convex    = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)
    λs_nonconvex = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)

    # Loping over the periods
    for t in 1:nPer
        gIₜ = gI[:,:,t]; gOₜ = gO[:,:,t]; bOₜ = bO[:,:,t]
        #bIₜ = nBI > 0 ? bI[:,:,t] : Array{Float64}(undef, 0,0)
        bIₜ = bI[:,:,t]

        # Solving for each dmu
        for z in 1:nDMUs
            gIₜ₀ = gIₜ[z,:]; gOₜ₀ = gOₜ[z,:]; bOₜ₀ = bOₜ[z,:]
            #bIₜ₀ = nBI > 0 ? bIₜ[z,:] : Float64[]
            bIₜ₀ = bIₜ[z,:]

            # CONVEX analysis....
            λs_convex[z,t] = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ;
                                     retToScale=retToScale,prodStructure=prodStructure,
                                     convexAssumption=true,
                                     directions=(dirGI,dirBI,dirGO,dirBO),
                                     startValues=(startθ,startμ,startλ))

            # NON Convex analysis
            λs_nonconvex[z,t] = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ;
                                        retToScale=retToScale,prodStructure=prodStructure,
                                        directions=(dirGI,dirBI,dirGO,dirBO),
                                        convexAssumption=false)
        end # end of each dmu
    end # endof each period

    nonConvTest_value = λs_nonconvex ./ λs_convex
    nonConvTest       = nonConvTest_value .< (1.0 - eps())
    λs = min.(λs_nonconvex,λs_convex)

    return λs, λs_convex, λs_nonconvex, nonConvTest_value, nonConvTest
end


function prodIndex(gI::Array{Float64,3},gO::Array{Float64,3},bO::Array{Float64,3},bI::Array{Float64,3}=Array{Float64}(undef, (size(gI,1),0,size(gI,3)));
                   retToScale="constant",prodStructure="multiplicative",convexAssumption=true,
                   startθ=0,startμ=0,startλ=1.1)

    nDMUs = size(gI,1)
    nPer = size(gI,3)
    startValues = (startθ,startμ,startλ)
    prodIndexes = Array{Union{Float64,Missing}}(undef, (size(gI,1),nPer-1))

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

            # t...
            idx_gi_t̃ = problem(gIᵤ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGI,startValues=startValues,crossTime=true)
            idx_gi_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGI,startValues=startValues)
            idx_bi_t̃ = problem(gIₜ₀,bIᵤ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBI,startValues=startValues,crossTime=true)
            idx_bi_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBI,startValues=startValues)
            idx_go_t̃ = problem(gIₜ₀,bIₜ₀,gOᵤ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGO,startValues=startValues,crossTime=true)
            idx_go_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGO,startValues=startValues)
            idx_bo_t̃ = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOᵤ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBO,startValues=startValues,crossTime=true)
            idx_bo_t = problem(gIₜ₀,bIₜ₀,gOₜ₀,bOₜ₀,gIₜ,bIₜ,gOₜ,bOₜ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBO,startValues=startValues)

            # u...
            idx_gi_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGI,startValues=startValues)
            idx_gi_ũ = problem(gIₜ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGI,startValues=startValues,crossTime=true)
            idx_bi_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBI,startValues=startValues)
            idx_bi_ũ = problem(gIᵤ₀,bIₜ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBI,startValues=startValues,crossTime=true)
            idx_go_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGO,startValues=startValues)
            idx_go_ũ = problem(gIᵤ₀,bIᵤ₀,gOₜ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirGO,startValues=startValues,crossTime=true)
            idx_bo_u = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOᵤ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBO,startValues=startValues)
            idx_bo_ũ = problem(gIᵤ₀,bIᵤ₀,gOᵤ₀,bOₜ₀,gIᵤ,bIᵤ,gOᵤ,bOᵤ,
                        retToScale=retToScale,prodStructure=prodStructure,
                        convexAssumption=convexAssumption,
                        directions=dirBO,startValues=startValues,crossTime=true)

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

# Dispatching to either convexProblem or nonConvexProblem
function problem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                 retToScale,prodStructure,convexAssumption=true,
                 directions=(),
                 startValues=(),crossTime=false)

     if convexAssumption
         return convexProblem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                          retToScale=retToScale,prodStructure=prodStructure,
                          directions=directions,
                          startValues=startValues,crossTime=crossTime)
     else
         return nonConvexProblem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                          retToScale=retToScale,prodStructure=prodStructure,
                          directions=directions)
     end

end

function convexProblem(inp₀,bInp₀,gO₀,bO₀,inp,bInp,gO,bO;
                       retToScale,prodStructure,
                       directions,
                       startValues,crossTime=false)

   nI, nGO, nBO, nDMUs, nbI = length(inp₀), length(gO₀), length(bO₀), size(inp,1), length(bInp₀)
   (dirGI,dirBI,dirGO,dirBO) = directions
   (startθ,startμ,startλ)    = startValues
   if prodStructure == "multiplicative"
       Ipopt.amplexe() do path
           # Model declaration (efficiency model)
           effmodel = Model(() -> AmplNLWriter.Optimizer(path, ["print_level=0"]))
           #effmodel = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
           #effmodel = Model(() -> AmplNLWriter.Optimizer("ipopt",["print_level=0"]))

           # Defining variables
           @variables effmodel begin
               θ[z in 1:nDMUs] >= 0, (start = startθ) # thetas (nDMUs)
               μ[z in 1:nDMUs] >= 0, (start = startμ) # mus (nDMUs)
           end
           if !crossTime
              @variable(effmodel, λ >= 1, start = startλ)
           else
              @variable(effmodel, λ, start = startλ)
           end

           # Defining constraints
           @NLconstraints effmodel begin
               cgInp1[i in 1:nI],   # constraint on good inputs
                  λ^dirGI * inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
               cbInp1[i in 1:nbI], # constraint on bad inputs
                  λ^dirBI * bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
               cgInp2[i in 1:nI],   # second constraint on good inputs
                  λ^dirGI * inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
               cbInp2[i in 1:nbI], # second constraint on bad inputs
                  λ^dirBI * bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
               cgO1[j in 1:nGO],   # first constraint on good outputs
                  λ^dirGO * gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
               cbO1[j in 1:nBO],   # first constraint on bad outputs
                  λ^dirBO * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
               cgO2[j in 1:nGO],   # second constraint on good outputs
                  λ^dirGO * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
               cbO2[j in 1:nBO],   # second constraint on bad outputs
                  λ^dirBO * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
           end
           if retToScale == "variable"
               @constraints effmodel begin
                   sumtheta,
                     sum(θ[z] for z in 1:1:nDMUs) == 1
                   summu,
                     sum(μ[z] for z in 1:1:nDMUs) == 1
               end
           end
           # Defining objective
           @objective effmodel Max begin
               λ
           end
           # Printing the model (for debugging)
           #print(effmodel) # The model in mathematical terms is printed

           # Solving the model
           oldstd = stdout
           #redirect_stdout(open("/dev/null", "w"))
           redirect_stdout(open("nul", "w"))
           #open("nul", "w")
           JuMP.optimize!(effmodel)
           #redirect_stdout((()->optimize!(effmodel)),open("/dev/null", "w"))
           redirect_stdout(oldstd) # recover original stdout
           # Copy the results to the output matrix...
           status = termination_status(effmodel)

           #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
           if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
               return value(λ)
           else
               return missing
           end
       end # end of do function
   else # prodStructure is then addictive...
       # Model declaration (efficiency model)
       effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
       # Defining variables
       @variables effmodel begin
           θ[z in 1:nDMUs] >= 0 # thetas (nDMUs)
           μ[z in 1:nDMUs] >= 0 # mus (nDMUs)
       end
       if !crossTime
          @variable(effmodel, λ >= 0)
       else
          @variable(effmodel, λ)
       end
       # Defining constraints
       # additive structure
       @constraints effmodel begin
           cgInp1[i in 1:nI],   # constraint on good inputs
              inp₀[i]- λ*dirGI * inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
           cbInp1[i in 1:nbI], # constraint on bad inputs
              bInp₀[i]- λ*dirBI * bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
           cgInp2[i in 1:nI],   # second constraint on good inputs
               inp₀[i] - λ*dirGI * inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
           cbInp2[i in 1:nbI], # second constraint on bad inputs
              bInp₀[i] - λ*dirBI * bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
           cgO1[j in 1:nGO],   # first constraint on good outputs
              gO₀[j] + λ * dirGO * gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
           cbO1[j in 1:nBO],   # first constraint on bad outputs
              bO₀[j] + λ * dirBO * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
           cgO2[j in 1:nGO],   # second constraint on good outputs
              gO₀[j] + λ*dirGO * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
           cbO2[j in 1:nBO],   # second constraint on bad outputs
              bO₀[j] + λ*dirBO * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
       end
       if retToScale == "variable"
           @constraints effmodel begin
               sumtheta,
                 sum(θ[z] for z in 1:1:nDMUs) == 1
               summu,
                 sum(μ[z] for z in 1:1:nDMUs) == 1
           end
       end
       # Defining objective
       @objective effmodel Max begin
           λ
       end
       # Printing the model (for debugging)
       #print(effmodel) # The model in mathematical terms is printed

       # Solving the model
       optimize!(effmodel)

       # Copy the results to the output matrix...
       status = termination_status(effmodel)

       #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
       if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
           return value(λ)
       else
           return missing
       end
   end # ending addictive  prodStructure
end # ending convexProblem

function nonConvexProblem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                       retToScale,prodStructure,
                       directions)

    nDMUs = size(gI,1)
    (dirGI,dirBI,dirGO,dirBO) = directions

    # TODO: consider bad inputs, consider directions, consider multiplicative/additive structure
    gI_ratio  =  gI₀' ./ gI
    gO_ratio  =  gO ./ gO₀'
    bO_ratio  =  bO ./ bO₀'


    if retToScale == "constant"
        normalFrontierDistances   = [minimum(gI_ratio[z,:]) * minimum(hcat(gO_ratio,bO_ratio)[z,:]) for z in 1:nDMUs]
        bFrontierDistances        = [minimum(gI_ratio[z,:]) * minimum(gO_ratio[z,:]) * (minimum(gO_ratio[z,:]) >= minimum(bO_ratio[z,:])) for z in 1:nDMUs]
    else
        normalFrontierDistances   = [minimum(gI_ratio[z,:]) * minimum(hcat(gO_ratio,bO_ratio)[z,:] * ( minimum(gI_ratio[z,:]) >= (1.0 - eps()) ) ) for z in 1:nDMUs]
        bFrontierDistances        = [minimum(gI_ratio[z,:]) * minimum(gO_ratio[z,:]) * (minimum(gO_ratio[z,:]) >= minimum(bO_ratio[z,:]) &&  minimum(gI_ratio[z,:]) >= (1.0 - eps())    ) for z in 1:nDMUs]
    end
    return min(maximum(normalFrontierDistances),maximum(bFrontierDistances))
end


end # module
