###############################################################
# Part of BDisposal.jl package
###############################################################

using LinearAlgebra

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
           #effmodel = Model(() -> AmplNLWriter.Optimizer(path, ["print_level=0"]))
           #effmodel = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
           effmodel = Model(() -> AmplNLWriter.Optimizer("ipopt",["print_level=0"]))
           #effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF))

           # Defining variables
           @variables effmodel begin
               0 <= θ[z in 1:nDMUs] <= Inf, (start = startθ) # thetas (nDMUs)
               0 <= μ[z in 1:nDMUs] <= Inf, (start = startμ) # mus (nDMUs)
           end
           if !crossTime
              @variable(effmodel, 1 <= λ <= Inf, start = startλ)
           else
              @variable(effmodel, 0 <= λ <= Inf, start = startλ)
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
   else # prodStructure is then additive...
       # Model declaration (efficiency model)
       effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
       # Defining variables
       @variables effmodel begin
           θ[z in 1:nDMUs] >= 0 # thetas (nDMUs)
           μ[z in 1:nDMUs] >= 0 # mus (nDMUs)
       end
       if !crossTime
          @variable(effmodel, 0 <= λ <= Inf)
       else
          @variable(effmodel, -Inf <= λ <= Inf)
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


"""
     dmuEfficiency(I₀,O₀,I,O)

Compute the efficiency score for a DMU using vanilla Data Envelope Analysis.

## Parameters:
  *  `I₀`: This DMU inputs (nI)
  *  `O₀`: This DMU outputs (n0)
  *  `I`: All DMUs inputs (nDMU x nI)
  *  `O`: All DMUs outputs (nDMU x n0)

## Returns:
- A scalar representing the efficiency score of the DMU
"""
function dmuEfficiency(I₀,O₀,I,O)
    (nDMU,nI,nO) = size(I,1), size(I,2), size(O,2)

    effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
    # Defining variables
    @variables effmodel begin
        wI[i in 1:nI] >= 0  # inputs weigths
        wO[o in 1:nO] >= 0  # outputs weigths
    end

    # Defining constraints
    @constraints effmodel begin
        effConstr[d in 1:nDMU],   #
           dot(O[d,:], wO) - dot(I[d,:], wI) <=  0.0
        iRegularisationConstr, # regularisation of inputs
            dot(I₀,wI) == 1.0

    end

    # Defining objective
    @objective effmodel Max begin
         dot(O₀, wO)
    end

    # Printing the model (for debugging)
    #print(effmodel) # The model in mathematical terms is printed

    # Solving the model
    optimize!(effmodel)

    status = termination_status(effmodel)

    #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
        wI = value.(wI)
        wO = value.(wO)
        duals = dual.(effConstr)
        obj = dot(O₀,wO) / dot(I₀,wI)
        refSet = Dict{Int64,Float64}()
        [refSet[i] = -f for (i,f) in enumerate(duals) if f != 0.0]
        eff  = (isapprox.(obj,1.0,atol=1e-15) && !any(isapprox.(wI,0.0,atol=1e-15)) && !any(isapprox.(wI,0.0,atol=1e-15))) # efficient if the objective is 1 AND no weigths are zero (otherwise it is on one periphal side of the frontier and shadowed by other points)
        return (eff, obj=obj, wI=wI, wO = wO, refSet=refSet)
    else
        return missing
    end
end

#cd(@__DIR__)
#using Pkg
#Pkg.activate("..")
