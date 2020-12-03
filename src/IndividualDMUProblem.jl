###############################################################
# Part of BDisposal.jl package
###############################################################


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