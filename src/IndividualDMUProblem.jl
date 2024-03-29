###############################################################
# Part of BDisposal.jl package
###############################################################

using LinearAlgebra

# Dispatching to either convexProblem or nonConvexProblem
function problem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                 retToScale,prodStructure,convexAssumption=true,
                 directions=(),
                 startValues=(),crossTime=false,forceLinearModel=false)

     if convexAssumption
         return convexProblem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                          retToScale=retToScale,prodStructure=prodStructure,
                          directions=directions,
                          startValues=startValues,crossTime=crossTime,forceLinearModel=forceLinearModel)
     else
         return nonConvexProblem(gI₀,bI₀,gO₀,bO₀,gI,bI,gO,bO;
                          retToScale=retToScale,prodStructure=prodStructure,
                          directions=directions)
     end

end

function convexProblem(inp₀,bInp₀,gO₀,bO₀,inp,bInp,gO,bO;
                       retToScale,prodStructure,
                       directions,
                       startValues,crossTime=false,forceLinearModel=false)

   nI, nGO, nBO, nDMUs, nbI = length(inp₀), length(gO₀), length(bO₀), size(inp,1), length(bInp₀)
   (dirGI,dirBI,dirGO,dirBO) = directions
   if !forceLinearModel
       (startθ,startμ,startλ)    = startValues
   end

   if prodStructure == "multiplicative"
       if !forceLinearModel
               # Redirecting the output to /deV/null as no way to suppress ipopt output
               oldstd = stdout
               redirect_stdout(open("nul", "w"))

               # Model declaration (efficiency model)
               #effmodel = Model(() -> AmplNLWriter.Optimizer(Ipopt.amplexe))
               effmodel = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
               #MOI.set(effmodel, MOI.RawOptimizerAttribute("print_level"), 0)
               #set_optimizer_attribute(effmodel, "print_level", 0)

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
               JuMP.optimize!(effmodel)
               redirect_stdout(oldstd) # recover original stdout

               # Copy the results to the output matrix...
               status = termination_status(effmodel)

               #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
               if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
                   return value(λ)
               else
                   return missing
               end
       else #forceLinearModel is then true

               # Model declaration (efficiency model)
               #effmodel = Model(() -> AmplNLWriter.Optimizer(path, ["print_level=0"]))
               #effmodel = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
               effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
                #effmodel = Model(optimizer_with_attributes(GLPK.Optimizer))
               #effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF))

               # Defining variables
               @variables effmodel begin
                   0 <= θ[z in 1:nDMUs] <= Inf# thetas (nDMUs)
                   0 <= μ[z in 1:nDMUs] <= Inf # mus (nDMUs)
               end

               if dirGO == 1
                   if !crossTime # we  keep the original λ
                  @variable(effmodel, 1 <= λ <= Inf)
                   else
                      @variable(effmodel, 0 <= λ <= Inf)
                   end
                   # @variable(effmodel, λ)


               else
                   if !crossTime # we use as variable its reverse 1/λ
                      @variable(effmodel, 0 <= invλ <= 1)
                   else
                      @variable(effmodel, 0 <= invλ <= Inf)
                   end
                   # @variable(effmodel, invλ)
               end

               if directions == (-1,0,0,0)
                   # Defining constraints
                   @constraints effmodel begin
                       cgInp1[i in 1:nI],   # constraint on good inputs
                          invλ * inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                       cbInp1[i in 1:nbI], # constraint on bad inputs
                           bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
                       cgInp2[i in 1:nI],   # second constraint on good inputs
                          invλ * inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                       cbInp2[i in 1:nbI], # second constraint on bad inputs
                           bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
                       cgO1[j in 1:nGO],   # first constraint on good outputs
                           gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                       cbO1[j in 1:nBO],   # first constraint on bad outputs
                           bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                       cgO2[j in 1:nGO],   # second constraint on good outputs
                           gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                       cbO2[j in 1:nBO],   # second constraint on bad outputs
                           bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                   end
              elseif directions == (0,-1,0,0)
                  # Defining constraints
                  @constraints effmodel begin
                      cgInp1[i in 1:nI],   # constraint on good inputs
                          inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp1[i in 1:nbI], # constraint on bad inputs
                         invλ * bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgInp2[i in 1:nI],   # second constraint on good inputs
                          inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp2[i in 1:nbI], # second constraint on bad inputs
                         invλ * bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgO1[j in 1:nGO],   # first constraint on good outputs
                          gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO1[j in 1:nBO],   # first constraint on bad outputs
                          bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                      cgO2[j in 1:nGO],   # second constraint on good outputs
                          gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO2[j in 1:nBO],   # second constraint on bad outputs
                          bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                  end
              elseif directions == (0,0,1,0)
                  @constraints effmodel begin
                      cgInp1[i in 1:nI],   # constraint on good inputs
                          inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp1[i in 1:nbI], # constraint on bad inputs
                          bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgInp2[i in 1:nI],   # second constraint on good inputs
                          inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp2[i in 1:nbI], # second constraint on bad inputs
                          bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgO1[j in 1:nGO],   # first constraint on good outputs
                         λ * gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO1[j in 1:nBO],   # first constraint on bad outputs
                          bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                      cgO2[j in 1:nGO],   # second constraint on good outputs
                         λ * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO2[j in 1:nBO],   # second constraint on bad outputs
                         bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                  end
              elseif directions == (0,0,0,-1)
                  # Defining constraints
                  @constraints effmodel begin
                      cgInp1[i in 1:nI],   # constraint on good inputs
                         inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp1[i in 1:nbI], # constraint on bad inputs
                         bInp₀[i] <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgInp2[i in 1:nI],   # second constraint on good inputs
                         inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                      cbInp2[i in 1:nbI], # second constraint on bad inputs
                          bInp₀[i] >= sum(μ[z]*bInp[z,i] for z in 1:nDMUs)
                      cgO1[j in 1:nGO],   # first constraint on good outputs
                          gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO1[j in 1:nBO],   # first constraint on bad outputs
                         invλ * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                      cgO2[j in 1:nGO],   # second constraint on good outputs
                         gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                      cbO2[j in 1:nBO],   # second constraint on bad outputs
                         invλ * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                  end

              else
                  @error "Directions not suppported when the option forceLinearModel is true"
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
               if dirGO == 1
                   @objective effmodel Max begin
                       λ
                   end
               else
                   @objective effmodel Min begin
                       invλ
                   end
               end
               # Printing the model (for debugging)
               #print(effmodel) # The model in mathematical terms is printed


               JuMP.optimize!(effmodel)

               # Copy the results to the output matrix...
               status = termination_status(effmodel)

               #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
               if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
                   if dirGO == 1
                       return value(λ)
                   else
                       return 1/value(invλ)
                   end
               else
                   return missing
               end


        end   # end of forCeLinearModel == true
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
              bO₀[j] - λ * dirBO * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
           cgO2[j in 1:nGO],   # second constraint on good outputs
              gO₀[j] + λ*dirGO * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
           cbO2[j in 1:nBO],   # second constraint on bad outputs
              bO₀[j] - λ*dirBO * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
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
       # print(effmodel) # The model in mathematical terms is printed

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

    if retToScale == "constant" && directions != (0,0,1,-1) @error "Not supported" end
    nGI, nGO, nBO, nDMUs, nBI = length(gI₀), length(gO₀), length(bO₀), size(gI,1), length(bI₀)
    (dirGI,dirBI,dirGO,dirBO) = directions


    if directions == (-1,0,0,0) && prodStructure == "multiplicative"
        gI_ratio  =  gI ./ gI₀' # nDMU x ngI
        # Normal dist function
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : minimum(maximum(gI_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : minimum(maximum(gI_ratio[globalContraint,:],dims=2))
        effscore = max(effScore_normal,effScore_bfrontier)
        return effscore
    elseif  directions == (1,0,0,0) && prodStructure == "additive"
        gI_ratio  =  1 .- (gI ./ gI₀') # nDMU x ngI
        # Normal dist function
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : maximum(minimum(gI_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : maximum(minimum(gI_ratio[globalContraint,:],dims=2))
        effscore = min(effScore_normal,effScore_bfrontier)
        return effscore

    elseif  directions == (0,-1,0,0) && prodStructure == "multiplicative"
        bI_ratio  =  bI ./ bI₀' # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : minimum(maximum(bI_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : minimum(maximum(bI_ratio[globalContraint,:],dims=2))
        effscore = max(effScore_normal,effScore_bfrontier)
        return effscore
    elseif  directions == (0,1,0,0) && prodStructure == "additive"
        bI_ratio  =  1 .- (bI ./ bI₀') # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : maximum(minimum(bI_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : maximum(minimum(bI_ratio[globalContraint,:],dims=2))
        effscore = min(effScore_normal,effScore_bfrontier)
        return effscore
    elseif  directions == (0,0,1,0) && prodStructure == "multiplicative"
        gO_ratio  =  gO₀' ./ gO  # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : minimum(maximum(gO_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .<=   bI, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : minimum(maximum(gO_ratio[globalContraint,:],dims=2))
        effscore = max(effScore_normal,effScore_bfrontier)
        return 1/effscore
    elseif  directions == (0,0,1,0) && prodStructure == "additive"
        gO_ratio  =  (gO₀' ./ gO ) .- 1 # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        bOConstraint = all(bO₀'  .<=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,bOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : maximum(minimum(gO_ratio[globalContraint,:],dims=2)) #todo check this
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .<=   bI, dims=2)
        bOConstraint = all(bO₀'  .>=   bO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,bOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : maximum(minimum(gO_ratio[globalContraint,:],dims=2)) # todo check this
        effscore = min(effScore_normal,effScore_bfrontier)
        return -effscore
    elseif  directions == (0,0,0,-1) && prodStructure == "multiplicative"
        bO_ratio  =  bO ./ bO₀' # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,gOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : minimum(minimum(bO_ratio[globalContraint,:],dims=2))
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .<=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,gOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : minimum(maximum(bO_ratio[globalContraint,:],dims=2))
        effscore = max(effScore_normal,effScore_bfrontier)
        return 1/effscore
    elseif  directions == (0,0,0,1) && prodStructure == "additive"
        bO_ratio  =  1 .- (bO ./ bO₀') # nDMU x ngI
        # Normal dist function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .>=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,gOConstraint), dims=2), dims=2)
        effScore_normal = ! any(globalContraint) ? missing : maximum(minimum(bO_ratio[globalContraint,:],dims=2)) #todo check this
        # Bdisposal Distance function
        gIConstraint = all(gI₀'  .>=   gI, dims=2)
        bIConstraint = all(bI₀'  .<=   bI, dims=2)
        gOConstraint = all(gO₀'  .<=   gO, dims=2)
        globalContraint = dropdims(all(hcat(gIConstraint,bIConstraint,gOConstraint), dims=2), dims=2)
        effScore_bfrontier = ! any(globalContraint) ? missing : minimum(minimum(bO_ratio[globalContraint,:];dims=2)) # todo check this
        effscore = min(effScore_normal,effScore_bfrontier)
        return effscore
    elseif  directions == (0,0,1,-1) && prodStructure == "multiplicative"
        gI_ratio  =  gI₀' ./ gI
        bI_ratio  =  bI₀' ./ bI
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

    else
        @error "Directions not supported"
        return missing
    end

end # end function


"""
     dmuEfficiency(I₀,O₀,I,O)

Compute the efficiency score for a DMU using vanilla Data Envelope Analysis.

## Parameters:
  *  `I₀`: This DMU inputs (nI)
  *  `O₀`: This DMU outputs (n0)
  *  `I`: All DMUs inputs (nDMU x nI)
  *  `O`: All DMUs outputs (nDMU x n0)

## Returns:
- A named tuple with:
  - `computed`: a boolean to indicase wheter the underlying optimisation succeeded
  - `eff`: a boolean to indicate if the DMU is efficient or not
  - `obj`: a scalar representing the efficiency score of the DMU
  - `wI`: a vector of the optimal input weights
  - `wO`: a vector of the optimal output weights
  - `refSet`: a dictionary having as key the numeral of the reference DMUs and
  values the relative constraint duals
"""
function dmuEfficiency(I₀,O₀,I,O) # # Cooper & oth., p.23-24 and 43
    (nDMU,nI,nO) = size(I,1), size(I,2), size(O,2)

    effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a null output
    # Defining variables
    @variables effmodel begin
        wI[i in 1:nI] >= 0  # inputs weigths
        wO[o in 1:nO] >= 0  # outputs weigths
    end

    # Defining constraints
    @constraints effmodel begin
        othDMUsEffConstr[d in 1:nDMU],   #
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

    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
        wI = value.(wI)
        wO = value.(wO)
        othDMUsEffConstrDuals = -dual.(othDMUsEffConstr)
        iRegConstrDual = -dual(iRegularisationConstr)
        obj = dot(O₀,wO) / dot(I₀,wI)
        refSet = Dict{Int64,Float64}()
        [refSet[i] = f for (i,f) in enumerate(othDMUsEffConstrDuals) if f != 0.0]
        eff  = (isapprox.(obj,1.0,atol=1e-15) && !any(isapprox.(wI,0.0,atol=1e-15)) && !any(isapprox.(wI,0.0,atol=1e-15))) # efficient if the objective is 1 AND no weigths are zero (otherwise it is on one periphal side of the frontier and shadowed by other points)
        return (computed=true,eff=eff, obj=obj, wI=wI, wO = wO, refSet=refSet, othDMUsEffConstrDuals=othDMUsEffConstrDuals,iRegConstrDual=iRegConstrDual)
    else
        return (computed=false,eff=missing,obj=missing,wI=missing,wO = missing,refSet=missing,othDMUsEffConstrDuals=missing,iRegConstrDual=missing)
    end
end

# Cooper & oth., p.43
"""
     dmuEfficiencyDual(I₀,O₀,I,O)

Compute the efficiency score for a DMU using vanilla Data Envelope Analysis in the dual representation using the two pass method.

## Parameters:
  *  `I₀`: This DMU inputs (nI)
  *  `O₀`: This DMU outputs (n0)
  *  `I`: All DMUs inputs (nDMU x nI)
  *  `O`: All DMUs outputs (nDMU x n0)

## Returns:
- A named tuple with:
  - `computed`: a boolean to indicase wheter the underlying optimisation succeeded
  - `eff`: a boolean to indicate if the DMU is efficient or not
  - `obj`: a scalar representing the efficiency score of the DMU
  - `λ`
  - `θ`
  - `dualsI`
  - `dualsO`
  - `s⁻`
  - `s⁺`
"""
function dmuEfficiencyDual(I₀,O₀,I,O)
    (nDMU,nI,nO) = size(I,1), size(I,2), size(O,2)

    effmodeldual = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
    # Defining variables
    @variables effmodeldual begin
        λ[i in 1:nDMU] >= 0  # the constrain dual in the primal
        θ                    # the costs per output ???
    end

    # Defining constraints
    @constraints effmodeldual begin
        inConstr[i in 1:nI],   #
           θ*I₀[i] - dot(I[:,i], λ) >=  0.0
        oRegularisationConstr[j in 1:nO], # regularisation of outputs ??
            dot(O[:,j],λ) >= O₀[j]
    end

    # Defining objective
    @objective effmodeldual Min begin
         θ
    end

    # Printing the model (for debugging)
    #print(effmodel) # The model in mathematical terms is printed

    # Solving the model
    optimize!(effmodeldual)

    status = termination_status(effmodeldual)

    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodeldual)
        λv = value.(λ)
        θv = value(θ)
        dualsI = dual.(inConstr)
        dualsO = dual.(oRegularisationConstr)
        obj = objective_value(effmodeldual)
        outPass2 = dmuPass2(θv,I₀,O₀,I,O)
        s⁻ = outPass2.s⁻ # input excesses
        s⁺ = outPass2.s⁺ # output shortfalls
        eff  = (isapprox.(obj,1.0,atol=1e-15) && all(isapprox.(s⁻,0.0,atol=1e-15)) && all(isapprox.(s⁺,0.0,atol=1e-15))) # efficient if the objective is 1 AND no weigths are zero (otherwise it is on one periphal side of the frontier and shadowed by other points)
        #refSet = Dict{Int64,Float64}()
        #[refSet[i] = -f for (i,f) in enumerate(duals) if f != 0.0]
        #eff  = (isapprox.(obj,1.0,atol=1e-15) && !any(isapprox.(wI,0.0,atol=1e-15)) && !any(isapprox.(wI,0.0,atol=1e-15))) # efficient if the objective is 1 AND no weigths are zero (otherwise it is on one periphal side of the frontier and shadowed by other points)
        #return (eff=eff, obj=obj, wI=wI, wO = wO, refSet=refSet)
        return (computed=true, eff=eff, obj=obj,  λ=λv, θ=θv, dualsI=dualsI, dualsO=dualsO, s⁻=s⁻, s⁺=s⁺)
    else
        #return (eff=missing,obj=missing,wI=missing,wO = missing,refSet=missing)
        return (computed=false, eff=eff, obj=missing, λ=missing, θ=missing, dualsI=missing , dualsO=missing , s⁻=missing, s⁺=missing)
    end
end

# Cooper & oth., p.44
function dmuPass2(θ,I₀,O₀,I,O)
    (nDMU,nI,nO) = size(I,1), size(I,2), size(O,2)

    effmodeldual = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
    # Defining variables
    @variables effmodeldual begin
        λ[i in 1:nDMU] >= 0  # the constrain dual in the primal
        s⁻[i in 1:nI]  >= 0  # input excesses
        s⁺[i in 1:nO]  >= 0  # output shortfalls
    end

    # Defining constraints
    @constraints effmodeldual begin
        inConstr[i in 1:nI],   #
           θ*I₀[i] - dot(I[:,i], λ) ==  s⁻[i]
        oRegularisationConstr[j in 1:nO], # regularisation of outputs ??
            dot(O[:,j],λ) - O₀[j] == s⁺[j]
    end

    # Defining objective
    @objective effmodeldual Max begin
         sum(s⁻)+sum(s⁺)
    end

    # Printing the model (for debugging)
    #print(effmodel) # The model in mathematical terms is printed

    # Solving the model
    optimize!(effmodeldual)

    status = termination_status(effmodeldual)

    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodeldual)
        λv = value.(λ)
        s⁻v = value.(s⁻)
        s⁺v = value.(s⁺)
        obj = objective_value(effmodeldual)
        return (computed=true, obj=obj, λ=λv, s⁻=s⁻v, s⁺=s⁺v )
    else
        #return (eff=missing,obj=missing,wI=missing,wO = missing,refSet=missing)
        return (computed=false, obj=missing, λ=missing,s⁻=missing, s⁺=missing )
    end
end

function problemFB(inp₀,bInp₀,gO₀,bO₀,inp,bInp,gO,bO; directions)

    nI, nGO, nBO, nDMUs, nbI = length(inp₀), length(gO₀), length(bO₀), size(inp,1), length(bInp₀)
    (dirGI,dirBI,dirGO,dirBO) = directions


    # Model declaration (efficiency model)
    effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output
    # Defining variables
    @variables effmodel begin
        θ[z in 1:nDMUs] >= 0 # thetas (nDMUs)
        μ[z in 1:nDMUs] >= 0 # mus (nDMUs)
    end

    @variable(effmodel, 0 <= λ <= Inf)

  
    # Defining constraints
    # additive structure
    @constraints effmodel begin
        cgInp2[i in 1:nI],   # second constraint on good inputs
        inp₀[i] - λ*dirGI  >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
        cbInp1[i in 1:nbI], # constraint on bad inputs
        bInp₀[i]- λ*dirBI  <= sum(θ[z]*bInp[z,i] for z in 1:nDMUs)
        cgO2[j in 1:nGO],   # second constraint on good outputs
        gO₀[j] + λ*dirGO  <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
        cbO1[j in 1:nBO],   # first constraint on bad outputs
        bO₀[j] - λ * dirBO  >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)

    end
  
    @constraints effmodel begin
        sumtheta,
        sum(θ[z] for z in 1:1:nDMUs) == 1
        summu,
        sum(μ[z] for z in 1:1:nDMUs) == 1
    end
      
    # Defining objective
    @objective effmodel Max begin
        λ
    end
    # Printing the model (for debugging)
    # print(effmodel) # The model in mathematical terms is printed

    # Solving the model
    optimize!(effmodel)

    # Copy the results to the output matrix...
    status = termination_status(effmodel)
    #println(effmodel)
    #println(status)
    

    #if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED) && has_values(effmodel)
        return value(λ)
    else
        return missing
    end
    
end # ending problemFB
