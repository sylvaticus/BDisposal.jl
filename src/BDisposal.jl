module BDisposal

using DataFrames, CSV, JuMP, GLPK

export efficiencyScores


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
function efficiencyScores(inputs,goodOutputs,badOutputs,data;
                                   retToScale="constant",effAnalysisType="static",formattedOutput=false)
    # Data processing
    sort!(data, [:period, :dmu]) # sort data by period and dmu
    periods = unique(data.period)
    dmus    = unique(data.dmu)
    nI, nGO, nBO, nPer, nDMUs = length(inputs), length(goodOutputs), length(badOutputs), length(periods),length(dmus)

    λs_convex    = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)
    λs_nonconvex = Array{Union{Missing,Float64,String}}(undef,nDMUs,nPer) # Output matrix hosting the results (nDMUs,nPer)

    # Loping over the periods
    for (pIdx,period) in enumerate(periods)
        #period = periods[1]
        #pIdx = 1
        periodData = data[data.period .== period,:]
        inp = convert(Matrix{Float64},periodData[:,inputs])
        gO  = convert(Matrix{Float64},periodData[:,goodOutputs])
        bO  = convert(Matrix{Float64},periodData[:,badOutputs])

        # Solving for each dmu
        for dmuIdx in 1:nDMUs
            #dmuIdx = 3
            inp₀ = inp[dmuIdx,:] # inputs of this dmu (x in the model)
            gO₀  = gO[dmuIdx,:] # good outputs of this dmu (y in the model)
            bO₀  = bO[dmuIdx,:] # bad outputs of this dmu (y in the model)

            # CONVEX analysis....

            # Model declaration (efficiency model)
            effmodel = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # we choose GLPK with a verbose output

            # Defining variables
            @variables effmodel begin
                θ[z in 1:nDMUs] >= 0 # thetas (nDMUs)
                μ[z in 1:nDMUs] >= 0 # mus (nDMUs)
                λ
            end

            # Defining constraints
            if retToScale == "constant"
                @constraints effmodel begin
                    cInp1[i in 1:nI],   # first constraint on inputs
                       inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                    cInp2[i in 1:nI],   # second constraint on inputs
                       inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                    cGO1[j in 1:nGO],   # first constraint on good outputs
                       λ * gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                    cBO1[j in 1:nBO],   # first constraint on bad outputs
                      λ * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                    cGO2[j in 1:nGO],   # second constraint on good outputs
                      λ * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                    cBO2[j in 1:nBO],   # second constraint on bad outputs
                      λ * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                end
            elseif retToScale == "variable"
            @constraints effmodel begin
                cInp1[i in 1:nI],   # first constraint on inputs
                   inp₀[i] >= sum(θ[z]*inp[z,i] for z in 1:nDMUs)
                cInp2[i in 1:nI],   # second constraint on inputs
                   inp₀[i] >= sum(μ[z]*inp[z,i] for z in 1:nDMUs)
                cGO1[j in 1:nGO],   # first constraint on good outputs
                   λ * gO₀[j] <= sum(θ[z]*gO[z,j] for z in 1:nDMUs)
                cBO1[j in 1:nBO],   # first constraint on bad outputs
                  λ * bO₀[j] >= sum(θ[z]*bO[z,j] for z in 1:nDMUs)
                cGO2[j in 1:nGO],   # second constraint on good outputs
                  λ * gO₀[j] <= sum(μ[z]*gO[z,j] for z in 1:nDMUs)
                cBO2[j in 1:nBO],   # second constraint on bad outputs
                  λ * bO₀[j] <= sum(μ[z]*bO[z,j] for z in 1:nDMUs)
                sumtheta,
                  sum(θ[z] for z in 1:1:nDMUs) == 1
                summu,
                  sum(μ[z] for z in 1:1:nDMUs) == 1
            end
            else
                @error "Unknow returns to scale type."
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
            #println("*******************************************************************")
            #println("** Solutions for the $(dmus[dmuIdx]) decision maker unit for period $(period) **")
            if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(effmodel)
                #if (status == MOI.OPTIMAL)
                #    println("** Problem solved correctly **")
                #else
                #    println("** Problem returned a (possibly suboptimal) solution **")
                #end
                #println("- Objective value: ", objective_value(effmodel))
                #println("- Optimal thetas:")
                #optθ = value.(θ)
                #[println("$(dmus[z]): $(optθ[z])") for z in 1:nDMUs]
                #println("- Optimal mus:")
                #optμ = value.(μ)
                #[println("$(dmus[z]): $(optμ[z])") for z in 1:nDMUs]
                #optλ = value(λ)
                #println("λ: $(optλ)")
                λs_convex[dmuIdx,pIdx] = value(λ)
                #println("- Dual of the first contraint on bad output:")
                #[println("$(badOutputs[j]) = $(shadow_price(cBO1[j]))") for j in 1:nBO]
            else
                #println("The model was not solved correctly.")
                #println(status)
                λs_convex[dmuIdx,pIdx] = missing
            end

            # NON Convex analysis
            gO_ratio  =  gO ./ gO₀'
            bO_ratio  =  bO ./ bO₀'
            inp_ratio = inp₀' ./ inp

            if retToScale == "constant"
                normalFrontierDistances   = [minimum(inp_ratio[z,:]) * minimum(hcat(gO_ratio,bO_ratio)[z,:]) for z in 1:nDMUs]
                bFrontierDistances        = [minimum(inp_ratio[z,:]) * minimum(gO_ratio[z,:]) * (minimum(gO_ratio[z,:]) >= minimum(bO_ratio[z,:])) for z in 1:nDMUs]
            else
                normalFrontierDistances   = [minimum(inp_ratio[z,:]) * minimum(hcat(gO_ratio,bO_ratio)[z,:] * ( minimum(inp_ratio[z,:]) >= (1.0 - eps()) ) ) for z in 1:nDMUs]
                bFrontierDistances        = [minimum(inp_ratio[z,:]) * minimum(gO_ratio[z,:]) * (minimum(gO_ratio[z,:]) >= minimum(bO_ratio[z,:]) &&  minimum(inp_ratio[z,:]) >= (1.0 - eps())    ) for z in 1:nDMUs]
            end
            λs_nonconvex[dmuIdx,pIdx] = min(maximum(normalFrontierDistances),maximum(bFrontierDistances))

        end # end of each dmu
    end # endof each period

    nonConvTest_value = λs_nonconvex ./ λs_convex
    nonConvTest       = nonConvTest_value .< (1.0 - eps())
    λs = min.(λs_nonconvex,λs_convex)

    if formattedOutput
        # Add periods as headers and decision making names as first column
        λs = vcat(periods',λs)
        λs = hcat(vcat("",dmus),λs)
        λs_convex = vcat(periods',λs_convex)
        λs_convex = hcat(vcat("",dmus),λs_convex)
        λs_nonconvex = vcat(periods',λs_nonconvex)
        λs_nonconvex = hcat(vcat("",dmus),λs_nonconvex)
        nonConvTest_value = vcat(periods',nonConvTest_value)
        nonConvTest_value = hcat(vcat("",dmus),nonConvTest_value)
        nonConvTest = vcat(periods',nonConvTest)
        nonConvTest = hcat(vcat("",dmus),nonConvTest)
    end
    return λs, λs_convex, λs_nonconvex, nonConvTest_value, nonConvTest
end

end # module
