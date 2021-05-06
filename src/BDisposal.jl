module BDisposal

using JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex, dmuEfficiency, dmuEfficiencyDual, dmuPass2
# for debug only:
export problem, convexProblem, nonConvexProblem

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("IndividualDMUProblem.jl")


end # module
