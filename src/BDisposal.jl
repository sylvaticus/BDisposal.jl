module BDisposal

using JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex, dmuEfficiency, dmuEfficiencyDual
# for debug only:
# export problem, convexProblem, nonConvexProblem, dmuPass2

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("IndividualDMUProblem.jl")


end # module
