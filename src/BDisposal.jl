module BDisposal

using Statistics, JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex, prodIndexPF, dmuEfficiency, dmuEfficiencyDual
# for debug only:
# export problem, convexProblem, nonConvexProblem, dmuPass2

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("ProdIndexPF.jl")
include("IndividualDMUProblem.jl")


end # module
