module BDisposal

using JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex, dmuEfficiency
# for debug only:
export problem, convexProblem

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("IndividualDMUProblem.jl")


end # module
