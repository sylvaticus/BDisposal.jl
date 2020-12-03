module BDisposal

using JuMP, Ipopt, GLPK, AmplNLWriter

export efficiencyScores, prodIndex

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("IndividualDMUProblem.jl")


end # module
