module BDisposal

using Statistics, JuMP, Ipopt, GLPK, AmplNLWriter, PrecompileTools

export efficiencyScores, prodIndex, prodIndexFB, dmuEfficiency, dmuEfficiencyDual
# for debug only:
# export problem, convexProblem, nonConvexProblem, dmuPass2

include("EfficiencyScores.jl")
include("ProdIndex.jl")
include("ProdIndexFB.jl")
include("IndividualDMUProblem.jl")
include("Precompilation.jl")


end # module
