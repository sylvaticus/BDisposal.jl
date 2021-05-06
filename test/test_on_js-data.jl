using Test, DataFrames, CSV, BDisposal
println("Testing BDisposal on js-data...")

# oecd-data test with bad inputs separated...
data = CSV.read(joinpath(dirname(pathof(BDisposal)),"..","test","data","js-data","oecd.txt"),DataFrame; delim=' ',ignorerepeated=true,copycols=true,header=false)
rename!(data,[:ccid,:year,:gdp,:co2,:capital,:labour,:energy])
goodInputsLabels  = ["capital","labour"]
badInputsLabels   = ["energy"]
goodOutputsLabels = ["gdp"]
badOutputsLabels  = ["co2"]
sort!(data, [:year, :ccid]) # sort data by period and dmu
periods = unique(data.year)
dmus    = unique(data.ccid)

nGI, nBI, nGO, nBO, nPer, nDMUs,  = length(goodInputsLabels), length(badInputsLabels), length(goodOutputsLabels), length(badOutputsLabels), length(periods),length(dmus)

gI = Array{Float64}(undef, (nDMUs,nGI,nPer))
bI = Array{Float64}(undef, (nDMUs,nBI,nPer))
gO = Array{Float64}(undef, (nDMUs,nGO,nPer))
bO = Array{Float64}(undef, (nDMUs,nBO,nPer))

myeps = 0.000000001
for (p,period) in enumerate(periods)
    periodData = data[data.year .== period,:]
    gI[:,:,p] = Matrix{Float64}(periodData[:,goodInputsLabels])
    if nBI > 0
        bI[:,:,p] = Matrix{Float64}(periodData[:,badInputsLabels])
    end
    gO[:,:,p] = Matrix{Float64}(periodData[:,goodOutputsLabels])
    bO[:,:,p] = Matrix{Float64}(periodData[:,badOutputsLabels])
end

t = 5
gIₜ = gI[:,:,t]
gIᵤ = gI[:,:,t+1]
bIₜ = bI[:,:,t]
bIᵤ = bI[:,:,t+1]
gOₜ = gO[:,:,t]
gOᵤ = gO[:,:,t+1]
bOₜ = bO[:,:,t]
bOᵤ = bO[:,:,t+1]

z = 1
gI₀ₜ = gI[z,:,t]
gI₀ᵤ = gI[z,:,t+1]
bI₀ₜ = bI[z,:,t]
bI₀ᵤ = bI[z,:,t+1]
gO₀ₜ = gO[z,:,t]
gO₀ᵤ = gO[z,:,t+1]
bO₀ₜ = bO[z,:,t]
bO₀ᵤ = bO[z,:,t+1]

for z in 1:nDMUs

println(z)

gI₀ₜ = gI[z,:,t]
gI₀ᵤ = gI[z,:,t+1]
bI₀ₜ = bI[z,:,t]
bI₀ᵤ = bI[z,:,t+1]
gO₀ₜ = gO[z,:,t]
gO₀ᵤ = gO[z,:,t+1]
bO₀ₜ = bO[z,:,t]
bO₀ᵤ = bO[z,:,t+1]

 # Everything at time t....
gI_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=(-1,0,0,0),startValues=(),forceLinearModel=true)
gI_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=(1,0,0,0),startValues=(),forceLinearModel=true)
@test gI_mult >= 1.0 - myeps
@test begin (gI_mult ≈ 1.0) ?  (0.0 - myeps) < gI_add < (0.0 + myeps) : gI_add > 0.0 - myeps end

bI_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=(0,-1,0,0),startValues=(),forceLinearModel=true)
bI_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=(0,1,0,0),startValues=(),forceLinearModel=true)
@test bI_mult >= 1.0 - myeps
@test begin (bI_mult ≈ 1.0) ?  (0.0 - myeps) < bI_add < (0.0 + myeps) : bI_add > 0.0 - myeps end

gO_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true)
gO_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true)
@test gO_mult >= 1.0 - myeps
@test begin (gO_mult ≈ 1.0) ?  (0.0 - myeps) < gO_add < (0.0 + myeps) : gO_add > 0.0 - myeps end

bO_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,0,-1),startValues=(),forceLinearModel=true)
bO_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,0,1),startValues=(),forceLinearModel=true)
@test bO_mult >= 1.0 - myeps
@test begin (bO_mult ≈ 1.0) ?  (0.0 - myeps) < bO_add < (0.0 + myeps) : bO_add > 0.0 - myeps end


# Observation at time t+1, all other dmu at time t....
gI_t̃_mult = convexProblem(gI₀ᵤ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(-1,0,0,0),startValues=(),forceLinearModel=true, crossTime=true)
gI_t̃_add = convexProblem(gI₀ᵤ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(1,0,0,0),startValues=(),forceLinearModel=true,crossTime=true)
@test  isapprox(gI_t̃_add,1 - 1/gI_t̃_mult, atol=0.000001)

bI_t̃_mult = convexProblem(gI₀ₜ,bI₀ᵤ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,-1,0,0),startValues=(),forceLinearModel=true, crossTime=true)
bI_t̃_add = convexProblem(gI₀ₜ,bI₀ᵤ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(0,1,0,0),startValues=(),forceLinearModel=true, crossTime=true)
@test bI_t̃_add ≈ 1 - 1/bI_t̃_mult

gO_t̃_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ᵤ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true, crossTime=true)
gO_t̃_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ᵤ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true, crossTime=true)
@test gO_t̃_add ≈ (gO_t̃_mult -1)

bO_t̃_mult = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,0,-1),startValues=(),forceLinearModel=true, crossTime=true)
bO_t̃_add = convexProblem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,0,1),startValues=(),forceLinearModel=true, crossTime=true)
@test bO_t̃_add ≈ 1 - 1/bO_t̃_mult


# Observation at time t, all other dmu at time t+1....
gI_ũ_mult = convexProblem(gI₀ₜ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(-1,0,0,0),startValues=(),forceLinearModel=true, crossTime=true)
gI_ũ_add = convexProblem(gI₀ₜ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(1,0,0,0),startValues=(),forceLinearModel=true,crossTime=true)
@test  isapprox(gI_ũ_add,1 - 1/gI_ũ_mult, atol=0.000001)

bI_ũ_mult = convexProblem(gI₀ᵤ,bI₀ₜ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,-1,0,0),startValues=(),forceLinearModel=true, crossTime=true)
bI_ũ_add = convexProblem(gI₀ᵤ,bI₀ₜ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(0,1,0,0),startValues=(),forceLinearModel=true, crossTime=true)
@test bI_ũ_add ≈ 1 - 1/bI_ũ_mult

gO_ũ_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ₜ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true, crossTime=true)
gO_ũ_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ₜ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true, crossTime=true)
@test gO_ũ_add ≈ (gO_ũ_mult -1)

bO_ũ_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,0,-1),startValues=(),forceLinearModel=true, crossTime=true)
bO_ũ_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,0,1),startValues=(),forceLinearModel=true, crossTime=true)
@test bO_ũ_add ≈ 1 - 1/bO_ũ_mult

# Everything at time t+1 (u)....

gI_u_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=(-1,0,0,0),startValues=(),forceLinearModel=true)
gI_u_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=(1,0,0,0),startValues=(),forceLinearModel=true)

bI_u_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=(0,-1,0,0),startValues=(),forceLinearModel=true)
bI_u_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=(0,1,0,0),startValues=(),forceLinearModel=true)

gO_u_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true)
gO_u_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,1,0),startValues=(),forceLinearModel=true)

bO_u_mult = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=(0,0,0,-1),startValues=(),forceLinearModel=true)
bO_u_add = convexProblem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=(0,0,0,1),startValues=(),forceLinearModel=true)


idx_gi_t = gI_mult
idx_bi_t = bI_mult
idx_go_t = gO_mult
idx_bo_t = bO_mult
idx_gi_t̃ = gI_t̃_mult
idx_bi_t̃ = bI_t̃_mult
idx_go_t̃ = gO_t̃_mult
idx_bo_t̃ = bO_t̃_mult
idx_gi_u = gI_u_mult
idx_bi_u = bI_u_mult
idx_go_u = gO_u_mult
idx_bo_u = bO_u_mult
idx_gi_ũ = gI_ũ_mult
idx_bi_ũ = bI_ũ_mult
idx_go_ũ = gO_ũ_mult
idx_bo_ũ = bO_ũ_mult


idx_i_t = (idx_gi_t̃/idx_gi_t) * (idx_bi_t̃/idx_bi_t)
idx_o_t = (idx_go_t/idx_go_t̃) * (idx_bo_t̃/idx_bo_t)
idx_t   = idx_o_t/idx_i_t
idx_i_u = (idx_gi_u/idx_gi_ũ) * (idx_bi_u/idx_bi_ũ)
idx_o_u = (idx_go_ũ/idx_go_u) * (idx_bo_u/idx_bo_ũ)
idx_u   = idx_o_u/idx_i_u
idx     = (idx_t * idx_u)^(1/2)


(idx_gi_t, idx_bi_t, idx_go_t, idx_bo_t) = (gI_add, bI_add, gO_add, bO_add)
(idx_gi_t̃, idx_bi_t̃, idx_go_t̃, idx_bo_t̃) = (gI_t̃_add, bI_t̃_add, gO_t̃_add, bO_t̃_add)
(idx_gi_u, idx_bi_u, idx_go_u, idx_bo_u) = (gI_u_add, bI_u_add, gO_u_add, bO_u_add)
(idx_gi_ũ, idx_bi_ũ, idx_go_ũ, idx_bo_ũ) = (gI_ũ_add, bI_ũ_add, gO_ũ_add, bO_ũ_add)

idx_i_t = (idx_gi_t̃ - idx_gi_t) + (idx_bi_t̃ - idx_bi_t)
idx_o_t = (idx_go_t - idx_go_t̃) + (idx_bo_t̃ - idx_bo_t)
idx_t   = idx_o_t - idx_i_t
idx_i_u = (idx_gi_u - idx_gi_ũ) + (idx_bi_u - idx_bi_ũ)
idx_o_u = (idx_go_ũ - idx_go_u) + (idx_bo_u - idx_bo_ũ)
idx_u   = idx_o_u - idx_i_u
idx     = (idx_t + idx_u) / 2





end






oecdAnalysis  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=true)




oecdAnalysis.prodIndexes

oecdAnalysisA  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="addittive",convexAssumption=true)













oecdAnalysisA.prodIndexes
