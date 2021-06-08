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
        bI[:,:,p] =
         Matrix{Float64}(periodData[:,badInputsLabels])
    end
    gO[:,:,p] = Matrix{Float64}(periodData[:,goodOutputsLabels])
    bO[:,:,p] = Matrix{Float64}(periodData[:,badOutputsLabels])
end

#=
using Distributions
ud = Uniform(1,100000)

gI = rand(ud,nDMUs,nGI,nPer)
bI = rand(ud,nDMUs,nBI,nPer)
gO = rand(ud,nDMUs,nGO,nPer)
bO = rand(ud,nDMUs,nBO,nPer)
=#


t = length(periods)-1
gIₜ = gI[:,:,t]
gIᵤ = gI[:,:,t+1]
bIₜ = bI[:,:,t]
bIᵤ = bI[:,:,t+1]
gOₜ = gO[:,:,t]
gOᵤ = gO[:,:,t+1]
bOₜ = bO[:,:,t]
bOᵤ = bO[:,:,t+1]

z = 4
gI₀ₜ = gI[z,:,t]
gI₀ᵤ = gI[z,:,t+1]
bI₀ₜ = bI[z,:,t]
bI₀ᵤ = bI[z,:,t+1]
gO₀ₜ = gO[z,:,t]
gO₀ᵤ = gO[z,:,t+1]
bO₀ₜ = bO[z,:,t]
bO₀ᵤ = bO[z,:,t+1]

convex = true

if convex == true
    (dirGIm,dirGIa) = (-1,0,0,0) , (1,0,0,0)
    (dirBIm,dirBIa) = (0,-1,0,0) , (0,1,0,0)
    (dirGOm,dirGOa) = (0,0,1,0)  , (0,0,1,0)
    (dirBOm,dirBOa) = (0,0,0,-1) , (0,0,0,1)
else
    (dirGIm,dirGIa) = (1,0,0,0)  , (1,0,0,0)
    (dirBIm,dirBIa) = (0,1,0,0)  , (0,1,0,0)
    (dirGOm,dirGOa) = (0,0,-1,0) , (0,0,1,0)
    (dirBOm,dirBOa) = (0,0,0,1)  , (0,0,0,-1)
end


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
gI_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=dirGIm,startValues=(),forceLinearModel=true,convexAssumption=convex)
gI_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=dirGIa,startValues=(),forceLinearModel=true,convexAssumption=convex)
@test gI_mult >= 1.0 - myeps
@test begin (gI_mult ≈ 1.0) ?  (0.0 - myeps) < gI_add < (0.0 + myeps) : gI_add > 0.0 - myeps end

bI_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=dirBIm,startValues=(),forceLinearModel=true,convexAssumption=convex)
bI_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=dirBIa,startValues=(),forceLinearModel=true,convexAssumption=convex)
@test bI_mult >= 1.0 - myeps
@test begin (bI_mult ≈ 1.0) ?  (0.0 - myeps) < bI_add < (0.0 + myeps) : bI_add > 0.0 - myeps end

gO_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true,convexAssumption=convex)
gO_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true,convexAssumption=convex)
@test gO_mult >= 1.0 - myeps
@test begin (gO_mult ≈ 1.0) ?  (0.0 - myeps) < gO_add < (0.0 + myeps) : gO_add > 0.0 - myeps end

bO_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true,convexAssumption=convex)
bO_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true,convexAssumption=convex)
@test bO_mult >= 1.0 - myeps
@test begin (bO_mult ≈ 1.0) ?  (0.0 - myeps) < bO_add < (0.0 + myeps) : bO_add > 0.0 - myeps end

# Observation at time t+1, all other dmu at time t....
gI_t̃_mult = problem(gI₀ᵤ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGIm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
gI_t̃_add = problem(gI₀ᵤ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirGIa,startValues=(),forceLinearModel=true,crossTime=true,convexAssumption=convex)
@test  isapprox(gI_t̃_add,1 - 1/gI_t̃_mult, atol=0.000001)
#=
gI_ratio  =  gIₜ ./ (gI₀ᵤ)'
bIConstraint = all((bI₀ₜ)'  .>=   bIₜ, dims=2)
gOConstraint = all((gO₀ₜ)'  .<=   gOₜ, dims=2)
bOConstraint = all((bO₀ₜ)'  .<=   bOₜ, dims=2)
globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
effScore_normal = minimum(maximum(gI_ratio[globalContraint,:],dims=2))
# Bdisposal Distance function
bIConstraint = all(bI₀'  .>=   bI, dims=2)
gOConstraint = all(gO₀'  .<=   gO, dims=2)
bOConstraint = all(bO₀'  .>=   bO, dims=2)
globalContraint = dropdims(all(hcat(bIConstraint,gOConstraint,bOConstraint), dims=2), dims=2)
effScore_bfrontier = minimum(gI_ratio[globalContraint,:])
effscore = max(effScore_normal,effScore_bfrontier)
return effscore

(bO₀ₜ)'  .<=   bOₜ
gIₜ
=#

bI_t̃_mult = problem(gI₀ₜ,bI₀ᵤ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBIm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
bI_t̃_add = problem(gI₀ₜ,bI₀ᵤ,gO₀ₜ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirBIa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test bI_t̃_add ≈ 1 - 1/bI_t̃_mult

gO_t̃_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ᵤ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
gO_t̃_add = problem(gI₀ₜ,bI₀ₜ,gO₀ᵤ,bO₀ₜ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test gO_t̃_add ≈ (gO_t̃_mult -1)

bO_t̃_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
bO_t̃_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test bO_t̃_add ≈ 1 - 1/bO_t̃_mult


# Observation at time t, all other dmu at time t+1....
gI_ũ_mult = problem(gI₀ₜ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGIm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
gI_ũ_add = problem(gI₀ₜ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirGIa,startValues=(),forceLinearModel=true,crossTime=true,convexAssumption=convex)
@test  isapprox(gI_ũ_add,1 - 1/gI_ũ_mult, atol=0.000001)

bI_ũ_mult = problem(gI₀ᵤ,bI₀ₜ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBIm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
bI_ũ_add = problem(gI₀ᵤ,bI₀ₜ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirBIa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test isapprox(bI_ũ_add, 1 - 1/bI_ũ_mult, atol=0.0000001)

gO_ũ_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ₜ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
gO_ũ_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ₜ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test gO_ũ_add ≈ (gO_ũ_mult -1)

bO_ũ_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
bO_ũ_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true, crossTime=true,convexAssumption=convex)
@test bO_ũ_add ≈ 1 - 1/bO_ũ_mult

# Everything at time t+1 (u)....

gI_u_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=dirGIm,startValues=(),forceLinearModel=true,convexAssumption=convex)
gI_u_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=dirGIa,startValues=(),forceLinearModel=true,convexAssumption=convex)

bI_u_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=dirBIm,startValues=(),forceLinearModel=true,convexAssumption=convex)
bI_u_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=dirBIa,startValues=(),forceLinearModel=true,convexAssumption=convex)

gO_u_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true,convexAssumption=convex)
gO_u_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true,convexAssumption=convex)

bO_u_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true,convexAssumption=convex)
bO_u_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true,convexAssumption=convex)

# DMU measures at time t+1 and frontier at time t...
gI_ut_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=dirGIm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
gI_ut_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=dirGIa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

bI_ut_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
  directions=dirBIm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
bI_ut_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
  directions=dirBIa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

gO_ut_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
gO_ut_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

bO_ut_mult = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
bO_ut_add = problem(gI₀ᵤ,bI₀ᵤ,gO₀ᵤ,bO₀ᵤ,gIₜ,bIₜ,gOₜ,bOₜ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

# DMU measures at time t and frontier at time t+1...
gI_tu_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=dirGIm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
gI_tu_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=dirGIa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

bI_tu_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
  directions=dirBIm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
bI_tu_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
  directions=dirBIa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

gO_tu_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirGOm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
gO_tu_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirGOa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)

bO_tu_mult = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="multiplicative",
 directions=dirBOm,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)
bO_tu_add = problem(gI₀ₜ,bI₀ₜ,gO₀ₜ,bO₀ₜ,gIᵤ,bIᵤ,gOᵤ,bOᵤ,retToScale="variable",prodStructure="addittive",
 directions=dirBOa,startValues=(),forceLinearModel=true,convexAssumption=convex,crossTime=true)


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

idx_gi_ut,idx_bi_ut,idx_go_ut,idx_bo_ut = gI_ut_mult,bI_ut_mult,gO_ut_mult,bO_ut_mult
idx_gi_tu,idx_bi_tu,idx_go_tu,idx_bo_tu = gI_tu_mult,bI_tu_mult,gO_tu_mult,bO_tu_mult

1/idx_gi_ut,1/idx_bi_ut,1/idx_go_ut,1/idx_bo_ut
1/idx_gi_tu,1/idx_bi_tu,1/idx_go_tu,1/idx_bo_tu

idx_gi_ũ = gI_ũ_mult
idx_bi_ũ = bI_ũ_mult
idx_go_ũ = gO_ũ_mult
idx_bo_ũ = bO_ũ_mult


idx_i_t = (idx_gi_t̃/idx_gi_t) * (idx_bi_t̃/idx_bi_t)
idx_o_t = (idx_go_t/idx_go_t̃) * (idx_bo_t/idx_bo_t̃)
idx_t   = idx_o_t/idx_i_t
idx_i_u = (idx_gi_u/idx_gi_ũ) * (idx_bi_u/idx_bi_ũ)
idx_o_u = (idx_go_ũ/idx_go_u) * (idx_bo_ũ/idx_bo_u)
idx_u   = idx_o_u/idx_i_u
idx     = (idx_t * idx_u)^(1/2)


# Disaggregation  Good/ bads

idx_Gi_t = (idx_gi_t̃/idx_gi_t)
idx_Go_t = (idx_go_t/idx_go_t̃)
idx_Gt   = idx_Go_t/idx_Gi_t
idx_Gi_u = (idx_gi_u/idx_gi_ũ)
idx_Go_u = (idx_go_ũ/idx_go_u)
idx_Gu   = idx_Go_u/idx_Gi_u
idx_G    = (idx_Gt * idx_Gu)^(1/2)

idx_Bi_t = (idx_bi_t̃/idx_bi_t)
idx_Bo_t = (idx_bo_t̃/idx_bo_t)
idx_Bt   = 1/(idx_Bo_t * idx_Bi_t)
idx_Bi_u = (idx_bi_u/idx_bi_ũ)
idx_Bo_u = (idx_bo_u/idx_bo_ũ)
idx_Bu   = 1/(idx_Bo_u*idx_Bi_u)
idx_B    = (idx_Bt * idx_Bu)^(1/2)

# Disaggregation T/E/S...
idx_T_G_O = ((idx_go_t/idx_go_tu) * (idx_go_ut/idx_go_u) )^-(1/2)
idx_T_B_O = ((idx_bo_t/idx_bo_tu) * (idx_bo_ut/idx_bo_u) )^-(1/2)
idx_T_O   = idx_T_G_O * idx_T_B_O
idx_T_G_I = ((idx_gi_t/idx_gi_tu) * (idx_gi_ut/idx_gi_u) )^-(1/2)
idx_T_B_I = ((idx_bi_t/idx_bi_tu) * (idx_bi_ut/idx_bi_u) )^-(1/2)
idx_T_I   = idx_T_G_I * idx_T_B_I
idx_T     = (idx_T_O * idx_T_I)

idx_E_G_O = (idx_go_u/idx_go_t)^-1
idx_E_B_O = (idx_bo_u/idx_bo_t)^-1
idx_E_O   = idx_E_G_O * idx_E_B_O
idx_E_G_I = (idx_gi_u/idx_gi_t)^-1
idx_E_B_I = (idx_bi_u/idx_bi_t)^-1
idx_E_I   = idx_E_G_I * idx_E_B_I
idx_E     = (idx_E_O * idx_E_I)

idx_S_I   = idx / (idx_T_I * idx_E_I)
idx_S_O   = idx / (idx_T_O * idx_E_O)
idx_S_G_O = ((idx_go_t̃ / idx_go_ut) * (idx_gi_t̃ / idx_gi_t) * (idx_go_tu / idx_go_ũ) * (idx_gi_u/idx_gi_ũ)  )^-(1/2)
idx_S_B_O = ((idx_bo_t̃ / idx_bo_ut) * (idx_bi_t̃ / idx_bi_t) * (idx_bo_tu / idx_bo_ũ) * (idx_bi_u/idx_bi_ũ)  )^-(1/2)
@test idx_S_O ≈ idx_S_G_O * idx_S_B_O
@test idx ≈ idx_T_O * idx_E_O * idx_S_O

idx_S_G_I = ((idx_gi_t̃ / idx_gi_ut) * (idx_go_t̃ / idx_go_t) * (idx_gi_tu / idx_gi_ũ) * (idx_go_u/idx_go_ũ)  )^-(1/2)
idx_S_B_I = ((idx_bi_t̃ / idx_bi_ut) * (idx_bo_t̃ / idx_bo_t) * (idx_bi_tu / idx_bi_ũ) * (idx_bo_u/idx_bo_ũ)  )^-(1/2)
@test idx_S_I ≈ idx_S_G_I * idx_S_B_I
@test idx ≈ idx_T_I * idx_E_I * idx_S_I


idx_S_G_O = idx_G /  (idx_T_G_O * idx_E_G_O)
idx_S_B_O = idx_B /  (idx_T_B_O * idx_E_B_O)
idx_S_O   = idx_S_G_O * idx_S_B_O
idx_S_G_I = idx_G /  (idx_T_G_I * idx_E_G_I)
idx_S_B_I = idx_B /  (idx_T_B_I * idx_E_B_I)
idx_S_I   = idx_S_G_I * idx_S_B_I
idx_S     = (idx_S_O * idx_S_I)^(1/2)


# Addittive
(idx_gi_t, idx_bi_t, idx_go_t, idx_bo_t) = (gI_add, bI_add, gO_add, bO_add)
(idx_gi_t̃, idx_bi_t̃, idx_go_t̃, idx_bo_t̃) = (gI_t̃_add, bI_t̃_add, gO_t̃_add, bO_t̃_add)
(idx_gi_u, idx_bi_u, idx_go_u, idx_bo_u) = (gI_u_add, bI_u_add, gO_u_add, bO_u_add)
(idx_gi_ũ, idx_bi_ũ, idx_go_ũ, idx_bo_ũ) = (gI_ũ_add, bI_ũ_add, gO_ũ_add, bO_ũ_add)
idx_gi_ut,idx_bi_ut,idx_go_ut,idx_bo_ut = gI_ut_add,bI_ut_add,gO_ut_add,bO_ut_add
idx_gi_tu,idx_bi_tu,idx_go_tu,idx_bo_tu = gI_tu_add,bI_tu_add,gO_tu_add,bO_tu_add



idx_i_t = (idx_gi_t̃ - idx_gi_t) + (idx_bi_t̃ - idx_bi_t)
idx_o_t = (idx_go_t - idx_go_t̃) - (idx_bo_t̃ - idx_bo_t)
idx_t   = idx_o_t - idx_i_t
idx_i_u = (idx_gi_u - idx_gi_ũ) + (idx_bi_u - idx_bi_ũ)
idx_o_u = (idx_go_ũ - idx_go_u) - (idx_bo_u - idx_bo_ũ)
idx_u   = idx_o_u - idx_i_u
idx     = (idx_t + idx_u) / 2




idx_T_G_O = ((-idx_go_t+idx_go_tu) + (-idx_go_ut+idx_go_u) )/2
idx_T_B_O = ((-idx_bo_t-idx_bo_tu) + (+idx_bo_ut+idx_bo_u) )/2
idx_T_O   = idx_T_G_O + idx_T_B_O
idx_T_G_I = ((-idx_gi_t+idx_gi_tu) + (-idx_gi_ut+idx_gi_u) )/2
idx_T_B_I = ((-idx_bi_t+idx_bi_tu) + (-idx_bi_ut+idx_bi_u) )/2
idx_T_I   = idx_T_G_I + idx_T_B_I
idx_T     = (idx_T_O + idx_T_I)

idx_E_G_O = -(idx_go_u-idx_go_t)
idx_E_B_O = -(idx_bo_u-idx_bo_t)
idx_E_O   = idx_E_G_O + idx_E_B_O
idx_E_G_I = -(idx_gi_u-idx_gi_t)
idx_E_B_I = -(idx_bi_u-idx_bi_t)
idx_E_I   = idx_E_G_I + idx_E_B_I
idx_E     = (idx_E_O + idx_E_I)

idx_S_I   = idx - (idx_T_I + idx_E_I)
idx_S_O   = idx - (idx_T_O + idx_E_O)
idx_S_G_O = -((idx_go_t̃ - idx_go_ut) + (idx_gi_t̃ - idx_gi_t) + (idx_go_tu - idx_go_ũ) + (idx_gi_u-idx_gi_ũ)  )/2
idx_S_B_O = -((idx_bo_t̃ + idx_bo_ut) + (idx_bi_t̃ - idx_bi_t) + (-idx_bo_tu - idx_bo_ũ) + (idx_bi_u-idx_bi_ũ)  )/2
idx_S_G_I = -((idx_gi_t̃ - idx_gi_ut) + (idx_go_t̃ - idx_go_t) + (idx_gi_tu - idx_gi_ũ) + (idx_go_u-idx_go_ũ)  )/2
idx_S_B_I = -((idx_bi_t̃ - idx_bi_ut) + (idx_bo_t̃ - idx_bo_t) + (idx_bi_tu - idx_bi_ũ) + (idx_bo_u-idx_bo_ũ)  )/2






oecdAnalysis  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=true)


isapprox(oecdAnalysis.prodIndexes_G .* oecdAnalysis.prodIndexes_B, oecdAnalysis.prodIndexes, atol=0.000001)
#isapprox(oecdAnalysis.prodIndexes_T .* oecdAnalysis.prodIndexes_E .* oecdAnalysis.prodIndexes_S, oecdAnalysis.prodIndexes, atol=0.000001)


oecdAnalysis.prodIndexes

oecdAnalysis.prodIndexes_T


oecdAnalysis.prodIndexes_E


oecdAnalysis.prodIndexes_S













oecdAnalysisA  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="addittive",convexAssumption=true)


oecdAnalysisA.prodIndexes_T
oecdAnalysisA.prodIndexes_E
oecdAnalysisA.prodIndexes_S





isapprox(oecdAnalysisA.prodIndexes_G .+ oecdAnalysisA.prodIndexes_B, oecdAnalysisA.prodIndexes, atol=0.000001)
#isapprox(oecdAnalysisA.prodIndexes_T .+ oecdAnalysisA.prodIndexes_E .+ oecdAnalysisA.prodIndexes_S, oecdAnalysisA.prodIndexes, atol=0.000001)


# Non convex test

oecdAnalysis_nc  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=false)


isapprox(oecdAnalysis_nc.prodIndexes_G .* oecdAnalysis_nc.prodIndexes_B, oecdAnalysis_nc.prodIndexes, atol=0.000001)
isapprox(oecdAnalysis_nc.prodIndexes_T .* oecdAnalysis_nc.prodIndexes_E .* oecdAnalysis_nc.prodIndexes_S, oecdAnalysis_nc.prodIndexes, atol=0.000001)




add = oecdAnalysisA.prodIndexes

#=
mBool = mult .> (1 - myeps)
aBool = add .> (0 - myeps)
mBool == aBool
mBool .== aBool
=#
