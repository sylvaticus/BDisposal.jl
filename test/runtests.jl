
using Test, DataFrames, CSV, BDisposal
println("Testing BDisposal...")

# Aitport data test with only one input category...
airportData = CSV.read(joinpath(@__DIR__,"data","airports.csv"),DataFrame; delim=';',copycols=true)
airportGoodInputs  = ["employees","totalCosts"]
airportBadInputs   = []
airportGoodOutputs = ["passengers"]
airportBadOutputs  = ["co2emissions"]
sort!(airportData, [:period, :dmu]) # sort data by period and dmu
periods = unique(airportData.period)
dmus    = unique(airportData.dmu)

nGI, nBI, nGO, nBO, nPer, nDMUs,  = length(airportGoodInputs), length(airportBadInputs), length(airportGoodOutputs), length(airportBadOutputs), length(periods),length(dmus)

gI = Array{Float64}(undef, (nDMUs,nGI,nPer))
bI = Array{Float64}(undef, (nDMUs,nBI,nPer))
gO = Array{Float64}(undef, (nDMUs,nGO,nPer))
bO = Array{Float64}(undef, (nDMUs,nBO,nPer))

for (p,period) in enumerate(periods)
    periodData = airportData[airportData.period .== period,:]
    gI[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodInputs])
    if nBI > 0
         bI[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadInputs])
    end
    gO[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodOutputs])
    bO[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadOutputs])
end

# Call the function to get the efficiency measurements for constant returns to scale
(λ_crs, λ_convex_crs, λ_nonconvex_crs, nonConvTest_value_crs, nonConvTest_crs) = efficiencyScores(
gI,gO,bO,bI,retToScale="constant", dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="multiplicative")
@test nonConvTest_value_crs[3,2]  ≈ 1.1587129278170434

(λ_vrs, λ_convex_vrs, λ_nonconvex_vrs, nonConvTest_value_vrs, nonConvTest_vrs) = efficiencyScores(
gI,gO,bO,bI,retToScale="variable", dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="additive")
@test nonConvTest_value_vrs[3,3]  ≈ 7.432043216538459

# Aitport data test with bad inputs separated...
airportData = CSV.read(joinpath(@__DIR__,"data","airports.csv"),DataFrame; delim=';',copycols=true)
airportGoodInputs  = ["employees"]
airportBadInputs   = ["totalCosts"]
airportGoodOutputs = ["passengers"]
airportBadOutputs  = ["co2emissions"]
sort!(airportData, [:period, :dmu]) # sort data by period and dmu
periods = unique(airportData.period)
dmus    = unique(airportData.dmu)

nGI, nBI, nGO, nBO, nPer, nDMUs,  = length(airportGoodInputs), length(airportBadInputs), length(airportGoodOutputs), length(airportBadOutputs), length(periods),length(dmus)

gI = Array{Float64}(undef, (nDMUs,nGI,nPer))
bI = Array{Float64}(undef, (nDMUs,nBI,nPer))
gO = Array{Float64}(undef, (nDMUs,nGO,nPer))
bO = Array{Float64}(undef, (nDMUs,nBO,nPer))

for (p,period) in enumerate(periods)
    periodData = airportData[airportData.period .== period,:]
    gI[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodInputs])
    if nBI > 0
        bI[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadInputs])
    end
    gO[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodOutputs])
    bO[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadOutputs])
end

prodIndices = prodIndex(gI,gO,bO,bI;
                   retToScale="constant",prodStructure="multiplicative",convexAssumption=true,
                   startθ=0,startμ=0,startλ=1.1)

@test prodIndices[3,2] ≈ 1.147467571896574

prodIndices = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=true,
                   startθ=0,startμ=0,startλ=1.1)

prodIndices = prodIndex(gI,gO,bO,bI;
                  retToScale="variable",prodStructure="additive",convexAssumption=true,
                  startθ=0,startμ=0,startλ=1.1)
@test prodIndices[3,2] ≈ -0.09534201521261434
