
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

# Setting empty containers for our data
# Each of them is a 3D matrix where the first dimension is the decision units, the second one is the individual input or output item and the third dimension is the period to which the data refer
gI = Array{Float64}(undef, (nDMUs,nGI,nPer)) # Good inputs
bI = Array{Float64}(undef, (nDMUs,nBI,nPer)) # Bad inputs (optional)
gO = Array{Float64}(undef, (nDMUs,nGO,nPer)) # Good outputs, aka "desiderable" outputs
bO = Array{Float64}(undef, (nDMUs,nBO,nPer)) # Bad outputs, aka "undesiderable" outputs

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
(λ_crs, λ_convex_crs, λ_nonconvex_crs, nonConvTest_crs, nonConvTest_value_crs) = efficiencyScores(
gI,gO,bO,bI,retToScale="constant", dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="multiplicative")
@test nonConvTest_value_crs[3,2]  ≈ 1.1587129278170434

(λ_vrs, λ_convex_vrs, λ_nonconvex_vrs, nonConvTest_vrs, nonConvTest_value_vrs) = efficiencyScores(
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

# Basic testing of dmuEfficiency


I  = [10 2; 8 4; 12 1.5; 24 3]
O =  [100;80;120;120]
I₀ = [24 3]
O₀ = [120]

dmuEfficiency(I₀,O₀,I,O)
I = [
3	5
2.5	4.5
4	6
6	7
2.3	3.5
4	6.5
7	10
4.4	6.4
3	5
5	7
5	7
2	4
5	7
4	4
2	3
3	6
7	11
4	6
3	4
5	6
]
O = [
40	55	30
45	50	40
55	45	30
48	20	60
28	50	25
48	20	65
80	65	57
25	48	30
45	64	42
70	65	48
45	65	40
45	40	44
65	25	35
38	18	64
20	50	15
38	20	60
68	64	54
25	38	20
45	67	32
57	60	40
]
nDMU = size(I,1)
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O) for  d in 1:nDMU]
efficiencies = hcat(1:nDMU,efficiencies)
efficiencies = efficiencies[sortperm(efficiencies[:, 2],rev=true), :]

I = [
4	140
5	90
6	36
10	300
11	66
8	36
9	12
5	210
5.5	33
8	288
10	80
8	8
]

O =[
2	28
1	22.5
6	12
8	60
7	16.5
6	12
7	6
3	30
4.4	5.5
4	72
2	20
1	4
]

nDMU = size(I,1)
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O) for  d in 1:nDMU]
efficiencies = hcat(1:nDMU,efficiencies)
efficiencies = efficiencies[sortperm(efficiencies[:, 2],rev=true), :]

I = [
1
1
1
1
1
1
1
1
]

O = [
1 7
2 7
4 6
6 4
7 2
7 1
2 2
5.3 5.3
]
nDMU = size(I,1)
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O) for  d in 1:nDMU]
#efficiencies = hcat(1:nDMU,efficiencies)
scatter(O[:,1],O[:,2])
