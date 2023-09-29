using Test, DataFrames, CSV, BDisposal
println("Testing BDisposal...")

println("Testing BDisposal on js-data...")
# Byung M. Jeon, Robin C. Sickles, "The Role of Environmental Factors in
# Growth Accounting: a Nonparametric Analysis", Journal of Applied
# Econometrics, Vol. 19, No. 5, 2004, pp. 567-591.

data = CSV.read(joinpath(dirname(pathof(BDisposal)),"..","test","data","js-data","oecd.txt"),DataFrame; delim=' ',ignorerepeated=true,copycols=true,header=false)
rename!(data,[:ccid,:year,:gdp,:co2,:capital,:labour,:energy])
dmuMap = Dict(1 => "Canada", 2 => "the United States", 3 => "Japan", 4 => "Austria",
              5 => "Belgium", 6 => "Denmark", 7 => "Finland", 8 => "France", 9 => "Germany",
              10 => "Greece", 11 => "Ireland", 12 => "Italy", 13 => "Norway", 14 => "Spain",
              15 => "Sweden", 16 => "U.K.", 17 => "Australia")
data.ccid = map(x->dmuMap[x], data.ccid)     
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

################################################################################
# Testing prodIndex

oecdAnalysis  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=true)

@test oecdAnalysis.prodIndexes[4,10] == 1.02891196267105


@test isapprox(oecdAnalysis.prodIndexes_G .* oecdAnalysis.prodIndexes_B, oecdAnalysis.prodIndexes, atol=0.000001)
decomp = isapprox.(oecdAnalysis.prodIndexes_T .* oecdAnalysis.prodIndexes_E .* oecdAnalysis.prodIndexes_S, oecdAnalysis.prodIndexes, atol=0.000001)
@test all(skipmissing(decomp)) == true


oecdAnalysisA  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="addittive",convexAssumption=true)

@test oecdAnalysisA.prodIndexes[4,10] == 0.02111291275337509
@test isapprox(oecdAnalysisA.prodIndexes_G .+ oecdAnalysisA.prodIndexes_B, oecdAnalysisA.prodIndexes, atol=0.000001)
decomp =  isapprox.(oecdAnalysisA.prodIndexes_T .+ oecdAnalysisA.prodIndexes_E .+ oecdAnalysisA.prodIndexes_S, oecdAnalysisA.prodIndexes, atol=0.000001)
@test all(skipmissing(decomp)) == true

oecdAnalysis_nc  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=false)

@test isapprox(oecdAnalysis_nc.prodIndexes_G .* oecdAnalysis_nc.prodIndexes_B, oecdAnalysis_nc.prodIndexes, atol=0.000001)

oecdAnalysis_ncA  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="additive",convexAssumption=false)

@test isapprox(oecdAnalysis_ncA.prodIndexes_G .+ oecdAnalysis_ncA.prodIndexes_B, oecdAnalysis_ncA.prodIndexes, atol=0.000001)

@test isapprox(log.(oecdAnalysis_nc.prodIndexes), oecdAnalysis_ncA.prodIndexes, atol=0.05)

oecdAnalysisFB  = prodIndexFB(gI,gO,bO,bI;remarcable_obs_dmu=17, remarcable_obs_period=10)

@test oecdAnalysisFB.prodIndexes[17,10] ≈ 1.0

@test oecdAnalysisFB.prodIndexes[1,1] ≈ 0.000990746731264865

oecdAnalysisFB2  = prodIndexFB(gI,gO,bO,bI;remarcable_obs_dmu=4, remarcable_obs_period=4)
@test oecdAnalysisFB2.prodIndexes[4,4] ≈ 1.0
@test oecdAnalysisFB2.prodIndexes[17,4] ≈ 94.60766200738495

@test all(.! ismissing.(oecdAnalysisFB2.prodIndexes)) # none is missing


################################################################################
# Testing efficiencyScore using the Airports dataset..


println("Testing BDisposal on airport data...")
airportData = CSV.read(joinpath(dirname(pathof(BDisposal)),"..","test","data","airports.csv"),DataFrame; delim=';',copycols=true)
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
    gI[:,:,p] = Matrix{Float64}(periodData[:,airportGoodInputs])
    if nBI > 0
         bI[:,:,p] = Matrix{Float64}(periodData[:,airportBadInputs])
    end
    gO[:,:,p] = Matrix{Float64}(periodData[:,airportGoodOutputs])
    bO[:,:,p] = Matrix{Float64}(periodData[:,airportBadOutputs])
end

(λ_crs, λ_convex_crs, λ_nonconvex_crs, nonConvTest_crs, nonConvTest_value_crs) = efficiencyScores(
  gI,gO,bO,bI,retToScale="variable", dirGI=0,dirBI=0,dirGO=1,dirBO=0, prodStructure="multiplicative")
@test nonConvTest_value_crs[3,2]  ≈ 0.8842007631310966


################################################################################
# Basic testing of dmuEfficiency

println("Testing of the basic dmuEfficiency function..")

I  = [10 2; 8 4; 12 1.5; 24 3]
O =  [100;80;120;120]
I₀ = [24 3]
O₀ = [120]

results = dmuEfficiency(I₀,O₀,I,O)


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
  5	6 ]
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
  57	60	40]

nDMU = size(I,1)
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O)[:obj] for  d in 1:nDMU]
efficiencies = hcat(1:nDMU,efficiencies)
efficiencies = efficiencies[sortperm(efficiencies[:, 2],rev=true), :]

# Test treating last output as bad
O2 = convert(Array{Float64,2},copy(O))
O2[:,3] = 1 ./ O2[:,3] # third is a bad output
efficiencies = [dmuEfficiency(I[d,:],O2[d,:],I,O2)[:obj] for  d in 1:nDMU]
efficiencies = hcat(1:nDMU,efficiencies)
efficiencies = efficiencies[sortperm(efficiencies[:, 2],rev=true), :]


#Table 13.1 Cooper 2006
I = ones(Float64,9)
O = [
  1 1
  2 1
  6 2
  8 4
  9 7
  5 2
  4 3
  6 4
  4 6
]
nDMU = size(I,1)
O2 = convert(Array{Float64,2},copy(O))
O2[:,2] = 1 ./ O2[:,2] # third is a bad output

efficiencies = [dmuEfficiency(I[d,:],O2[d,:],I,O2)[:obj] for  d in 1:nDMU]
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
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O)[:obj] for  d in 1:nDMU]
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
 1]

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
efficiencies = [dmuEfficiency(I[d,:],O[d,:],I,O)[:obj] for  d in 1:nDMU]
#efficiencies = hcat(1:nDMU,efficiencies)
#scatter(O[:,1],O[:,2])
@test efficiencies == [1.0,1.0,1.0,1.0,1.0,1.0,0.37735849056603776,1.0]

# From Example 2.2 of Cooper et oth. 2006 "Data envelopment Analysis. A Comprehensive Text with...."

X = [4 7 8 4 2 10; 3 3 1 2 4 1]
Y = [1 1 1 1 1 1]

nDMU = size(X,2)
out = [dmuEfficiency(X[:,d],Y[:,d],X',Y') for  d in 1:nDMU]
@test out[1].obj ≈ 0.8571428571428571
@test out[1].wI ≈ [0.14285714285714285, 0.14285714285714285]
@test out[1].wO ≈ [0.8571428571428571]
@test collect(keys(out[1].refSet)) ≈ [5,4] || collect(keys(out[1].refSet)) ≈ [4,5]
@test collect(values(out[1].refSet)) ≈ [0.2857142857142857,0.7142857142857143] ||  collect(values(out[1].refSet)) ≈ [0.7142857142857143,0.2857142857142857]
@test [i[:eff] for i in out] == [false,false,true,true,true,false]

outDual = [dmuEfficiencyDual(X[:,d],Y[:,d],X',Y') for  d in 1:nDMU]
@test [i[:eff] for i in outDual] == [false,false,true,true,true,false]


# From Example 2.2 of Cooper et oth. 2006 "Data envelopment Analysis. A Comprehensive Text with...."

X = [4 7 8 4 2 10; 3 3 1 2 4 1]
Y = [1 1 1 1 1 1]

nDMU = size(X,2)
out = [dmuEfficiency(X[:,d],Y[:,d],X',Y') for  d in 1:nDMU]
@test out[1].obj ≈ 0.8571428571428571
@test out[1].wI ≈ [0.14285714285714285, 0.14285714285714285]
@test out[1].wO ≈ [0.8571428571428571]
@test collect(keys(out[1].refSet)) ≈ [5,4] || collect(keys(out[1].refSet)) ≈ [4,5]
@test collect(values(out[1].refSet)) ≈ [0.2857142857142857,0.7142857142857143] ||  collect(values(out[1].refSet)) ≈ [0.7142857142857143,0.2857142857142857]
@test [i[:eff] for i in out] == [false,false,true,true,true,false]

outDual = [dmuEfficiencyDual(X[:,d],Y[:,d],X',Y') for  d in 1:nDMU]
@test [i[:eff] for i in outDual] == [false,false,true,true,true,false]
