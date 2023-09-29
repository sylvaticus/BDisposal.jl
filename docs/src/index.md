# ![BDisposalLogos](assets/BDisposal_logo_30x30.png) BDisposal.jl Documentation

Welcome to the documentation for the [_Non parametric productivity analysis under the B-disposal assumption_](https://github.com/sylvaticus/BDisposal.jl) package.

The BDisposal package proposes a serie of environmental efficiency and productivity algorithms for non-parametric modelling when we relax the disposability assumption of some of the outputs and/or inputs (e.g. pollution). These efficiency and productivity measures are implemented through convex and non convex [Data Envelopment Analysis (DEA)](https://en.wikipedia.org/wiki/Data_envelopment_analysis) (aka _Frontier Efficiency Analysis_) models.


## Installation

Install the BDisposal package with:

* `] add BDisposal.jl`

## Loading the module(s)

```julia
using BDisposal
```

## Usage

The BDisposal package contains two prominent components.
First, the BDisposal package defines environmental efficiency indicators for a set of Decision Making Units (`efficiencyScores()`), computing the distance to the best environmental production procedures, i.e., to the efficient production frontier.
Second, the BDisposal package displays environmental productivity indices (`prodIndex()`). These productivity measures are implemented for different time periods (eg., years, months etc.) or spatial units (eg., countries, cities etc.), based on the aforementioned environmental efficiency indicators.

Both components can consider constant or variable returns to scale, convex or non-convex production frontiers and additive or multiplicative distance functions, and work with multiple "good" ("desirable") and "bad" ("undesirable") inputs and outputs (with undesirable inputs being optional) (not all combinations are implemented, see the individual functions for details and limitations).

They are detailed in their respective pages:

- [`efficiencyScores`](@ref): Compute efficiency indicators and convexity test results;
- [`prodIndex`](@ref): Compute productivity indexes;
- [`prodIndexFB`](@ref): Compute productivity indexes under a provided fixed base;

The package provides also a couple of functions to compute individual DMU problem using "vanilla" DEA (without considering the beta disposability assumption), [`dmuEfficiency`](@ref) and [`dmuEfficiencyDual`](@ref).

## Examples

### Airport example

The following example (with data from _Abad and Briec, 2019_) show how to compute the efficiency indicators for 14 main French airports, where the (good) inputs considered are _employees_ and _totalCosts_ and the (single in this case) bad and good outputs are respectively (airport) _co2emissions_ and _passengers_.


```julia
using DataFrames, CSV, BDisposal
println("Estimating French airports efficiency and productivity indexes...")

# Loading airport data and preparing the data..
# Data is stored as a  CSV files with 6 columns: period, dmu, co2emissions, passengers, employees, totalCosts
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
# Transferring data to the containers
for (p,period) in enumerate(periods)
    periodData = airportData[airportData.period .== period,:]
    gI[:,:,p] = Matrix(periodData[:,airportGoodInputs])
    if nBI > 0
         bI[:,:,p] = Matrix(periodData[:,airportBadInputs])
    end
    gO[:,:,p] = Matrix(periodData[:,airportGoodOutputs])
    bO[:,:,p] = Matrix(periodData[:,airportBadOutputs])
end

# Call the function to get the efficiency measurements for constant returns to scale
(λ, λ_convex, λ_nonconvex, nonConvTest, nonConvTest_value) = efficiencyScores(
  gI,gO,bO,bI,retToScale="constant", dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="multiplicative")


# Add periods as headers and decision making names as first column in order to show the data
# Efficiency Indexes
λ = hcat(dmus,λ)
λdf = DataFrame(λ , Symbol.(vcat("DMU",periods)))
```

```
│ Row │ DMU                      │ 2007    │ 2008    │ 2009    │ 2010    │ 2011    │
│     │ Any                      │ Any     │ Any     │ Any     │ Any     │ Any     │
├─────┼──────────────────────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ Beauvais                 │ 1.0     │ 1.0     │ 1.0     │ 1.0     │ 1.0     │
│ 2   │ Bordeaux-Mérignac        │ 1.0     │ 1.0     │ 1.0     │ 1.0     │ 1.0     │
│ 3   │ Bâle-Mulhouse            │ 1.08295 │ 1.08777 │ 1.17694 │ 1.15278 │ 1.0     │
│ 4   │ Lille                    │ 1.13183 │ 1.15032 │ 1.06293 │ 1.0     │ 1.0     │
│ 5   │ Lyon-Saint Exupéry       │ 1.09567 │ 1.07909 │ 1.06901 │ 1.10468 │ 1.00911 │
│ 6   │ Marseille-Provence       │ 1.00071 │ 1.04287 │ 1.0     │ 1.0     │ 1.0     │
│ 7   │ Montpellier-Méditerranée │ 1.0     │ 1.0     │ 1.0     │ 1.10249 │ 1.14045 │
│ 8   │ Nantes-Atlantique        │ 1.09003 │ 1.05685 │ 1.0013  │ 1.05295 │ 1.07531 │
│ 9   │ Nice-Côte d'azur         │ 1.0     │ 1.02201 │ 1.02562 │ 1.08272 │ 1.02654 │
│ 10  │ Paris CDG                │ 1.14695 │ 1.14915 │ 1.15329 │ 1.17032 │ 1.17137 │
│ 11  │ Paris ORY                │ 1.21574 │ 1.22947 │ 1.21526 │ 1.22424 │ 1.21444 │
│ 12  │ Strasbourg-Entzheim      │ 1.0     │ 1.14408 │ 1.13534 │ 1.23498 │ 1.27062 │
│ 13  │ Toulouse-Blagnac         │ 1.0     │ 1.0     │ 1.0     │ 1.01143 │ 1.0     │
```
```julia
# Non-convexity test
nc_test = hcat(dmus,nonConvTest)
nc_test_df = DataFrame(nc_test , Symbol.(vcat("DMU",periods)))
```

```
│ Row │ DMU                      │ 2007 │ 2008 │ 2009 │ 2010 │ 2011 │
│     │ Any                      │ Any  │ Any  │ Any  │ Any  │ Any  │
├─────┼──────────────────────────┼──────┼──────┼──────┼──────┼──────┤
│ 1   │ Beauvais                 │ 1    │ 1    │ 1    │ 1    │ 1    │
│ 2   │ Bordeaux-Mérignac        │ 1    │ 1    │ 1    │ 1    │ 1    │
│ 3   │ Bâle-Mulhouse            │ 0    │ 0    │ 0    │ 0    │ 1    │
│ 4   │ Lille                    │ 0    │ 0    │ 0    │ 1    │ 1    │
│ 5   │ Lyon-Saint Exupéry       │ 1    │ 1    │ 1    │ 0    │ 1    │
│ 6   │ Marseille-Provence       │ 1    │ 0    │ 1    │ 1    │ 1    │
│ 7   │ Montpellier-Méditerranée │ 1    │ 1    │ 1    │ 1    │ 1    │
│ 8   │ Nantes-Atlantique        │ 0    │ 0    │ 1    │ 1    │ 1    │
│ 9   │ Nice-Côte d'azur         │ 1    │ 1    │ 1    │ 1    │ 1    │
│ 10  │ Paris CDG                │ 0    │ 0    │ 0    │ 0    │ 0    │
│ 11  │ Paris ORY                │ 0    │ 0    │ 0    │ 0    │ 0    │
│ 12  │ Strasbourg-Entzheim      │ 1    │ 1    │ 1    │ 1    │ 1    │
│ 13  │ Toulouse-Blagnac         │ 1    │ 1    │ 1    │ 0    │ 1    │
```

### OECD Example
Here we compute the productivity indexes of various OECD countries (with data from
[Jeon and Sickles (2004)](https://doi.org/10.1002/jae.769)) in terms of gdp growth over traditional production inputs as capital and labour, but also considering CO2 emissions and energy use.

```julia
using DataFrames, CSV, BDisposal

# Loading data and formatting them in the way required by `prodIndex`
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

# Performing the analysis
oecdAnalysis  = prodIndex(gI,gO,bO,bI;
                   retToScale="variable",prodStructure="multiplicative",convexAssumption=true)

# Showing production indexes for all countries..
pIdx = oecdAnalysis.prodIndexes

# Add periods as headers and decision making names as first column in order to show the data
# Efficiency Indexes
pIdx  = hcat(dmus,pIdx)
pIdxDf = DataFrame(pIdx, Symbol.(vcat("Country",periods[2:end])))
```
```
│ Row │ Country           │ 1981     │ 1982     │ 1983     │ 1984     │ 1985     │ 1986     │ 1987     │ 1988     │ 1989     │ 1990     │
│     │ Any               │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │
├─────┼───────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│ 1   │ Australia         │ 0.980615 │ 0.874619 │ 1.07314  │ 0.946456 │ 0.913121 │ 0.971135 │ 0.952049 │ 0.963533 │ 0.891746 │ 0.869411 │
│ 2   │ Austria           │ 0.950153 │ 1.13918  │ 1.08687  │ 0.861064 │ 0.996057 │ 0.96336  │ 1.04583  │ 1.02559  │ 1.00891  │ 0.930929 │
│ 3   │ Belgium           │ 0.997489 │ 1.12907  │ 1.18585  │ 0.964121 │ 0.957716 │ 0.985223 │ 1.0288   │ 1.06768  │ 0.972852 │ 0.957646 │
│ 4   │ Canada            │ 1.08493  │ 1.02772  │ 1.01386  │ 0.964103 │ 0.931114 │ 1.02906  │ 0.959126 │ 0.877363 │ 0.938174 │ 1.02891  │
│ 5   │ Denmark           │ 1.13549  │ 1.17356  │ 1.13089  │ 0.961468 │ 0.811909 │ 1.03924  │ 0.949642 │ 1.08736  │ 1.00549  │ 1.08056  │
│ 6   │ Finland           │ 1.03238  │ 1.24164  │ 1.06111  │ 1.02071  │ 0.911113 │ 0.979476 │ 0.919278 │ 1.06038  │ 0.987295 │ 0.947106 │
│ 7   │ France            │ 1.1395   │ 1.09728  │ 0.999966 │ 1.03334  │ 0.960895 │ 1.08548  │ 1.0046   │ 1.10412  │ 0.90306  │ 0.9627   │
│ 8   │ Germany           │ 1.07873  │ 1.02427  │ 1.00385  │ 0.987006 │ 0.983202 │ 0.967011 │ 1.03386  │ 1.03729  │ 1.03669  │ 0.981721 │
│ 9   │ Greece            │ 0.975428 │ 1.04479  │ 0.880086 │ 0.941675 │ 0.910409 │ 0.937702 │ 0.851516 │ 0.885416 │ 0.99565  │ 0.892712 │
│ 10  │ Ireland           │ 1.07683  │ 1.06534  │ 0.974727 │ 1.07845  │ 0.902097 │ 0.824995 │ 1.08184  │ 1.00502  │ 0.974711 │ 0.948907 │
│ 11  │ Italy             │ 0.966255 │ 1.02578  │ 1.02878  │ 0.980468 │ 0.952348 │ 1.00805  │ 0.894706 │ 1.021    │ 0.931195 │ 1.00194  │
│ 12  │ Japan             │ 1.02906  │ 1.11366  │ 1.04336  │ 0.869509 │ 1.02776  │ 1.02243  │ 0.971448 │ 0.90699  │ 0.944852 │ 0.94903  │
│ 13  │ Norway            │ 1.13115  │ 1.1108   │ 0.949576 │ 0.936553 │ 0.904088 │ 0.956691 │ 0.859323 │ 1.16245  │ 0.898815 │ 1.02487  │
│ 14  │ Spain             │ 0.915789 │ 0.883303 │ 0.996596 │ 1.01146  │ 1.03484  │ 0.991454 │ 1.01437  │ 0.921896 │ 0.849651 │ 1.08936  │
│ 15  │ Sweden            │ 1.13215  │ 1.2188   │ 1.04962  │ 1.0498   │ 0.878779 │ 1.0206   │ 0.958626 │ 1.00718  │ 1.09563  │ 1.03445  │
│ 16  │ U.K.              │ 1.04728  │ 1.09075  │ 0.996706 │ 0.998393 │ 0.966998 │ 1.00421  │ 1.00086  │ 1.04238  │ 0.944179 │ 0.9494   │
│ 17  │ the United States │ 1.05193  │ 1.04846  │ 1.02621  │ 0.955054 │ 1.00271  │ 0.995934 │ 0.940216 │ 0.929503 │ 0.980136 │ 0.964614 │
```

```julia
# Focusing on Austria decomposition
AustriaPIdx = vcat(oecdAnalysis.prodIndexes[2,:]',
     oecdAnalysis.prodIndexes_G[2,:]',
     oecdAnalysis.prodIndexes_B[2,:]',
     oecdAnalysis.prodIndexes_T[2,:]',
     oecdAnalysis.prodIndexes_E[2,:]',
     oecdAnalysis.prodIndexes_S[2,:]')

AustriaPIdx   = hcat(["Overall production indexes",
              "Decomposition for \"good\" inputs and outputs",
              "Decomposition for \"bad\" inputs and outputs",
              "Decomposition for the technological component",
              "Decomposition for the efficiency component",
              "Decomposition for the scale (residual) component",
              ],AustriaPIdx)


AustriaPIdxDf  = DataFrame(AustriaPIdx, Symbol.(vcat("Item",periods[2:end])))
```

```
│ Row │ Item                                             │ 1981     │ 1982    │ 1983     │ 1984     │ 1985     │ 1986     │ 1987     │ 1988     │ 1989     │ 1990     │
│     │ Any                                              │ Any      │ Any     │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │ Any      │
├─────┼──────────────────────────────────────────────────┼──────────┼─────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│ 1   │ Overall production indexes                       │ 0.950153 │ 1.13918 │ 1.08687  │ 0.861064 │ 0.996057 │ 0.96336  │ 1.04583  │ 1.02559  │ 1.00891  │ 0.930929 │
│ 2   │ Decomposition for "good" inputs and outputs      │ 0.969768 │ 0.99197 │ 1.01057  │ 1.0042   │ 1.01544  │ 1.00682  │ 1.01592  │ 1.03277  │ 1.03173  │ 1.02555  │
│ 3   │ Decomposition for "bad" inputs and outputs       │ 0.979773 │ 1.1484  │ 1.0755   │ 0.857465 │ 0.980912 │ 0.956835 │ 1.02945  │ 0.99305  │ 0.977884 │ 0.907733 │
│ 4   │ Decomposition for the technological component    │ 1.0192   │ 1.0416  │ 1.06671  │ 1.15647  │ 0.96198  │ 1.03854  │ 1.02924  │ 1.1261   │ 1.13558  │ 1.07684  │
│ 5   │ Decomposition for the efficiency component       │ 0.868255 │ 1.07801 │ 1.06762  │ 0.75353  │ 1.08088  │ 0.955791 │ 1.09237  │ 0.948266 │ 0.976689 │ 0.944898 │
│ 6   │ Decomposition for the scale (residual) component │ 1.0737   │ 1.01453 │ 0.954357 │ 0.988098 │ 0.957945 │ 0.970511 │ 0.930198 │ 0.960434 │ 0.909664 │ 0.914914 │
```

# Production indexes considering a fixed base

Similarly, we can compute the production index with reference to a fixed base, both in term of specific DMU and period.
For example, if we consider as reference base the US observation for 1989 we have:

```julia
# Performing the analysis
oecdAnalysisFB  = prodIndexFB(gI,gO,bO,bI;remarcable_obs_dmu=17,remarcable_obs_period=10);
pIdx = oecdAnalysisFB.prodIndexes;

# Add periods as headers and decision making names as first column in order to show the data
# Efficiency Indexes
pIdx  = hcat(dmus,pIdx);
pIdxDf = DataFrame(pIdx, Symbol.(vcat("Country",periods)))
```
```
 Row │ Country            1980         1981         1982         1983         1984         1985         1986         1987         1988         1989         1990        
     │ Any                Any          Any          Any          Any          Any          Any          Any          Any          Any          Any          Any         
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Australia          0.000990747  0.00102238   0.00104535   0.00106707   0.0011966    0.00134927   0.00137213   0.00153062   0.00170403   0.00190516   0.00202583
   2 │ Austria            0.000103987  0.000103968  9.0544e-5    8.53889e-5   0.00010007   0.000103614  0.000109087  0.000107653  0.000111959  0.000118265  0.000133616
   3 │ Belgium            0.000548718  0.000512234  0.000458746  0.000381597  0.000410967  0.000426637  0.000444979  0.000451222  0.000465334  0.000512454  0.000569875
   4 │ Canada             0.00832702   0.00790575   0.00662924   0.00675101   0.00772083   0.00882244   0.00891513   0.00979034   0.0118275    0.0126891    0.0115645
   5 │ Denmark            8.21952e-5   6.84335e-5   6.17971e-5   5.65799e-5   6.41921e-5   8.4313e-5    8.54746e-5   8.76725e-5   7.98611e-5   7.8859e-5    7.46455e-5
   6 │ Finland            8.12688e-5   7.86115e-5   6.68769e-5   6.49423e-5   6.62465e-5   7.61118e-5   7.96103e-5   9.31865e-5   9.68827e-5   0.00011073   0.00011347
   7 │ France             0.00754285   0.00644225   0.00596744   0.00580938   0.00561904   0.00592596   0.00565579   0.0057622    0.00555473   0.00638032   0.00672812
   8 │ Germany            0.0154787    0.0137312    0.0126993    0.0129362    0.0134858    0.0139653    0.0151107    0.0149917    0.0153983    0.0142087    0.0155073
   9 │ Greece             4.32324e-5   4.27438e-5   4.09761e-5   4.48322e-5   4.8327e-5    5.60026e-5   5.97699e-5   6.86058e-5   8.28546e-5   8.98524e-5   0.000101018
  10 │ Ireland            8.83073e-6   8.32965e-6   7.71712e-6   7.35993e-6   7.1506e-6    8.16134e-6   9.54632e-6   9.53161e-6   1.03376e-5   1.22625e-5   1.50621e-5
  11 │ Italy              0.00413459   0.00411261   0.00388937   0.00372799   0.00393073   0.00423347   0.00436165   0.00508875   0.00529599   0.00587679   0.00595866
  12 │ Japan              0.0157031    0.0156941    0.014342     0.0138899    0.0167466    0.0173988    0.0171764    0.0184454    0.0221081    0.0245773    0.0273242
  13 │ Norway             7.61648e-5   6.70166e-5   5.88881e-5   6.68082e-5   7.83479e-5   9.38574e-5   0.000103544  0.000120039  9.98515e-5   0.000107198  0.000106517
  14 │ Spain              0.00101654   0.00103208   0.00115287   0.00114513   0.00109003   0.00106861   0.00112689   0.00119822   0.00137574   0.00168567   0.00156616
  15 │ Sweden             0.000302108  0.000255463  0.000209698  0.000203095  0.000207004  0.000241194  0.000241541  0.000256412  0.000254944  0.000233936  0.000221067
  16 │ U.K.               0.00809758   0.00716331   0.00660083   0.00705684   0.00729171   0.0078856    0.00842427   0.00907481   0.00949474   0.0102033    0.0103832
  17 │ the United States  0.75476      0.725291     0.628489     0.638697     0.747291     0.759032     0.774467     0.844417     0.951598     1.0          1.01465
```

### Vanilla DEA Example

`BDisposal` provides also a simple function for "vanilla" DEA computation.
In this example we have 4 DMU with two inputs and one output, and we compute the
efficiency of the latest DMU:

```julia
I  = [10   2;
       8   4;
      12 1.5;
      24   3]
O =  [100;80;120;120]
I₀ = [24 3]
O₀ = [120]

results = dmuEfficiency(I₀,O₀,I,O)
```
```
(computed = true, eff = false, obj = 0.5, wI = [0.041666666666666664, 0.0], wO = [0.004166666666666667], refSet = Dict(3 => 1.0), othDMUsEffConstrDuals = [0.0, 0.0, 1.0, 0.0], iRegConstrDual = 0.5)
```
Here we see that the last DMU is not efficient (indeed it's coefficient is only 0.5) as it is "dominated" by the third DMU that is the only efficient DMU in this set.


## Other packages

This is a list of other Julia packages that use classical DEA method, although without relaxing the disposability assumption.

- [DataEnvelopmentAnalysis.jl](https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl)
- [FrontierEfficiencyAnalysis.jl](https://github.com/wen-chih/FrontierEfficiencyAnalysis.jl)
- [SearchRef.jl](https://github.com/wen-chih/SearchRef.jl)
- [JuMP4DEA.jl](https://github.com/henry8527/JuMP4DEA.jl)


## References

- Abad, A. (2020) [Environmental Efficiency and Productivity Analysis](https://hal.inrae.fr/hal-03032038), _HAL_, 03032038
- Abad, A., P. Ravelojaona (2020a) [Pollution-adjusted Productivity Analysis: The Use of Malmquist and Luenberger Productivity Measures](https://doi.org/10.1002/mde.3260), _Managerial and Decision Economics_
- Abad, A., P. Ravelojoana (2020b) [A Generalization of Environmental Productivity Analysis](https://hal.inrae.fr/hal-02964799), _HAL_, 02964799
- Abad, A., W. Briec (2019) [On the Axiomatic of Pollution-generating Technologies: a Non-Parametric Approach](https://doi.org/10.1016/j.ejor.2019.02.027), _European Journal of Operational Research_, **277**(1), 377-390
- Abad, A. (2018) [Les Enseignements de la Micro-économie de la Production face aux Enjeux Environnementaux: Etude des Productions Jointes. Théorie et Applications](https://tel.archives-ouvertes.fr/tel-01963415), _Ph.D dissertation_, University of Perpignan.
- Abad, A., P. Ravelojaona (2017) [Exponential environmental productivity index and indicators](https://doi.org/10.1007/s11123-017-0513-7), _Journal of Productivity Analysis_, **48**(2), 147-166.
- Abad, A. (2015) [An environmental generalised Luenberger-Hicks-Moorsteen productivity indicator and an environmental generalised Hicks-Moorsteen productivity index](https://doi.org/10.1016/j.jenvman.2015.06.055), _Journal of Environmental Management_, **161**, 325-334.






## Acknowledgements

The development of this package at the _Bureau d'Economie Théorique et Appliquée_ (BETA, Nancy) was supported by the French National Research Agency through the [Laboratory of Excellence ARBRE](http://mycor.nancy.inra.fr/ARBRE/), a part of the “Investissements d'Avenir” Program (ANR 11 – LABX-0002-01).

[![BLogos](assets/logos_betaumr.png)](http://www.beta-umr7522.fr/)
