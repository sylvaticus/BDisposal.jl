# ![BDisposalLogos](assets/BDisposal_logo_30x30.png) BDisposal.jl Documentation

Welcome to the documentation for the [_Non parametric productivity analysis under the B-disposal assumption_](https://github.com/sylvaticus/BDisposal.jl) package.

The BDisposal package proposes a serie of environmental efficiency and productivity algorithms for non-parametric modelling when we relax the disposability assumption of some of the outputs and/or inputs (e.g. pollution). These efficiency and productivity measures are implemented through convex and non convex [Data Envelopment Analysis (DEA)](https://en.wikipedia.org/wiki/Data_envelopment_analysis) (aka _Frontier Efficiency Analysis_) models.


## Installation

Until the BDisposal package is registered, please use:
* `] add git@github.com/sylvaticus/BDisposal.jl.git`
Once the BDisposal package will be included in the standard Julia register, install it with:
* `] add BDisposal.jl`

## Loading the module(s)

```julia
using BDisposal
```

## Usage

The BDisposal package contains two prominent components.
First, the BDisposal package defines environmental efficiency indicators for a set of Decision Making Units (`efficiencyScores()`), computing the distance to the best environmental production procedures, i.e., to the efficient production frontier.
Second, the BDisposal package displays environmental productivity indices (`prodIndex()`). These productivity measures are implemented for different time periods (eg., years, months etc.) or spatial units (eg., countries, cities etc.), based on the aforementioned environmental efficiency indicators.

Both components can consider constant or variable returns to scale, convex or non-convex production frontiers and additive or multiplicative distance functions, and work with multiple "good" ("desirable") and "bad" ("undesirable") inputs and outputs (with undesirable inputs being optional).

They are detailed in their respective pages:

- [**`efficiencyScores()`**](efficiencyScores.html): Compute efficiency indicators and convexity test results;
- [**`prodIndex()`**](prodIndex.html): Compute productivity indexes;

## Examples

The following example (with data from _Abad and Briec, 2019_) show how to compute the efficiency indicators and productivity indexes for 14 main French airports, where the (good) inputs considered are _employees_ and _totalCosts_ and the (single in this case) bad and good outputs are respectively (airport) _co2emissions_ and _passengers_.


```julia
using DataFrames, CSV, BDisposal
println("Estimating French airports efficiency and productivity indexes...")

# Loading airport data and preparing the data..
# Data is stored as a  CSV files with 6 columns: period, dmu, co2emissions, passengers, employees, totalCosts
airportData = CSV.read(joinpath(abspath("BDisposal"),"test","data","airports.csv"),DataFrame; delim=';',copycols=true)
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
    gI[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodInputs])
    if nBI > 0
         bI[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadInputs])
    end
    gO[:,:,p] = convert(Matrix{Float64},periodData[:,airportGoodOutputs])
    bO[:,:,p] = convert(Matrix{Float64},periodData[:,airportBadOutputs])
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
```julia
prodIndexes = prodIndex(gI,gO,bO;
                   retToScale="constant",prodStructure="multiplicative",convexAssumption=true,
                   startθ=0,startμ=0,startλ=1.1)
```
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
