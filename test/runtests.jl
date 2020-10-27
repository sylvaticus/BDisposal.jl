
using Test, DataFrames, CSV, BDisposal
println("Testing BDisposal...")

airportInputs      = ["employees","totalCosts"]
airportGoodOutputs = ["passengers"]
airportBadOutputs  = ["co2emissions"]
airportData        = CSV.read(joinpath(@__DIR__,"data","airports.csv"),DataFrame; delim=';',copycols=true)

# Call the function to get the efficiency measurements for constant returns to scale
(λ_crs, λ_convex_crs, λ_nonconvex_crs, nonConvTest_value_crs, nonConvTest_crs) = efficiencyScores(
airportInputs,airportGoodOutputs,airportBadOutputs,airportData,retToScale="constant", formattedOutput=false,
dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="multiplicative")
@test nonConvTest_value_crs[3,2]  ≈ 1.1587129278170434

(λ_vrs, λ_convex_vrs, λ_nonconvex_vrs, nonConvTest_value_vrs, nonConvTest_vrs) = efficiencyScores(
airportInputs,airportGoodOutputs,airportBadOutputs,airportData,retToScale="variable",formattedOutput=false,
dirGI=0,dirBI=0,dirGO=1,dirBO=-1, prodStructure="additive")
@test nonConvTest_value_vrs[3,3]  ≈ 7.432043216538459
