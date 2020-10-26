
using Test, DataFrames, CSV, BDisposal
println("Testing BDisposal...")

airportInputs      = ["employees","totalCosts"]
airportGoodOutputs = ["passengers"]
airportBadOutputs  = ["co2emissions"]
airportData        = CSV.read(joinpath(@__DIR__,"data","airports.csv"),DataFrame; delim=';',copycols=true)

# Call the function to get the efficiency measurements for constant returns to scale
(λ_crs, λ_convex_crs, λ_nonconvex_crs, nonConvTest_value_crs, nonConvTest_crs) = efficiencyScores(
airportInputs,airportGoodOutputs,airportBadOutputs,airportData,retToScale="constant", formattedOutput=false)

(λ_vrs, λ_convex_vrs, λ_nonconvex_vrs, nonConvTest_value_vrs, nonConvTest_vrs) = efficiencyScores(
airportInputs,airportGoodOutputs,airportBadOutputs,airportData,retToScale="variable",formattedOutput=false)

@test nonConvTest_value_crs[3,2]  ≈ 1.1084313044037957
@test nonConvTest_value_vrs[3,3]  ≈ 0.9614904634202474
