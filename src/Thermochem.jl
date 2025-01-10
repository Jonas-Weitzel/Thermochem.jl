module Thermochem

using DataFrames, CSV

# Write your package code here.
include("Transport_Coeffs.jl")

# Export publicly available functions
export calc_D12

end

