module Thermochem

using DataFrames, CSV

path = "C:/Users/jonas/OneDrive/Desktop/Programming/Julia/Chemical_Engineering/Packages/Data/Component_Data.csv"
gen_dat = DataFrame(CSV.File(path))
path = "C:/Users/jonas/OneDrive/Desktop/Programming/Julia/Chemical_Engineering/Packages/Data/Cp_Constants.csv"
Cp_dat = DataFrame(CSV.File(path))
path = "C:/Users/jonas/OneDrive/Desktop/Programming/Julia/Chemical_Engineering/Packages/Data/Fuller_Params.csv"
fuller_df = DataFrame(CSV.File(path));

# Write your package code here.
include("Transport_Coeffs.jl")
include("Miscellaneous.jl")

# Export publicly available functions

end

