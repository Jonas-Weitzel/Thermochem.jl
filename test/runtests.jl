using Thermochem
using Test

T = 175 + 273; p = 1e5
comp1 = "Ar"; comp2 = "H2"
D12 = calc_D12(T, p, comp1, comp2)

@testset "Thermochem.jl" begin

    # Write your tests here.
    # @test calc_D12(T, p, comp1, comp2; fuller_df=fuller_df) == 4          # 2^2 = 4

end
