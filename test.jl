using Test
include("lib/simulator.jl")

@testset "basic single variable integration" begin
    idx = x -> x
    @test Simulator.rk4(0.0, 0.0, 0.0, 0.1, (x, c, t) -> 2 * x + 1, idx) ≈ 0.11 atol = 0.01
    @test Simulator.rk4(0.0, 0.0, 3.0, 3.1, (x, c, t) -> sin(x), idx) ≈ 0.009142653672834007 atol = 0.01
end;