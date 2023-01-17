using Test
using LinearAlgebra
using Plots
include("lib/GNCTestServer.jl")

@testset "basic single variable integration" begin
    idx = x -> x
    @test Simulator.rk4(0.0, 0.0, 0.1, (x, t) -> 2 * x + 1) ≈ 0.11 atol = 0.01
    @test Simulator.rk4(0.0, 3.0, 0.1, (x, t) -> sin(x)) ≈ 0.009142653672834007 atol = 0.01
end;

@testset "detumbling" begin
    function control_law(state, params, t)
        ω = state.ω
        b = params.b

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return Simulator.Control(
            m
        )
    end


    data = Simulator.sim(control_law)
    display(plot(data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))


end