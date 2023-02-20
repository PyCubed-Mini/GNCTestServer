using Test
using LinearAlgebra
using Plots
using Random
include("src/GNCTestServer.jl")

@testset "basic single variable integration" begin
    idx = x -> x
    @test GNCTestServer.rk4(0.0, 0.0, 0.1, (x, t) -> 2 * x + 1) ≈ 0.11 atol = 0.01
    @test GNCTestServer.rk4(0.0, 3.0, 0.1, (x, t) -> sin(x)) ≈ 0.009142653672834007 atol = 0.01
end;

@testset "detumbling" begin
    Random.seed!(1234)
    function control_law(state, params, t)
        ω = state.ω
        b = params.b

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return GNCTestServer.Control(
            m
        )
    end


    (data, time) = GNCTestServer.simulate(control_law, max_iterations=10000)
    display(plot(time, data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))


end

@testset "io" begin
    Random.seed!(1234)
    (data, time) = GNCTestServer.simulate(`sh runfakesat.sh`)
    display(plot(time, data, title="Socket DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end