using Test
using LinearAlgebra
using Plots
using Random
include("../src/GNCTestServer.jl")

@testset "basic single variable integration" begin
    idx = x -> x
    @test GNCTestServer.rk4(0.0, 0.0, 0.1, (x, t) -> 2 * x + 1) ≈ 0.11 atol = 0.01
    @test GNCTestServer.rk4(0.0, 3.0, 0.1, (x, t) -> sin(x)) ≈ 0.009142653672834007 atol = 0.01
end;

@testset "detumbling" begin
    Random.seed!(1234)
    function control_law(measurement, t)
        ω = measurement[1].angular_velocity
        b = measurement[2].b

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return GNCTestServer.Control(
            m
        )
    end

    function measurement_id(state, params)
        return (state, params)
    end
    struct basicSat
        state
        params
    end
    initial_state = basicSat(GNCTestServer.initialize_orbit(), GNCTestServer.initialize_params())
    (data, time) = GNCTestServer.simulate(control_law, measurement_id, initial_state, max_iterations=10000)
    display(plot(time, data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end

@testset "io" begin
    Random.seed!(1234)
    function measurement(state, params)
        return (state, params)
    end
    initial_state = basicSat(GNCTestServer.initialize_orbit(), GNCTestServer.initialize_params())
    (data, time) = GNCTestServer.simulate(`sh runfakesat.sh`, measurement, initial_state, max_iterations=2000)
    display(plot(time, data, title="Socket DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end