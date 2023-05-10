using Test
using LinearAlgebra
using Plots
using Random
using SatelliteDynamics
include("../src/SatellitePlayground.jl")
SP = SatellitePlayground

@testset "basic single variable integration" begin
    idx = x -> x
    @test SP.rk4(0.0, 0.0, 0.1, (x, t) -> 2 * x + 1) ≈ 0.11 atol = 0.01
    @test SP.rk4(0.0, 3.0, 0.1, (x, t) -> sin(x)) ≈ 0.009142653672834007 atol = 0.01
end;

function no_control(measurement)
    return zero(SP.Control)
end

@testset "Distance from center of earth" begin
    earth_radius = 6378.1 * 1000
    function log_step(hist, state)
        point = norm(state.position) - earth_radius
        push!(hist, point)
    end
    (data, time) = SP.simulate(no_control, max_iterations=10000, log_step=log_step, dt=10.0)
    data = SP.vec_to_mat(data)
    data /= 1000
    display(plot(time, data, title="Distance from earth", xlabel="Time (s)", ylabel="Distance from earth (km)", labels="r"))
end

function energy(state; J=[0.3 0 0; 0 0.3 0; 0 0 0.3], μ=3.9860044188e14)
    r = state.position
    ω = state.angular_velocity
    v = state.velocity

    KE = 0.5 * norm(v)^2
    PE = -μ / norm(r)
    Eᵣ = 0.5 * ω' * J * ω
    return Eᵣ + KE + PE
end

@testset "Conservation of energy" begin
    function log_step(hist, state)
        point = energy(state)
        push!(hist, point)
    end
    (data, time) = SP.simulate(no_control, max_iterations=100000, log_step=log_step, dt=0.5)
    data = SP.vec_to_mat(data)
    display(plot(time, data, title="Energy", xlabel="Time (s)", ylabel="Energy", labels="E"))
end

@testset "detumbling" begin
    Random.seed!(1234)
    function control_law(measurement)
        (ω, b) = measurement

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return SP.Control(
            m
        )
    end

    @inline function measure(state, env)
        return (state.angular_velocity, env.b)
    end

    @time (data, time) = SP.simulate(control_law, max_iterations=10_000, measure=measure, log_step=SP.angular_log_step)
    # 4.344277 seconds (15.07 M allocations: 795.512 MiB, 10.87% gc time, 2.39% compilation time)
    # 4.517416 seconds (14.91 M allocations: 780.852 MiB, 7.04% gc time, 2.35% compilation time)
    # 3.893083 seconds (15.89 M allocations: 776.579 MiB, 8.71% gc time, 1.71% compilation time)
    # 3.608502 seconds (13.51 M allocations: 711.263 MiB, 8.77% gc time, 1.50% compilation time)
    data = SP.vec_to_mat(data)
    display(plot(time, data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end

@testset "io" begin
    Random.seed!(1234)
    (data, time) = SP.simulate(`sh runfakesat.sh`, max_iterations=2000, log_step=SP.angular_log_step)
    data = SP.vec_to_mat(data)
    display(plot(time, data, title="Socket DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end

@testset "Simulate multiple detumbling satellites" begin
    Random.seed!(1235)
    function control_law(measurement)
        (ω, b) = measurement

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return SP.Control(
            m
        )
    end

    @inline function measure(state, env)
        return (state.angular_velocity, env.b)
    end

    function log_step(hist, states)
        ω = [norm(state.angular_velocity) for state in states]
        push!(hist, ω)
    end

    function terminate(states, env, time, i)
        return all(norm(state.angular_velocity) < 0.1 for state in states)
    end

    N = 4
    controls = [control_law for _ in 1:N]
    measures = [measure for _ in 1:N]

    @time (data, time) = SP.simulate_multiple(controls, max_iterations=50_000,
        measures=measures, log_step=log_step, terminal_condition=terminate, dt=0.5)
    data = SP.vec_to_mat(data)
    display(plot(time, data, title="DeTumbling Multiple Satellites", xlabel="Time (s)", ylabel="|Angular Velocity| (rad/s)", labels=["ω1" "ω2" "ω3" "ω4"]))
end
