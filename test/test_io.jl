begin
    using LinearAlgebra
    include("deps.jl")
end

function default_detumble_controller(measurement)
    (ω, b) = measurement

    b̂ = b / norm(b)
    k = 7e-4
    M = -k * (I(3) - b̂ * b̂') * ω
    m = 1 / (dot(b, b)) * cross(b, M)
    return SP.Control(
        m
    )
end

@testset "io" begin
    SIMULATION_ITERATIONS = 100

    model = copy(SP.pqmini_model)
    model.control_limit = [Inf, Inf, Inf]
    @time (data, time) = SP.simulate(`python3 python_satellites/fakecubesat.py`, max_iterations=SIMULATION_ITERATIONS, log_step=SP.angular_log_step,
        model=model, initial_condition=default_data.state)
    ω_sil = data[end]

    @inline function measure(state, env)
        return (state.angular_velocity, env.b)
    end
    @time (data, time) = SP.simulate(default_detumble_controller, max_iterations=SIMULATION_ITERATIONS, log_step=SP.angular_log_step,
        model=model, initial_condition=default_data.state, measure=measure)
    ω_sim = data[end]

    @test ω_sim ≈ ω_sil
end

@testset "log_position" begin
    SIMULATION_ITERATIONS = 10

    model = copy(SP.pqmini_model)
    model.control_limit = [Inf, Inf, Inf]

    function log_step(hist, step)
        push!(hist, step.position)
    end

    @inline function measure(state, env)
        return Dict(
            :r => state.position,
        )
    end

    @time (data, time) = SP.simulate(`python3 python_satellites/log_position.py `, max_iterations=SIMULATION_ITERATIONS, log_step=log_step,
        model=model, initial_condition=default_data.state, measure=measure)
end