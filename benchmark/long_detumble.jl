using Test
using LinearAlgebra
using SatelliteDynamics
using Plots
using SatellitePlayground

SP = SatellitePlayground

function measure(state, env)
    return (state.angular_velocity, env.b)
end

function control_law(measurement)
    (ω, b) = measurement

    b̂ = b / norm(b)
    k = 7e-6
    M = -k * (I(3) - b̂ * b̂') * ω
    m = 1 / (dot(b, b)) * cross(b, M)
    return SP.Control(
        m
    )
end

x_osc_0 = [400e3 + SatelliteDynamics.R_EARTH, 0.0, deg2rad(50), deg2rad(-1.0), 0.0, 0.0] # a, e, i, Ω, ω, M
q0 = [1.0, 0.0, 0.0, 0.0]
ω0 = 0.1 * [0.3, 0.1, -0.2]
ω0 = ω0 / norm(ω0) * deg2rad(5.0)

x0 = SP.state_from_osc(x_osc_0, q0, ω0)

function log_step(hist, state)
    point = [state.angular_velocity; norm(state.angular_velocity)]
    push!(hist, point)
end

function terminal_condition(state, env, t, i)
    return norm(state.angular_velocity) < deg2rad(0.1)
end

day = 60 * 60 * 24
time_step = 0.1
@time (data, time) = SP.simulate(control_law, max_iterations=day / time_step, dt=time_step,
    log_step=log_step, initial_condition=x0, terminal_condition=terminal_condition,
    measure=measure)

down_sample_rate = 100

data = data[1:down_sample_rate:end]
data = SP.vec_to_mat(data)

time = time[1:down_sample_rate:end]
time = time[1:size(data)[1]]
time /= 60

data = rad2deg.(data)


display(plot(time, data, title="DeTumbling", xlabel="Time (minutes)", ylabel="Angular Velocity (deg/s)", labels=["ω1" "ω2" "ω3" "ω"]))
