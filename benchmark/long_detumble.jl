using Test
using LinearAlgebra
using SatelliteDynamics
using Plots
using GNCTestServer
using ProfileView

py4_dipole_limits = [0.06997731147540984,
    0.053130000000000004,
    0.06976756111111111]

function control_law(m, t)
    ω = m[1].angular_velocity
    b = m[2].b

    b̂ = b / norm(b)
    k = 7e-4
    M = -k * (I(3) - b̂ * b̂') * ω
    m = 1 / (dot(b, b)) * cross(b, M)
    return GNCTestServer.Control(
        clamp.(m, -py4_dipole_limits, py4_dipole_limits)
    )
end

x_osc_0 = [400e3 + SatelliteDynamics.R_EARTH, 0.0, deg2rad(50), deg2rad(-1.0), 0.0, 0.0] # a, e, i, Ω, ω, M
q0 = [1.0, 0.0, 0.0, 0.0]
ω0 = 0.1 * [0.3, 0.1, -0.2]
ω0 = ω0 / norm(ω0) * deg2rad(50.0)

x0 = GNCTestServer.state_from_osc(x_osc_0, q0, ω0)

function measurement(state)
    return [state.angular_velocity; norm(state.angular_velocity)]
end

mutable struct log_state
    cur_count::Int
    cur_measure::Vector{Float64}
    measures::Vector{Vector{Float64}}
end

function log_init(state)
    measure = measurement(state)
    return log_state(1, measure, [])
end

down_sample_rate = 100

function log_step(hist, state)
    hist.cur_count += 1
    hist.cur_measure += measurement(state)
    if hist.cur_count == down_sample_rate
        hist.cur_measure /= down_sample_rate
        push!(hist.measures, hist.cur_measure)
        hist.cur_count = 0
        hist.cur_measure = zeros(4)
    end
end

function log_end(hist)
    println(hist.measures)
    return GNCTestServer.default_log_end(hist.measures)
end

function terminal_condition(state, params, t, i)
    return norm(state.angular_velocity) < deg2rad(0.1)
end

day = 60 * 60 * 24
time_step = 0.1
@time (data, time) = GNCTestServer.simulate(control_law, max_iterations=day / time_step, dt=time_step,
    log_init=log_init, log_step=log_step, log_end=log_end, initial_condition=x0, terminal_condition=terminal_condition)
time = time[1:down_sample_rate:end]
time = time[1:size(data)[1]]
time /= 60

data = rad2deg.(data)


# display(plot(time, data, title="DeTumbling", xlabel="Time (minutes)", ylabel="Angular Velocity (deg/s)", labels=["ω1" "ω2" "ω3" "ω"]))

#  18.903280 seconds (99.61 M allocations: 4.767 GiB, 7.28% gc time, 0.90% compilation time: 77% of which was recompilation)