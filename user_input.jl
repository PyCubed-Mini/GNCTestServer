# where the user will input their own functions to be used
# and where some basic defaults will be

using GNCTestServer
using LinearAlgebra
using Plots
using Random

struct Py4State
    position::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
end

struct Py4
    state::Py4State
    params
end

struct basicSat
    state::RBState
    params::GNCTestServer.Parameters
end
# struct PycubedMiniState


function PyM_measurement(state, params)
    return state.position
end

function PyM_control_law(measurement, t)
    return GNCTestServer.Control(
        [0, 0, 0]
    )
end

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


initial_state = basicSat(GNCTestServer.initialize_orbit(), GNCTestServer.initialize_params())

(data, time) = GNCTestServer.simulate(control_law, measurement_id, initial_state)
display(plot(time, data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
