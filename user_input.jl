# where the user will input their own functions to be used
# and where some basic defaults will be

using GNCTestServer

struct trueState
    state::Any
    params::Any
end

struct PycubedMiniState
    position :: Vector{Float64}
    velocity :: Vector{Float64}
end

struct PycubedMini <: trueState
    state::PycubedMiniState
    params::Any
ends

function PyM_measurement(x::trueState)
    return x.state.position
end

function PyM_control_law(measurement, t)
    return Dict(position=>measurement.state.position)
end

GNCTestServer.simulate(PyM_control_law, PyM_measurement, )

