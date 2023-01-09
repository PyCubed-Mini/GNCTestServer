module Simulator
using SatelliteDynamics
import Base: *, +

export AbstractState, AbstractControl, r4k

*(v::NamedTuple, s::Number) = map(y -> s * y, v)
*(s::Number, v::NamedTuple) = map(y -> s * y, v)
+(a::NamedTuple, b::NamedTuple) = map(+, a, b)

function rk4(state, control, t::Float64, dt::Float64, dynamics, normalize)
    """
    rk4(state, control, t, dt, normalize, dynamics)

    Rung-Kutta 4th order integrator, with the state normalized at the end.
    """
    k₁ = dt * dynamics(state, control, t)
    k₂ = dt * dynamics(state + k₁ * 0.5, control, t + dt * 0.5)
    k₃ = dt * dynamics(state + k₂ * 0.5, control, t + dt * 0.5)
    k₄ = dt * dynamics(state + k₃, control, t + dt)

    x⁺ = state + (1 / 6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄)

    return normalize(x⁺)
end

end