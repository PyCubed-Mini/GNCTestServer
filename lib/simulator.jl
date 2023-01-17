module Simulator
using SatelliteDynamics
using LinearAlgebra
import Base: *, +
include("quaternions.jl")
include("mag_field.jl")

export sim, Control, Parameters

*(v::NamedTuple, s::Number) = map(y -> s * y, v)
*(s::Number, v::NamedTuple) = map(y -> s * y, v)
+(a::NamedTuple, b::NamedTuple) = map(+, a, b)

mutable struct Parameters
    J::Array{Float64,2} # Inertia matrix
    b::Array{Float64,1} # Magnetic field
end

mutable struct Control
    m::Array{Float64,1} # Control input
end

"""
    rk4(x, derivative, t, dt)

    Rung-Kutta 4th order integrator
"""
function rk4(x, t, dt, derivative)
    k₁ = dt * derivative(x, t)
    k₂ = dt * derivative(x + k₁ * 0.5, t + dt * 0.5)
    k₃ = dt * derivative(x + k₂ * 0.5, t + dt * 0.5)
    k₄ = dt * derivative(x + k₃, t + dt)

    x⁺ = x + (1 / 6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄)

    return x⁺
end

function dynamics(state, parameters, control, t)
    u = cross(control.m, parameters.b)
    return (
        r=state.v,
        v=accel_perturbations(t, state.r, state.v),
        q=qdot(state.q, state.ω),
        ω=parameters.J \ (u - cross(state.ω, parameters.J * state.ω))
    )
end

function integrate_state(state, parameters, control, t::Epoch, dt::Float64)
    """
    rk4(state, control, t, dt, normalize, dynamics)

    Rung-Kutta 4th order integrator, with the state normalized at the end.
    """
    state⁺ = rk4(state, t, dt, (state, t) -> dynamics(state, parameters, control, t))

    return (
        r=state⁺.r,
        v=state⁺.v,
        q=state⁺.q / norm(state⁺.q), # normalize the quaternion
        ω=state⁺.ω
    )
end

"""
    accel_perturbations(epc, r, v; mass, area_drag, coef_drag, area_srp, coef_srp, n_grav, m_grav, third_body)

    Generates the acceleration for a spacecraft in LEO. Accounts for a variety of factors, 
    including spherical harmonics, atmospheric drag, SRP, and thirdbody from sun and moon

    ForwardDiff friendly (written by Kevin)
"""
function accel_perturbations(epc::Epoch, r, v;
    mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3, area_srp::Real=1.0,
    coef_srp::Real=1.8, n_grav::Integer=10, m_grav::Integer=10, third_body::Bool=true)

    # These functions don't like static, so we gotta adjust
    x = [r; v]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E = earth_rotation(epc)
    W = polar_motion(epc)
    R = W * E * PN

    # Compute sun and moon position
    r_sun = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration (eltype(x) makes this forward diff friendly)
    a = zeros(eltype(x), 3)

    # spherical harmonic gravity
    a += accel_gravity(x, R, n_grav, m_grav)

    # atmospheric drag
    ρ = density_harris_priester(epc, r)
    a += accel_drag([r; v], ρ, mass, area_drag, coef_drag, Array{Real,2}(PN))

    # SRP
    nu = eclipse_cylindrical(x, r_sun)  # conical doesn't work correctly
    a += nu * accel_srp(x, r_sun, mass, area_srp, coef_srp)

    if third_body
        # third body sun
        a += accel_thirdbody_sun(x, r_sun)

        # third body moon
        a += accel_thirdbody_moon(x, r_moon)
    end

    return (a)
end

""" initialize_orbit(; r, v, a, q, ω, Rₑ, μ)

Initializes a random, viable orbit given a few different terms, usually 
a position 'r' in Cartesian coordinates. Initial velocity may be specified, but 
if specified it will not necessarily result in a stable orbit. 

The initial starting position, velocity, semi-major axis, orientation, and angular 
velocity may be either specified or determined randomly. 

Arguments:
- r:  (Optional) Height above ground that the satellite will start its orbit at    |  Scalar 
- v:  (Optional) Magnitude of initial velocity                                     |  Scalar 
- a:  (Optional) Semi-major axis                                                   |  Scalar 
- q:  (Optional) Initial attitude, as a unit quaternion                            |  [4,]
- ω:  (Optional) Initial angular velocity                                          |  [3,]

Returns:
- x:  Initial state, as 
          x = [r, v, q, ω]
"""
function initialize_orbit(; r=nothing, v=nothing, a=nothing, q=nothing, ω=nothing, Rₑ=6378.1363e3, μ=3.9860044188e14)

    ### POSITION ###

    # If unspecified, generate a random altitude for the satellite 
    _r = !isnothing(r) ? r : (rand(100:1000) * 1000) + Rₑ # Random height, in km 

    # If unspecified, pick some amount of eccentricity 
    _a = !isnothing(a) ? a * _r : rand(1.0:0.05:1.5) * _r   # Semi-major axis. Results in a circular orbit for a == r

    # If unspecified, calculate the necessary initial velocity to maintain a valid orbit 
    _v = !isnothing(v) ? v : sqrt(μ * ((2 / _r) - (1 / _a)))

    # Now we have to pick directions for r and v... 
    axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    ax = rand(1:3)
    r₀ = _r * splice!(axes, ax)  # removes the selected axis from axes 
    v₀ = _v * rand(axes)         # Remaining options are both perpendicular and valid 

    #   allow for negative axes...
    r₀ = (randn() > 0.0) ? r₀ : -r₀
    v₀ = (randn() > 0.0) ? v₀ : -v₀


    ### ORIENTATION ###

    # If unspecified, generate a random starting orientation
    q₀ = !isnothing(q) ? q : randn(4)
    q₀ /= norm(q₀)  # Ensure it is unit 

    # If unspecified, generate a random starting angular velocity
    ω₀ = !isnothing(ω) ? ω : 0.5 * randn(3)

    return (r=r₀, v=v₀, q=q₀, ω=ω₀)
end

function control_step(state, params::Parameters, control_fn, t)
end

""" sim_step (state, params, control_fn, t, dt)

steps the simulation forward by one time step, dt.
given the current state, free parameters, control function, and time.

Arguments:
- state:      Current state of the system, as a State struct                 |  State
- params:     Free parameters of the system, as a Parameters struct          |  Parameters
- control_fn: Function of state, parameters, and time that yields a control  |  Function
- t:          Current time, as an Epoch struct                               |  Epoch
- dt:         Time step, in seconds                                          |  Scalar
"""
function sim_step(state, params::Parameters, control_fn, t, dt)
    params.b = IGRF13(state.r, t)
    control = control_fn(state, params, t)
    state = integrate_state(state, params, control, t, dt)
    t += dt                      # Don't forget to update time (not that it really matters...)
    return (state, t)
end

function sim(control_fn, log_init=default_log_init, log_step=default_log_step)
    x = initialize_orbit()
    println("intialized orbit!")

    J = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
    params = Parameters(J, [0.0, 0.0, 0.0])

    t = Epoch(2020, 11, 30)          # Starting time is Nov 30, 2020
    dt = 0.5                         # Time step, in seconds

    N = 100000
    hist = log_init(x, N)
    for i = 1:N-1
        (x, t) = sim_step(x, params, control_fn, t, dt)
        log_step(hist, x, i)
        if norm(x.ω) < 0.001
            hist = hist[1:i, :]
            break
        end
    end

    return hist
end

function default_log_init(state, iterations)
    hist = zeros(iterations, 4)
    hist[1, 1:3] .= state.ω
    hist[1, 4] = norm(state.ω)
    return hist
end

function default_log_step(hist, state, i)
    hist[i+1, 1:3] .= state.ω
    hist[i+1, 4] = norm(state.ω)
end

end