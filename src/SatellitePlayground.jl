# Version 0.1.0

module SatellitePlayground
import Base: *, +
using LinearAlgebra
using MsgPack
using Printf
using SatelliteDynamics
using Sockets
include("quaternions.jl")
include("mag_field.jl")
include("communication.jl")
include("RBState.jl")
include("initial_conditions.jl")

export simulator, Control, Parameters, RBState, world_to_body

"""
    Parameters

"""
mutable struct Parameters
    J::Array{Float64,2} # Inertia matrix
end

Base.copy(p::Parameters) = Parameters(copy(p.J))

mutable struct Control
    m::Array{Float64,1} # Control input
end

function Base.zero(::Type{Control})
    return Control(zeros(3))
end

J = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
default_parameters = Parameters(J)

"""
    Default enviorment type for the simulation

Should be in the body frame
"""
mutable struct DefaultEnvironment
    time::Epoch
    b::Array{Float64,1} # Magnetic field
end
Base.copy(e::DefaultEnvironment) = DefaultEnvironment(e.time, copy(e.b))

default_environment = DefaultEnvironment(
    Epoch(2020, 11, 30),
    [0.0, 0.0, 0.0]
)

"""
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

""" 
Returns the derivative of the state, given the current state, parameters, control input, and time.

Arguments:
- state: the current state of the spacecraft                             | State
- parameters: the parameters of the spacecraft                           | Parameters
- control: the control input to apply to the spacecraft                  | Control
- t: the current time                                                    | Epoch

Returns:
- state: the derivative of the state                                     | State
"""
function dynamics(state::RBState, parameters::Parameters, environment::DefaultEnvironment, control::Control, t::Epoch)::RBState
    u = cross(control.m, environment.b)
    return RBState(
        state.velocity,
        accel_perturbations(t, state.position, state.velocity),
        qdot(state.attitude, state.angular_velocity),
        parameters.J \ (u - cross(state.angular_velocity, parameters.J * state.angular_velocity))
    )
end

""" 
Integrate the state forward by dt seconds, using the given control input.

Arguments
- state: the current state of the spacecraft                             | State
- parameters: the parameters of the spacecraft                           | Parameters
- control: the control input to apply to the spacecraft                  | Control
- t: the current time                                                    | Epoch
- dt: the time step to integrate forward by                              | Float64

Returns
- state: the state of the spacecraft after the integration               | State
"""
function integrate_state(state::RBState, parameters::Parameters, environment::DefaultEnvironment, control::Control, dt::Float64)::RBState
    function time_dynamics(state::RBState, t)
        return dynamics(state, parameters, environment, control, t)
    end
    new_state = rk4(state, environment.time, dt, time_dynamics)
    return renorm(new_state)
end

"""
Generates the acceleration for a spacecraft in LEO. Accounts for a variety of factors, 
including spherical harmonics, atmospheric drag, SRP, and thirdbody from sun and moon

ForwardDiff friendly (written by Kevin)
"""
function accel_perturbations(epc::Epoch, r, v;
    mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3, area_srp::Real=1.0,
    coef_srp::Real=1.8, n_grav::Integer=10, m_grav::Integer=10, third_body::Bool=true)

    # These functions don't like static, so we gotta adjust
    r = Vector(r)
    v = Vector(v)
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

    return a
end

"""
    initialize_orbit(; r=nothing, v=nothing, a=nothing, q=nothing, ω=nothing, Rₑ=6378.137e3, μ=3.9860044188e14)

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
- x:  Initial state, as (r, v, q, ω)                                               |  State
"""
function initialize_orbit(; r=nothing, v=nothing, a=nothing, q=nothing, ω=nothing, Rₑ=6378.1363e3, μ=3.9860044188e14)::RBState

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

    return RBState(r₀, v₀, q₀, ω₀)
end

function update_environment(state::RBState, env::DefaultEnvironment)
    env.b = world_to_body(state, IGRF13(state.position, env.time))
end

mutable struct FunctionSim
    dt::Float64
    control::Control
end

"""
    simulate(control::Function)
    simulate(launch::Cmd)
    simulate(control::Function, log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000,
    dt=0.5, initial_condition=nothing, measure=default_measure)

Runs a simulation from a random initial condition (or from `initial_condition`) if given.
The simulation runs for `max_iterations` steps.

The control input to the magnetorquer coils at each time step is set either by the given control 
function, or by the GNCTestClient launched by the launch command.

By default the controller recieves the (state, parameters, time) in the body frame.
But this can be changed by setting the measurement function `measure`.

The simulation logs the angular velocity and its mangitude by default.
However, by setting the log_* functions one can log arbitrary data.
"""
function simulate(control::Function; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000, dt=0.5,
    initial_condition=nothing, measure=default_measure, initial_parameters=default_parameters,
    initial_environment=default_environment, silent=false)
    function setup()
        return FunctionSim(dt, Control([0.0, 0.0, 0.0]))
    end
    function step(sim, measurement, i)
        sim.control = control(measurement)
    end
    function cleanup(sim)
        return
    end
    return simulate_helper(setup, step, cleanup,
        log_init, log_step,
        terminal_condition, max_iterations,
        initial_condition, measure,
        initial_parameters, initial_environment,
        silent)
end

mutable struct SocketSim
    dt::Float64
    satellite_process::Base.Process
    uplink::Array{UInt8}
    uplink_ptr::Any # Kept to prevent garbage collection
    uplink_sem::Any
    downlink::Array{UInt8}
    downlink_ptr::Any # Kept to prevent garbage collection
    downlink_sem::Any
    control::Control
end

function simulate(launch::Cmd; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000,
    initial_condition=nothing, measure=default_measure,
    initial_parameters=default_parameters, initial_environment=default_environment,
    silent=false)
    function setup()
        println("Creating shared memory and semaphores...")
        uplink, uplink_ptr = mk_shared("gnc_uplink", 128)
        uplink_sem = mk_semaphore(67)
        downlink, downlink_ptr = mk_shared("gnc_downlink", 128)
        downlink_sem = mk_semaphore(68)
        println("Launching satellite...")
        sleep(0.1)
        satellite_process = run(launch, wait=false)
        println("Launched")
        return SocketSim(
            0.5,
            satellite_process,
            uplink,
            uplink_ptr,
            uplink_sem,
            downlink,
            downlink_ptr,
            downlink_sem,
            Control([0.0, 0.0, 0.0])
        )
    end
    function step(sim, measurement, i)
        if sim.satellite_process.exitcode >= 0
            throw(error("Satellite process exited with code $(sim.satellite_process.exitcode), check /tmp/satlog.txt for details"))
        end
        downlink(sim.downlink, sim.downlink_sem, measurement)
        uplink_data = uplink(sim.uplink, sim.uplink_sem, i)
        sim.control = Control(uplink_data["m"])
        sim.dt = uplink_data["dt"]
    end
    function cleanup(sim)
        kill(sim.satellite_process)
        sim.uplink_sem.remove()
        sim.downlink_sem.remove()
        println("Killed satellite process")
    end
    return simulate_helper(setup, step, cleanup,
        log_init, log_step,
        terminal_condition, max_iterations,
        initial_condition, measure,
        initial_parameters, initial_environment,
        silent)
end

function print_iteration(i, max_iterations, state, env, sim)
    @printf("[%d/%d]: norm(ω)=%.3f r=<%.3f %.3f %.3f> b=<%.3f %.3f %.3f> dt=%.3f",
        i, max_iterations, norm(state.angular_velocity), state.position[1], state.position[2], state.position[3], env.b[1],
        env.b[2], env.b[3], sim.dt)
end

function print_iteration(i, max_iterations)
    @printf("[%d/%d]", i, max_iterations)
end

function simulate_multiple(controls; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000, dt=0.1,
    initial_conditions=nothing, measures=nothing,
    initial_parameters=nothing, initial_environment=default_environment)
    N = length(controls)
    @assert (N > 0)
    @assert (isnothing(initial_conditions) || length(initial_conditions) == N)
    @assert (isnothing(measures) || length(measures) == N)

    if isnothing(initial_conditions)
        states = [initialize_orbit() for _ = 1:N]
    else
        states = initial_conditions
    end
    if isnothing(initial_parameters)
        params = [copy(default_parameters) for _ = 1:N]
    else
        params = initial_parameters
    end
    if isnothing(measures)
        measure_funcs = [default_measure for _ = 1:N]
    else
        measure_funcs = measures
    end
    println("intialized orbit!")


    env = copy(initial_environment)
    start_time = env.time

    hist = log_init(state)
    time_hist = []

    for i = 1:max_iterations

        for i = 1:N
            update_environment(states[i], env)
            u = controls[i](measure_funcs[i](states[i], env))
            states[i] = integrate_state(states[i], params[i], env, u, dt)
        end
        env.time += dt

        log_step(hist, states)
        append!(time_hist, env.time - start_time)
        print("\r\033[K")
        print_iteration(i, max_iterations)
        if terminal_condition(states, env, time, i)
            break
        end
    end
    print("\n")
    println("Simulation complete!")

    return (hist, time_hist)
end

function simulate_helper(setup::Function, step::Function, cleanup::Function,
    log_init::Function, log_step::Function,
    terminal_condition::Function, max_iterations, initial_condition, measure::Function,
    initial_parameters::Parameters, initial_environment, silent::Bool)
    if isnothing(initial_condition)
        state = initialize_orbit()
    else
        state = initial_condition
    end
    if !silent
        println("intialized orbit!")
    end

    params = copy(initial_parameters)

    sim = setup()

    env = copy(initial_environment)
    start_time = env.time

    hist = log_init(state)
    time_hist = []

    try
        for i = 1:max_iterations
            update_environment(state, env)

            step(sim, measure(state, env), i)

            state = integrate_state(state, params, env, sim.control, sim.dt)
            env.time += sim.dt

            log_step(hist, state)
            append!(time_hist, env.time - start_time)
            print("\r\033[K")
            if !silent
                print_iteration(i, max_iterations, state, env, sim)
            end
            if terminal_condition(state, params, time, i)
                break
            end
        end
        if !silent
            print("\n")
            println("Simulation complete!")
        end
        cleanup(sim)

        return (hist, time_hist)
    catch e
        println("Simulation failed: $e")
        cleanup(sim)
        throw(e)
    end

end

"""
Default function to initialize the log, only logs angular velocity and magnitude.

Arguments:
- state:       Current state of the system, as a State struct            |  State
- iterations:  Number of iterations the simulation is run                |  Scalar

Returns:
- hist:        Initialized log of the simulation                         |  Matrix
"""
function default_log_init(state)
    return []
end

"""
Logs angular velocity and magnitude.

Arguments:
- hist:        Log of the simulation                                     |  Matrix
- state:       Current state of the system, as a State struct            |  State
"""
function angular_log_step(hist, state)
    point = [state.angular_velocity; norm(state.angular_velocity)]
    push!(hist, point)
end

function default_log_step(hist, state)
    push!(hist, state)
end

"""
    vec_to_mat(hist)

Converts a vector of vectors to a matrix.
Useful for converting the output of a simulation to a format that can be plotted.
"""
function vec_to_mat(hist)
    if hist[1] isa RBState || hist[1] isa SVector
        hist = Vector.(hist)
    end
    hist = reduce(hcat, hist)
    hist = hist'
    return hist
end

"""
    default_terminate(state, params, time, i)

default termination condition for `simulate`, always returns `false`.
"""
function default_terminate(state, params, time, i)
    return false
end

"""
    default_measure(state, params, t)

Default measurement function used by `simulate`, returns the state and parameters.
"""
@inline function default_measure(state, env)
    return (state, env)
end

"""
    world_to_body(::RBState, ::Vector)

Returns the vector in the body frame
"""
@inline function world_to_body(state, vec)
    return quaternionToMatrix(state.attitude)' * vec
end

end