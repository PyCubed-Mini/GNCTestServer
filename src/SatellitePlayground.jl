module SatellitePlayground
import Base: *, +
using LinearAlgebra
using Printf
using SatelliteDynamics
using Sockets
include("quaternions.jl")
include("mag_field.jl")
include("communication.jl")
include("RBState.jl")
include("initial_conditions.jl")
include("environment.jl")
include("model.jl")
include("dynamics.jl")

export simulator, Control, Model, RBState, world_to_body

function Base.zero(::Type{Control})
    return Control(zeros(3))
end

"""
Rung-Kutta 4th order integrator
"""
@inline function rk4(x, t, dt, derivative::F) where F
    k₁ = dt * derivative(x, t)
    k₂ = dt * derivative(x + k₁ * 0.5, t + dt * 0.5)
    k₃ = dt * derivative(x + k₂ * 0.5, t + dt * 0.5)
    k₄ = dt * derivative(x + k₃, t + dt)

    x⁺ = x + (1.0 / 6.0) * (k₁ + 2.0 * k₂ + 2.0 * k₃ + k₄)

    return x⁺
end

"""
Rung-Kutta 4th order integrator, with full simulator parameters
"""
@inline function rk4(x, dt, parameters, environment, control)
    t = environment.time
    k₁ = dt * dynamics(x, parameters, environment, control)
    environment.time = t + dt * 0.5
    k₂ = dt * dynamics(x + k₁ * 0.5, parameters, environment, control)
    k₃ = dt * dynamics(x + k₂ * 0.5, parameters, environment, control)
    environment.time = t + dt
    k₄ = dt * dynamics(x + k₃, parameters, environment, control)
    environment.time = t

    return x + (1.0 / 6.0) * (k₁ + 2.0 * k₂ + 2.0 * k₃ + k₄)
end


""" 
Integrate the state forward by dt seconds, using the given control input.

Arguments
- state: the current state of the spacecraft                             | State
- parameters: the parameters of the spacecraft                           | Model
- control: the control input to apply to the spacecraft                  | Control
- t: the current time                                                    | Epoch
- dt: the time step to integrate forward by                              | Float64

Returns
- state: the state of the spacecraft after the integration               | State
"""
function integrate_state(state::RBState, parameters::Model, environment::Environment, control::Control, dt::Float64)::RBState
    new_state = rk4(state, dt, parameters, environment, control)
    return renorm(new_state)
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

mutable struct FunctionSim
    dt::Float64
    control::Control
end

"""
    simulate(control::Function)
    simulate(launch::Cmd)
    simulate(control::Function; log_init=default_log_init, log_step=default_log_step,
        terminal_condition=default_terminate, max_iterations=1000, dt=0.5,
        initial_condition=nothing, measure=default_measure, model=default_model,
        environment=default_environment, silent=false)

Runs a simulation from a random initial condition (or from `initial_condition`) if given.
The simulation runs for `max_iterations` steps.

The model specifies things like the mass, inertia, drag coefficient of the satellite and control type (magnetorquer, reaction wheel, etc).
The environment specifies things like the magnetic field.

The control input is computed via `control(measurement)` where `measurement` is the output of the measurement function.
Or it is computed by an external process using `GNCTestClient.py` if `launch` is given.

By default the controller recieves the (state, environment) in the body frame.
But this can be changed by setting the measurement function `measure`.

The simulation logs the state of the satellite at each time step by default.
This can be changed by setting the log_* functions.
"""
function simulate(control::Function; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000, dt=0.5,
    initial_condition=nothing, measure=default_measure, model=default_model,
    environment=default_environment, silent=false)
    function setup()
        return FunctionSim(dt, Control([0.0, 0.0, 0.0]))
    end
    function step(sim, measurement, i)
        sim.control = control(measurement)
    end
    function cleanup(sim)
        return
    end
    if initial_condition !== nothing
        initial_condition = copy(initial_condition)
    end
    return simulate_helper(setup, step, cleanup,
        log_init, log_step,
        terminal_condition, max_iterations,
        initial_condition, measure,
        copy(model), copy(environment),
        silent)
end

mutable struct SocketSim
    dt::Float64
    satellite_process::Base.Process
    socket::Any
    control::Control
end

function simulate(launch::Cmd; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000, dt=0.5,
    initial_condition=nothing, measure=sil_measure,
    model=default_model, environment=default_environment,
    silent=false)
    function setup()
        println("Initializing ZMQ socket...")
        socket = init_zmq_socket("tcp://*:5555")
        println("Launching satellite...")
        satellite_process = run(launch, wait=false)
        println("Launched")
        return SocketSim(
            dt,
            satellite_process,
            socket,
            Control([0.0, 0.0, 0.0])
        )
    end
    function step(sim, measurement, i)
        if sim.satellite_process.exitcode >= 0
            throw(error("Satellite process exited with code $(sim.satellite_process.exitcode)"))
        end
        uplink(measurement, sim.socket)
        downlink_data = downlink(sim.socket)
        sim.control = Control(downlink_data["m"])
    end
    function cleanup(sim)
        kill(sim.satellite_process)
        ZMQ.close(sim.socket)
        println("Killed satellite process")
    end

    if initial_condition !== nothing
        initial_condition = copy(initial_condition)
    end

    return simulate_helper(setup, step, cleanup,
        log_init, log_step,
        terminal_condition, max_iterations,
        initial_condition, measure,
        copy(model), copy(environment),
        silent)
end

function print_iteration(i, max_iterations, state, env, sim)
    b = 1e9 * env.b
    @printf("[%d/%d]: norm(ω)=%.3f r=<%.3f %.3f %.3f> b=<%.3f %.3f %.3f>nT dt=%.3f",
        i, max_iterations, norm(state.angular_velocity), state.position[1], state.position[2], state.position[3], b[1],
        b[2], b[3], sim.dt)
end

function print_iteration(i, max_iterations)
    @printf("[%d/%d]", i, max_iterations)
end

"""
    simulate_multiple(controls::AbstractArray{Function})
    simulate_multiple(controls::AbstractArray{Function}; log_init=default_log_init, log_step=default_log_step,
        terminal_condition=default_terminate, max_iterations=1000, dt=0.1,
        initial_conditions=nothing, measures=nothing,
        models=nothing, initial_environment=default_environment)
"""
function simulate_multiple(controls; log_init=default_log_init, log_step=default_log_step,
    terminal_condition=default_terminate, max_iterations=1000, dt=0.1,
    initial_conditions=nothing, measures=nothing,
    models=nothing, initial_environment=default_environment)
    N = length(controls)
    @assert (N > 0)
    @assert (isnothing(initial_conditions) || length(initial_conditions) == N)
    @assert (isnothing(measures) || length(measures) == N)

    if isnothing(initial_conditions)
        states = [initialize_orbit() for _ = 1:N]
    else
        states = initial_conditions
    end
    if isnothing(models)
        models = [copy(default_model) for _ = 1:N]
    else
        models = models
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

        for j = 1:N
            local_env = state_view_environment(states[j], env)
            u = controls[j](measure_funcs[j](states[j], local_env))
            states[j] = integrate_state(states[j], models[j], local_env, u, dt)
        end
        update_environment(env, dt)

        log_step(hist, states)
        append!(time_hist, env.time - start_time)
        print("\r\033[K")
        print_iteration(i, max_iterations)
        if terminal_condition(states, env, i)
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
    model::Model, environment::Environment, silent::Bool)
    if isnothing(initial_condition)
        state = initialize_orbit()
    else
        state = initial_condition
    end
    if !silent
        println("intialized orbit!")
    end

    model = copy(model)

    sim = setup()

    env = copy(environment)
    start_time = env.time

    hist = log_init(state)
    time_hist = []

    try
        for i = 1:max_iterations
            local_env = state_view_environment(state, env)

            step(sim, measure(state, local_env), i)

            state = integrate_state(state, model, local_env, sim.control, sim.dt)

            update_environment(env, sim.dt)

            log_step(hist, state)
            append!(time_hist, env.time - start_time)
            if !silent
                print("\r\033[K")
                print_iteration(i, max_iterations, state, local_env, sim)
            end
            if terminal_condition(state, env, i)
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
    default_terminate(state::RBState, env::Environment, i::Int)

default termination condition for `simulate`, always returns `false`.
"""
function default_terminate(state::RBState, env::Environment, i::Int)
    return false
end

"""
    default_measure(state::RBState, env::Environment)

Default measurement function used by `simulate`, returns the state and environment.
"""
@inline function default_measure(state::RBState, env::Environment)
    return (state, env)
end

"""
    sil_measure(state::RBState, env::Environment)

Default measurement function used by `simulate` when launching satellite with a command, returns the state and environment.
"""
@inline function sil_measure(state::RBState, env::Environment)
    return Dict(
        :ω => state.angular_velocity,
        :b => env.b,
        :t => 1694097430,
    )
end

"""
    world_to_body(::RBState, ::Vector)

Returns the vector in the body frame
"""
@inline function world_to_body(state, vec)
    return quaternionToMatrix(state.attitude)' * vec
end

end