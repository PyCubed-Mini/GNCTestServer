using StaticArrays
using SatelliteDynamics

""" 
    dynamics(state::RBState, parameters::Parameters, environment::Environment, control::Control, t::Epoch)::RBState

You can configure which types of effects are taken into account (J2, SRP, etc...) by tweaking the environment.
You can configure the satellite's parameters such as intertia by tweaking the parameters.
"""
function dynamics(state::RBState, model::Model, environment::Environment, control::Control)::RBState
    a, τ = cartesian_acceleration_torque(state, control.m, model, environment)
    ω = state.angular_velocity
    J = model.inertia
    return RBState{eltype(state)}(
        state.velocity,
        a,
        qdot(state.attitude, state.angular_velocity),
        J \ (τ - cross(ω, J * ω))
    )
end

function control(x::RBState, ::Val{:torque}, u::AbstractArray, model::Model, env::Environment)
    return zeros(3), u
end

function control(x::RBState, ::Val{:thrust}, u::AbstractArray, model::Model, env::Environment)
    return u, zeros(3)
end

function control(x::RBState, ::Val{:dipole}, u::AbstractArray, model::Model, env::Environment)
    τ_mag = magnetic_torque(x, u, model, env)
    return zeros(3), τ_mag
end

"""
    cartesian_acceleration_torque(x::RBSTATE, params::OrbitSimulationParameters, epc::Epoch)

Compute the acceleration and torque experienced by the satellite
"""
function cartesian_acceleration_torque(x::RBState, u::AbstractArray{<:Real}, model::Model, env::Environment)
    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = SatelliteDynamics.bias_precession_nutation(env.time)
    E = SatelliteDynamics.earth_rotation(env.time)
    W = SatelliteDynamics.polar_motion(env.time)
    R = W * E * PN

    body_to_world = quaternionToMatrix(x.attitude)

    # eltype(x) makes this forward diff friendly
    a = @MVector zeros(eltype(x), 3) # inertial frame acceleration
    τ = @MVector zeros(eltype(x), 3) # body frame torque

    # spherical harmonic gravity
    a += SatelliteDynamics.accel_gravity(x.position, R, env.config.n_gravity, env.config.m_gravity)

    # add control torque and acceleration
    ca, cτ = control(x, model.control_type, u, model, env)
    a += ca
    τ += cτ

    # drag
    if env.config.include_drag
        a_drag_body, τ_drag = drag_acceleration_torque(x, model, env) # body frame
        a_drag_inertial = body_to_world * a_drag_body
        τ += τ_drag
        a += a_drag_inertial
    end

    # Solar radiation pressure
    if env.config.include_solar_radiation_pressure
        a_srp_body, τ_srp = solar_acceleration_torque(x, model, env)
        a_srp_inertial = body_to_world * a_srp_body
        a += a_srp_inertial
        τ += τ_srp
    end

    # gravity gradient torque
    if env.config.include_gravity_gradient_torque
        τ_gg = gravity_gradient_torque(x, model)
        τ += τ_gg
    end

    # third body sun
    if env.config.include_sun_gravity
        r_sun = SatelliteDynamics.sun_position(env.time)
        a += SatelliteDynamics.accel_thirdbody_sun(x, r_sun)
    end

    # third body moon
    if env.config.include_moon_gravity
        r_moon = SatelliteDynamics.moon_position(env.time)
        a += SatelliteDynamics.accel_thirdbody_moon(x, r_moon)
    end

    return a, τ
end

"""
    solar_acceleration_torque(x::Vector{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)

Compute the acceleration and torque due to solar radiation pressure on each face.
Assumes the SRP center of pressure is the geometric center of each face.
Does not account for self-shielding of faces, assumes constant distance between earth and sun (P_sun is constant).
See Montenbruk and Gill eq 3.73 and Markley and Crassidis eq 3.167
"""
function solar_acceleration_torque(x::RBState, model::Model, env::Environment)

    r_sun = SatelliteDynamics.sun_position(env.time)
    r_sun_body = quaternionToMatrix(x.attitude)'r_sun
    s = r_sun_body / norm(r_sun_body)

    a = zeros(eltype(x), 3)
    τ = zeros(eltype(x), 3)

    nu = SatelliteDynamics.eclipse_conical([x.position; x.velocity], r_sun)
    P_sun = SatelliteDynamics.P_SUN

    for f in model.faces
        cθ = s'f.normal
        ϵ = f.reflectivity_coefficient_ϵ
        if cθ > 0.0
            force_srp = -nu * P_sun * f.area * cθ * ((1 - ϵ) * s + 2 * ϵ * cθ * f.normal)
            τ += cross(f.position + model.center_of_mass, force_srp)
            a += force_srp / model.mass
        end
    end

    return a, τ
end

"""
    drag_acceleration_torque(airspeed, epc, params)

Compute the resultant forces and torques on a spacecraft due to drag.
See Markley and Crassidis pg 108

All units are SI kg-m-s, return values are in the body-frame
"""
function drag_acceleration_torque(x::RBState, model::Model, env::Environment)
    q = x.attitude
    r = x.position

    airspeed_inertial = airspeed_from_state(x)
    airspeed = quaternionToMatrix(q)' * airspeed_inertial

    force = zeros(eltype(airspeed), 3)
    torque = zeros(eltype(airspeed), 3)

    ρ = SatelliteDynamics.density_harris_priester(env.time, r)
    CD = model.drag_coefficient
    CoM = model.center_of_mass

    for f in model.faces
        if airspeed'f.normal > 0.0
            f_force = -0.5 * ρ * CD * f.area * (airspeed'f.normal) * airspeed
            force += f_force
            torque += cross(f.position + CoM, f_force)
        end
    end
    accel = force / model.mass
    return accel, torque
end

"""
    airspeed_from_state(x)

Compute the inertial-frame airspeed wrt to the atmosphere given the satellite state `x`
"""
function airspeed_from_state(x::RBState)
    ω_atmosphere = [0, 0, SatelliteDynamics.OMEGA_EARTH]
    airspeed_inertial = x.velocity + cross(ω_atmosphere, x.position)

    return airspeed_inertial
end

"""
    magnetic_torque()

Returns the magnetic torque on the satellite due to the magnetic field of the earth.
"""
function magnetic_torque(x::RBState, u::AbstractArray{<:Real}, model::Model, env::Environment)
    τ_mag = cross(u, env.b)
    return τ_mag
end

"""
    gravity_gradient_torque

    Computes the gravity gradient torque on the satellite.
    See "Fundementals of Spacecraft Attitude Determination and Control" by Markley and Crassidis 3.3.6.1
"""
function gravity_gradient_torque(x::AbstractArray{<:Real}, model::Model)
    μ = SatelliteDynamics.GM_EARTH

    r_inertial = x.position 
    q_b2i = x.attitude

    r_body = quaternionToMatrix(q_b2i)'r_inertial

    r_mag = norm(r_body)
    n = -r_body / r_mag # body frame nadir pointing vector

    J = model.inertia

    τ_gg = (3 * μ / r_mag^3) * hat(n) * J * n

    return τ_gg
end