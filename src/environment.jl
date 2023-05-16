struct EnvironmentConfig
  n_gravity::Int
  m_gravity::Int
  include_solar_radiation_pressure::Bool
  include_sun_gravity::Bool
  include_moon_gravity::Bool
  include_gravity_gradient_torque::Bool
end

"""
  EnvironmentState

Information about the environment state in the ECI frame.
"""
mutable struct Environment
  time::Epoch
  b::Array{Float64,1} # Magnetic field
  config::EnvironmentConfig
end


Base.copy(e::Environment) = Environment(e.time, copy(e.b), copy(e.config))
Base.copy(e::EnvironmentConfig) = EnvironmentConfig(e.n_gravity, e.m_gravity,
  e.include_solar_radiation_pressure, e.include_sun_gravity,
   e.include_moon_gravity, e.include_gravity_gradient_torque)

default_environment_config = EnvironmentConfig(10, 10, true, true, true, true)

default_environment = Environment(
  Epoch(2020, 11, 30),
  [0.0, 0.0, 0.0],
  default_environment_config
)

function update_environment(env::Environment, dt)
  env.time += dt
end

"""
  state_view_environment(state::RBState, env::Environment)

Returns the environment from the pov of the satellite.
"""
function state_view_environment(state::RBState, env::Environment)
  return Environment(
    env.time,
    world_to_body(state, IGRF13(state.position, env.time)),
    env.config
  )
end
