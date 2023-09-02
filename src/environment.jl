struct EnvironmentConfig
  n_gravity::Int
  m_gravity::Int
  include_drag::Bool
  include_solar_radiation_pressure::Bool
  include_sun_gravity::Bool
  include_moon_gravity::Bool
  include_gravity_gradient_torque::Bool
  EnvironmentConfig(; n_gravity=10, m_gravity=10, include_drag=true,
    include_solar_radiation_pressure=true, include_sun_gravity=true,
    include_moon_gravity=true, include_gravity_gradient_torque=true) =
    new(n_gravity, m_gravity, include_drag, include_solar_radiation_pressure,
      include_sun_gravity, include_moon_gravity, include_gravity_gradient_torque)
  EnvironmentConfig(n_gravity, m_gravity, include_drag, include_solar_radiation_pressure,
    include_sun_gravity, include_moon_gravity, include_gravity_gradient_torque) =
    new(n_gravity, m_gravity, include_drag, include_solar_radiation_pressure,
      include_sun_gravity, include_moon_gravity, include_gravity_gradient_torque)
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
  e.include_drag, e.include_solar_radiation_pressure, e.include_sun_gravity,
  e.include_moon_gravity, e.include_gravity_gradient_torque)

default_environment_config = EnvironmentConfig()

default_environment = Environment(
  Epoch(2020, 11, 30),
  [0.0, 0.0, 0.0],
  default_environment_config
)

spherical_only_environment_config = EnvironmentConfig(
  n_gravity=1,
  m_gravity=1,
  include_drag=false,
  include_solar_radiation_pressure=false,
  include_sun_gravity=false,
  include_moon_gravity=false,
  include_gravity_gradient_torque=false
)

spherical_only_environment = Environment(
  Epoch(2020, 11, 30),
  [0.0, 0.0, 0.0],
  spherical_only_environment_config
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
