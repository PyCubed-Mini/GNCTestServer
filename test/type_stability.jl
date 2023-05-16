using Test
include("../src/SatellitePlayground.jl")
SP = SatellitePlayground

@testset "qdot" begin
  state = SP.initialize_orbit()
  Test.@inferred SP.qdot(state.attitude, state.angular_velocity)
end


@testset "cartesian_acceleration_torque" begin
  state = SP.initialize_orbit()
  env = copy(SP.default_environment)
  model = copy(SP.default_model)
  Test.@inferred SP.cartesian_acceleration_torque(state, [0.0, 0.0, 0.0], model, env)
end


@testset "dynamics" begin
  state = SP.initialize_orbit()
  params = copy(SP.default_parameters)
  env = copy(SP.default_environment)
  u = SP.Control([0.0, 0.0, 0.0])
  Test.@inferred SP.dynamics(state, params, env, u)
end

@testset "rk4" begin
  state = SP.initialize_orbit()
  params = copy(SP.default_parameters)
  env = copy(SP.default_environment)
  u = SP.Control([0.0, 0.0, 0.0])
  dt = 0.1

  Test.@inferred SP.rk4(state, env.time, dt, params, env, u)
end

@testset "Integrate state" begin
  state = SP.initialize_orbit()
  params = copy(SP.default_parameters)
  env = copy(SP.default_environment)
  u = SP.Control([0.0, 0.0, 0.0])
  dt = 0.1
  Test.@inferred SP.integrate_state(state, params, env, u, dt)

end
