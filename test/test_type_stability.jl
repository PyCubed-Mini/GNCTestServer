begin
  include("deps.jl")
end
begin
  x = SP.RBState([6.7751363e6
    0.0
    0.0
    0.0
    8155.162619089651
    0.0
    -0.9442165245970416
    -0.11229349782424022
    0.28649111413323164
    0.117337830843173
    0.013209720675656584
    0.36976717864576536
    -0.548490925366443]
  )
  model = SP.default_model
  env = SP.default_environment
  env = SP.state_view_environment(x, env)
end

@testset "qdot" begin
  state = SP.initialize_orbit()
  Test.@inferred SP.qdot(state.attitude, state.angular_velocity)
end

@testset "dynamics" begin
  state = SP.initialize_orbit()
  u = SP.Control([0.0, 0.0, 0.0])
  Test.@inferred SP.dynamics(state, model, env, u)
end

@testset "rk4" begin
  state = SP.initialize_orbit()
  u = SP.Control([0.0, 0.0, 0.0])
  dt = 0.1

  Test.@inferred SP.rk4(state, dt, model, env, u)
end

@testset "Integrate state" begin
  state = SP.initialize_orbit()
  u = SP.Control([0.0, 0.0, 0.0])
  dt = 0.1
  Test.@inferred SP.integrate_state(state, model, env, u, dt)
end

@testset "Airspeed from state" begin
  Test.@inferred SP.airspeed_from_state(default_data.state)
end

# @testset "drag_acceleration_torque" begin
#   Test.@inferred SP.drag_acceleration_torque(default_data.state, default_data.model, default_data.env)
# end

@testset "density_harris_priester" begin
  Test.@inferred SP.SatelliteDynamics.density_harris_priester(default_data.env.time, default_data.state.position)
end
