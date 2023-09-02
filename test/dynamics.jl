begin
    using LinearAlgebra
    using SatelliteDynamics
end

@testset "Airspeed from state" begin
    expected = [0.0, 8649.213359431995, 0.0]
    # 0.000017 seconds (7 allocations: 304 bytes)
    @time got = SP.airspeed_from_state(default_data.state)
    @test got ≈ expected
end

@testset "Drag acceleration torque" begin
    a_expected = [1.9645203449589757e-6, -6.508239236619761e-6, 9.950592930420614e-7]
    τ_expected= [1.3610453269277392e-9, 1.0339757656912846e-25, -2.6870772966569952e-9]
    @time a, τ = SP.drag_acceleration_torque(default_data.state, default_data.model, default_data.env)
    # 0.000035 seconds (111 allocations: 6.891 KiB)
    @test a ≈ a_expected
    @test τ ≈ τ_expected
end

@testset "Cartesian acceleration & torque" begin
    u = [0.0, 0.0, 0.0]
    @time a, τ = SP.cartesian_acceleration_torque(default_data.state, u, default_data.model, default_data.env)
    # 0.000059 seconds (383 allocations: 25.406 KiB)
    a_expected = [-8.696182355824059, -4.582971748165271e-5, -0.00010804973126869058]
    τ_expected = [-2.0364643283794252e-9, 2.6084352385283605e-9, -6.775772747899396e-9]
    @test a ≈ a_expected
    @test τ ≈ τ_expected

    u = [4.0, 3.0, 2.0]
    @time a, τ = SP.cartesian_acceleration_torque(default_data.state, u, default_data.model, default_data.env)
    # 0.000127 seconds (383 allocations: 25.406 KiB)
    a_expected = [-8.696182355824059, -4.582971748165271e-5, -0.00010804973126869058]
    τ_expected= [1.8457691095622258e-7, -1.3719520438198725e-7, -5.3585508411985804e-8]
    @test a ≈ a_expected
    @test τ ≈ τ_expected
end

@testset "spherical gravity" begin
    u = [0.0,0.0,0.0]
    r = 6.7751363e6
    x = SP.RBState([r
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
    env = SP.spherical_only_environment
    env = SP.state_view_environment(x, env)
    @time a, τ = SP.cartesian_acceleration_torque(x, u, model, env)
    a_expected = [-SatelliteDynamics.GM_EARTH / r^2, 0.0, 0.0]
    τ_expected = [0.0, 0.0, 0.0]
    @test τ ≈ τ_expected
    @test a ≈ a_expected
end
