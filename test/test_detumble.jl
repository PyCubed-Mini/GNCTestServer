using LinearAlgebra

@testset "detumbling" begin
    function control_law(measurement)
        (ω, b) = measurement

        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂ * b̂') * ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return SP.Control(
            m
        )
    end

    @inline function measure(state, env)
        return (state.angular_velocity, env.b)
    end
    
    model = copy(SP.pqmini_model)
    model.control_limit = [Inf, Inf, Inf]

    @time (data, time) = SP.simulate(control_law, max_iterations=10_000, 
        measure=measure, log_step=SP.angular_log_step, model=model, initial_condition=default_data.state)
    # 4.344277 seconds (15.07 M allocations: 795.512 MiB, 10.87% gc time, 2.39% compilation time)
    # 4.517416 seconds (14.91 M allocations: 780.852 MiB, 7.04% gc time, 2.35% compilation time)
    # 3.893083 seconds (15.89 M allocations: 776.579 MiB, 8.71% gc time, 1.71% compilation time)
    # 3.608502 seconds (13.51 M allocations: 711.263 MiB, 8.77% gc time, 1.50% compilation time)
    # 2.770832 seconds (13.51 M allocations: 711.281 MiB, 7.51% gc time, 1.98% compilation time)
    # 3.839175 seconds (16.68 M allocations: 1.065 GiB, 13.51% gc time, 2.75% compilation time)
    # 3.799204 seconds (16.92 M allocations: 1.079 GiB, 12.32% gc time, 3.54% compilation time)
    # 3.870619 seconds (16.91 M allocations: 1.078 GiB, 14.30% gc time, 1.44% compilation time)
    data = SP.vec_to_mat(data)
    time /= 60
    ω_final = data[end, 1:3]
    @test ω_final ≈ [-0.11135501930818302, 0.1342855195461621, -0.3456467150817579]
    # display(plot(time, data, title="DeTumbling", xlabel="Time (minutes)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))
end