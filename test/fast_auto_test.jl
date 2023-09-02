begin
    using LinearAlgebra
end

function no_control(measurement)
    return zero(SP.Control)
end

@testset "logging states" begin
    function log_step(hist, state)
        push!(hist, state)
    end
    (data, time) = SP.simulate(no_control, max_iterations=100, log_step=log_step, dt=10.0, initial_condition=default_data.state)
    @test length(data) == 100
    @test length(time) == 100
    @test data[1] isa SP.RBState

    (data, time) = SP.simulate(no_control, max_iterations=100, log_step=log_step, dt=10.0, initial_condition=default_data.state)
    data = SP.vec_to_mat(data)
    @test size(data) == (100, 13)
    @test length(time) == 100
end

@testset "default logging" begin
    (data, time) = SP.simulate(no_control, max_iterations=100, initial_condition=default_data.state)
    @test length(data) == 100
    @test length(time) == 100
end

@testset "logging vectors" begin
   function log_step(hist, state) 
    push!(hist, norm(state.angular_velocity))
   end
   (data, time) = SP.simulate(no_control, max_iterations=100, log_step=log_step, dt=10.0, initial_condition=default_data.state)
   @test length(data) == 100
   @test length(time) == 100
end