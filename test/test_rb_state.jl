using Test
using StaticArrays
using BenchmarkTools
using LinearAlgebra

include("../src/satellite_models.jl")

r = @SVector rand(3)
v = @SVector rand(3)
q = @SVector rand(4)
q /= norm(q)
ω = @SVector rand(3)

x = RBState{Float64}(r, v, q, ω)
@test x[1:3] ≈ r ≈ x.position
@test x[4:6] ≈ v ≈ x.velocity
@test x[7:10] ≈ q ≈ x.attitude
@test x[11:13] ≈ ω ≈ x.angular_velocity
@test length(x) == 13

x32 = RBState{Float32}(r, v, q, ω)
@test x32 isa RBState{Float32}
@test eltype(x32[1:3]) == Float32
@test x32.position isa SVector{3,Float32}
@test x32.attitude isa SVector{4,Float32}

# Pass in other types of vectors
x = RBState{Float64}(Vector(r), v, q, ω)
@test x.position isa SVector{3,Float64}
x = RBState{Float64}(Vector(r), Vector(v), q, Vector(ω))
@test x.position isa SVector{3,Float64}
@test x.velocity isa SVector{3,Float64}
@test x.angular_velocity isa SVector{3,Float64}

x = RBState{Float64}(view(r, :), view(Vector(v), 1:3), q, Vector(ω))
@test x.position isa SVector{3,Float64}

@test RBState{Float32}(x) isa RBState{Float32}
@test RBState{Float64}(x32) isa RBState{Float64}
@test RBState(x) isa RBState{Float64}
@test RBState(x32) isa RBState{Float32}

# Let the constructor figure out data type
@test RBState(r, v, q, ω) isa RBState{Float64}
@test RBState(Float32.(r), Float32.(v), Float32.(q), Float32.(ω)) isa RBState{Float32}
@test RBState(Float32.(r), Float32.(v), q, Float32.(ω)) isa RBState{Float64}
@test RBState(r, Float32.(v), Float32.(q), Float32.(ω)) isa RBState{Float64}


# # Test comparison (with double-cover)
# q = rand(UnitQuaternion)
# x1 = RBState(r, v, q, ω)
# x2 = RBState(r, v, -q, ω)
# @test x1[7:11] ≈ -x2[7:11]
# @test x1 ≈ x2
# @test !(SVector(x1) ≈ SVector(x2))

# Test indexing
x = RBState(r, v, q, ω)
@test x[1:3] ≈ r
@test x[4:6] ≈ v
@test x[7:10] ≈ q
@test x[11:13] ≈ ω

@test getindex(x, 2) == r[2]
@test getindex(x, 4) == v[1]
@test getindex(x, 13) == ω[3]

# Renorm
x2 = RBState(r, v, 2 * q, ω)
@test norm(x2.attitude) ≈ 2
x = renorm(x2)
@test norm(x2.attitude) ≈ 2
@test norm(x.attitude) ≈ 1


# Addition and Subtraction
x1 = RBState(rand(13))
x2 = RBState(rand(13))
x = x1 + x2
@test x.position ≈ x1.position + x2.position
@test x.attitude ≈ x1.attitude + x2.attitude
@test x.velocity ≈ x1.velocity + x2.velocity
@test x.angular_velocity ≈ x1.angular_velocity + x2.angular_velocity
x = x1 - x2
@test x.position ≈ x1.position - x2.position
@test x.attitude ≈ x1.attitude - x2.attitude
@test x.velocity ≈ x1.velocity - x2.velocity
@test x.angular_velocity ≈ x1.angular_velocity - x2.angular_velocity