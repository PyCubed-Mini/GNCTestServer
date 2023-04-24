using StaticArrays
using LinearAlgebra

export RBState

"""
Based on: https://github.com/RoboticExplorationLab/RobotDynamics.jl/blob/master/src/rbstate.jl


    RBState{T} <: StaticVector{13,T}

Represents the state of a rigid body in 3D space, consisting of position, linear_velocity, attitude, 
    and angular velocity, respresented as a vector stacked in that order, with
    the rotation represented as the 4 elements unit quaternion.

Implements the `StaticArrays` interface so can be treated as an `SVector` with additional
    methods.

# Constructors
    RBState{T}(r, v, q, ω)
    RBState{T}(x)
    RBState(r, v, q, ω)
    RBState(x)

where `r`, `v`, and `ω` are three-dimensional vectors, `q` is either a `Rotation` or a
    four-dimenional vector representing the parameters of unit quaternion, and `x` is a
    13-dimensional vector (or tuple),
"""
struct RBState{T} <: StaticVector{13,T}
    position::SVector{3,T}
    velocity::SVector{3,T}
    attitude::SVector{4,T}
    angular_velocity::SVector{3,T}
    function RBState{T}(r, v, q, ω) where {T}
        @assert length(r) == 3
        @assert length(v) == 3
        @assert length(q) == 4
        @assert length(ω) == 3
        new{T}(r, v, q, ω)
    end
    @inline function RBState{T}(x::RBState) where {T}
        RBState{T}(x.position, x.velocity, x.attitude, x.angular_velocity)
    end
end

function RBState(r::AbstractVector, v::AbstractVector, q::AbstractVector, ω::AbstractVector)
    T = promote_type(eltype(r), eltype(v), eltype(q), eltype(ω))
    RBState{T}(r, v, q, ω)
end

function RBState(x::AbstractVector)
    @assert length(x) == 13
    T = eltype(x)
    RBState{T}(x[1:3], x[4:6], x[7:10], x[11:13])
end

@inline RBState(x::RBState) = x

@generated function RBState(::Type{R}, x::AbstractVector) where {R<:AbstractVector{3}}
    ir = SA[1, 2, 3]
    iv = SA[4, 5, 6]
    iq = SA[7, 8, 9]
    iω = SA[10, 11, 12]
    if R <: UnitQuaternion
        iq = SA[7, 8, 9, 10]
        iω = iω .+ 1
    end
    quote
        q = UnitQuaternion(R(x[$iq]))
        RBState(x[$ir], x[$iv], q, x[$iω])
    end
end


# Static Arrays interface
function (::Type{RB})(x::NTuple{13}) where {RB<:RBState}
    RB(
        SA[x[1], x[2], x[3]],
        SA[x[4], x[5], x[6]],
        x[7], x[8], x[9], x[10],
        SA[x[11], x[12], x[13]]
    )
end
Base.@propagate_inbounds function Base.getindex(x::RBState, i::Int)
    if i < 4
        x.position[i]
    elseif i < 7
        x.velocity[i-3]
    elseif i < 11
        x.attitude[i-6]
    else
        x.angular_velocity[i-10]
    end
end
Base.Tuple(x::RBState) = (
    x.position[1], x.position[2], x.position[3],
    x.velocity[1], x.velocity[2], x.velocity[3],
    x.attitude.w, x.attitude.x, x.attitude.y, x.attitude.z,
    x.angular_velocity[1], x.angular_velocity[2], x.angular_velocity[3]
)

"""
    renorm(x::RBState)

Re-normalize the unit quaternion.
"""
@inline renorm(x::RBState) = RBState(x.position, x.velocity, x.attitude / norm(x.attitude), x.angular_velocity)

function Base.isapprox(x1::RBState, x2::RBState; kwargs...)
    isapprox(x1.position, x2.position; kwargs...) &&
        isapprox(x1.velocity, x2.velocity; kwargs...) &&
        isapprox(x1.angular_velocity, x2.angular_velocity; kwargs...) &&
        isapprox(x1.attitude, x2.attitude; kwargs...)
end

"""
    +(::RBState, ::RBState)

Add two rigid body states, which adds the position, linear and angular velocities, and
    composes the orientations.
"""
function Base.:+(s1::RBState, s2::RBState)
    RBState(s1.position + s2.position, s1.velocity + s2.velocity, s1.attitude + s2.attitude, s1.angular_velocity + s2.angular_velocity)
end

"""
    -(::RBState, ::RBState)

Substract two rigid body states, which substracts the position, linear and angular velocities,
    and composes the inverse of the second orientation with the first, i.e. `inv(q2)*q1`.
"""
function Base.:-(s1::RBState, s2::RBState)
    RBState(s1.position - s2.position, s1.velocity - s2.velocity, s1.attitude - s2.attitude, s1.angular_velocity - s2.angular_velocity)
end


function Base.:*(s::RBState, a::Number)::RBState
    RBState(s.position * a, s.velocity * a, s.attitude * a, s.angular_velocity * a)
end

function Base.:*(a::Number, s::RBState)::RBState
    RBState(s.position * a, s.velocity * a, s.attitude * a, s.angular_velocity * a)
end
