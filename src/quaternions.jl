# [Simulator/quaternions.jl]

""" quaternions.jl 

    Contains common functions required when working with unit quaternions for attitude 
"""

using LinearAlgebra     # For I, norm 

""" hat(v)

    Converts a 3-element vector into a cross-product matrix.

    Arguments:
     - v:  3-element vector                   |  [3,] 

    Returns:
     - M:  A [3 × 3] skew-symmetric matrix    |  [3, 3]
"""
function hat(v::AbstractArray{<:Real,1})
    z = zero(eltype(v))
    M = [z -v[3] v[2]
        v[3] z -v[1]
        -v[2] v[1] z]

    return M
end

""" L(q) 

      Converts a scalar-first unit quaternion into the left-side matrix for 
    quaternion multiplication, as described in "Planning with Attiude" (Jackson)

    Arguments:
     - q:  A scalar-first unit quaternion                         |  [4,]

    Returns: 
     - M:  Left-side matrix representing the given quaternion     |  [4, 4]
"""
function L(q::AbstractArray{<:Real,1})
    qₛ, qᵥ = q[1], q[2:4]

    temp = qₛ * I(3) + hat(qᵥ)
    M = zeros(eltype(q), 4, 4)
    M[1, 1] = qₛ
    M[2:4, 1] = qᵥ
    M[1, 2:4] = -qᵥ
    M[2:4, 2:4] = temp
    # M = [qₛ -qᵥ'
    #      qᵥ temp]

    return M
end


""" R(q) 

      Converts a scalar-first unit quaternion into the right-side matrix for 
    quaternion multiplication, as described in "Planning with Attiude" (Jackson)

    Arguments:
     - q:  A scalar-first unit quaternion                         |  [4,]

    Returns: 
     - M:  Right-side matrix representing the given quaternion     |  [4, 4]
"""
function R(q::AbstractArray{<:Real,1})
    qₛ, qᵥ = q[1], q[2:4]

    M = [qₛ -qᵥ'
        qᵥ qₛ*I(3)-hat(qᵥ)]

    return M
end

""" qdot(q, ω)

    Provides the derivative of a quaternion, given the angular velocity 

    Arguments:
     - q:  Scalar-first unit quaternion                              |  [4,]
     - ω:  Angular velocity                                          |  [3,]

    Returns:
     - q̇:  Derivative of attitude, parameterized as a quaternion     |  [4,]
"""
function qdot(q::AbstractArray{<:Real,1}, ω::AbstractArray{<:Real,1})
    q̇ = 0.5 * L(q) * H * ω
    return q̇
end


const H = [zeros(1, 3); I(3)]     # Converts from a 3-element vector to a 4-element vector with 0 scalar part 
const T = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]   # Forms the conjugate of q, i.e. q† = Tq   


⊙(q₁, q₂) = L(q₁) * q₂   # Hamiltonian product 

""" quaternionToMatrix (q)
    Arguments:
     - q: Scalar-first unit quaternion                               | [4,]

    Returns:
     - Q: Rotation matrix representing the same rotation
"""
function quaternionToMatrix(q::AbstractArray{<:Real,1})
    s, v = q[1], q[2:4]
    return I(3) + 2 * hat(v) * (s * I(3) + hat(v))
end

function qErr(q₁::AbstractArray{<:Real,1}, q₂::AbstractArray{<:Real,1})
    return norm((L(q₁)'*q₂)[2:4])
end

function eulerError(e1, e2)
    return acos(dot(e1, e2) / (norm(e1) * norm(e2)))
end