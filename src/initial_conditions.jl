using SatelliteDynamics

function state_from_osc(xosc::Vector{<:Real}, q::Vector{<:Real}, ω::Vector{<:Real})
    rv = SatelliteDynamics.sOSCtoCART(xosc, use_degrees=false)
    return RBState(rv[1:3], rv[4:6], q, ω)
end