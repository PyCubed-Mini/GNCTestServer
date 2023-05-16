struct SatelliteFace
    area::Real # (m^2)
    normal::AbstractArray{<:Real}
    position::AbstractArray{<:Real} # (m) relative to geometric center
    reflectivity_coefficient_ϵ::Real # between 0 and 1, 0 = totally absorbing, 1 = totally reflecting
end

function SatelliteFace(
    area,
    normal,
    position;
    reflectivity_coefficient_ϵ=0.1
)
    return SatelliteFace(
        area,
        normal,
        position,
        reflectivity_coefficient_ϵ
    )
end

"""
    Model

Defines how to model the satellite:
    - mass
    - inertia matrix
    - actuation type (magnetic torque coils or reaction wheels)
    - faces (for drag)
"""
struct Model
    mass::Real # kg
    inertia::AbstractArray{Float64,2} # Inertia matrix,
    control_type::Val # Type of control
    center_of_mass::AbstractArray{Float64,1} # Center of mass
    faces::AbstractArray{SatelliteFace,1} # Faces for drag
    drag_coefficient::Real # C_D, dimensionless
    max_dipoles::AbstractArray{Float64,1} # Max dipole moment for each axis
    Model(mass, inertia, control_type, center_of_mass, faces, drag_coefficient, max_dipoles) = new(mass, inertia, control_type, center_of_mass, faces, drag_coefficient, max_dipoles)
    Model(;mass, inertia, control_type, center_of_mass, faces, drag_coefficient, max_dipoles) = new(mass, inertia, control_type, center_of_mass, faces, drag_coefficient, max_dipoles) 
end


Base.copy(m::Model) = Model(
    m.mass,
    copy(m.inertia),
    m.control_type,
    copy(m.center_of_mass),
    copy(m.faces),
    copy(m.drag_coefficient),
    copy(m.max_dipoles)
)

mutable struct Control
    m::Array{Float64,1} # Control input
end

function coil_resistance(N, length, width, trace_area)
    copper_resistivity = 1.77e-8 # Ohm-meters
    perimeter = 2 * (length + width)
    return copper_resistivity * N * perimeter / trace_area
end

sat_dipole_magnitude(N, I, A) = N * I * A # N turns, I amps, A area

pqmini_faces = [
    SatelliteFace(0.06 * 0.066, [0.0, 1.0, 0.0], [0.0, 0.05 / 2, 0.0]) # Top +Y - projected area from -Y board
    SatelliteFace(0.06 * 0.066, [0.0, -1.0, 0.0], [0.0, -0.05 / 2, 0.0]) # Bottom -Y
    SatelliteFace(0.05^2, [0.0, 0.0, 1.0], [0.0, 0.0, 0.05 / 2]) # Side +Z
    SatelliteFace(0.05^2, [0.0, 0.0, -1.0], [0.0, 0.0, -0.05 / 2]) # Side -Z
    SatelliteFace(0.05^2, [1.0, 0.0, 0.0], [0.05 / 2, 0.0, 0.0]) # Side +X
    SatelliteFace(0.05^2, [-1.0, 0.0, 0.0], [-0.05 / 2, 0.0, 0.0]) # Side -X
]

pqmini_dipole_magnitude = let
    pqmini_voltage = 3.7
    pqmini_trace_width = 0.19e-3 # 7.5mil
    pqmini_trace_thickness = 0.036e-3 # 1oz copper
    pqmini_trace_area = pqmini_trace_thickness * pqmini_trace_width

    pqmini_N_turns = 154
    pqmini_coil_R = coil_resistance(pqmini_N_turns, 0.1, 0.1, pqmini_trace_area)
    pqmini_coil_I = pqmini_voltage / pqmini_coil_R

    sat_dipole_magnitude(pqmini_N_turns, pqmini_coil_I, 0.05 * 0.05)
end

pqmini_inertia_matrix = let
    inertia_CAD = [
        3.28E+05 -3853.031 1060.576
        -3853.031 3.18E+05 22.204
        1060.576 22.204 3.33E+05] # g*mm^2, determined from Fusion360 CAD model

    mass_CAD = 642.629 # g
    mass_measured = 191 # g
    g2kg = 1e3
    mm2m = 1e3
    inertia_matrix_kgm2 = (mass_measured .* inertia_CAD ./ mass_CAD) / g2kg / mm2m
end

pqmini_model = Model(
    mass=0.191,
    inertia=pqmini_inertia_matrix, # kg*m^2 determined from CAD
    control_type=Val(:dipole),
    center_of_mass=[0.0, 0.0, 0.0],
    faces=pqmini_faces,
    drag_coefficient=2.2, # CD
    max_dipoles=pqmini_dipole_magnitude * ones(3),
)
default_model = pqmini_model