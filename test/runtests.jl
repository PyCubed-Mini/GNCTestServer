begin
    # Generic test dependencies
    using Test
    include("../src/SatellitePlayground.jl")
    SP = SatellitePlayground
    include("./utils/data.jl")
end

include("./dynamics.jl")
include("./fast_auto_test.jl")
include("./type_stability.jl")
include("./test_rb_state.jl")
# include("./test_detumble.jl")
