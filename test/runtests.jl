include("deps.jl")

# for some reason this test fails when launched from runtests.jl, and passes when run via `julia --project=. test/test_io.jl`  
# This likely is due to something with our test environment.
# include("./test_io.jl")

include("./dynamics.jl")
include("./fast_auto_test.jl")
include("./type_stability.jl")
include("./test_rb_state.jl")
# include("./test_detumble.jl")
