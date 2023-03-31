using GNCTestServer
using ProfileView
using Cthulhu
using StaticArrays
using SatelliteDynamics

const position = SVector(-7.1371363e6, -0.0, -0.0)
const time = SatelliteDynamics.Epoch("2020-03-31")

function test_igrf(iter=1000)
    for _ in 1:iter
        GNCTestServer.IGRF13(position, time)
    end
end

@assert [-8.529074752279366, 2.6700487774761283, 18.698375602749916] ≈ GNCTestServer.IGRF13(position, time)
test_igrf(1)
@time test_igrf(1000)
#   0.023724 seconds (454.00 k allocations: 16.876 MiB)
