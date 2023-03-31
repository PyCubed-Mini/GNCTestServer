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

test_igrf(1)
@profview test_igrf(1000)