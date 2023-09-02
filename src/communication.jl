using ZMQ
using MsgPack

function init(location)
    socket = ZMQ.Socket(ZMQ.REP)
    ZMQ.bind(socket, location)
    return socket
end

"""
    Send message from server to satellite
"""
function downlink(measurement, socket)
    sensors = Dict(
        :ω => measurement[1].angular_velocity,
        :b => measurement[2].b,
    )
    payload = MsgPack.pack(sensors)
    ZMQ.send(socket, payload)
end

"""
    Recieve message from satellite
"""
function uplink(socket)
    raw_msg = ZMQ.recv(socket)
    return MsgPack.unpack(raw_msg)
end