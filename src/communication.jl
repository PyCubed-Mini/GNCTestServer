using ZMQ
using MsgPack

function init_zmq_socket(location)
    socket = ZMQ.Socket(ZMQ.PAIR)
    socket.rcvtimeo = 1000
    ZMQ.bind(socket, location)
    return socket
end

"""
    Send message from server to satellite
"""
function uplink(measurement, socket)
    payload = MsgPack.pack(measurement)
    ZMQ.send(socket, payload)
end

"""
    Recieve message from satellite
"""
function downlink(socket)
    raw_msg = ZMQ.recv(socket)
    return MsgPack.unpack(raw_msg)
end