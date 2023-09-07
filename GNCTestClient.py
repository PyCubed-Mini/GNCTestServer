import msgpack
import zmq

class GNCTestClient:

    def __init__(self, port):
        context = zmq.Context()
        self.socket = context.socket(zmq.PAIR)
        self.socket.connect(f"tcp://localhost:{port}")

    def uplink(self):
        payload = self.socket.recv()
        return msgpack.unpackb(payload)

    def downlink(self, control):
        payload = {
            "m": control,
        }
        payload = msgpack.packb(payload, use_bin_type=True)
        self.socket.send(payload)
