import time
from multiprocessing import shared_memory
import random
import msgpack
import threading
try:
    from ulab.numpy import eye as identity, array, linalg, cross, dot as matmul, isfinite, all
except Exception:
    from numpy import identity, array, linalg, cross, matmul, isfinite, all

UPLINK_FILE = "/tmp/satuplink"
DOWNLINK_FILE = "/tmp/satdownlink"

print('Starting fake cubesat')
uplink = open(UPLINK_FILE, "wb")
print('Uplink established')

MAGIC_PACKET_SIZE = 43


class cubesat:

    def __init__(self):
        self.gyro = [0.0, 0.0, 0.0]
        self.magnetic = [0.0, 0.0, 0.0]


Satellite = cubesat()


def read_data():
    try:
        read_file = open(DOWNLINK_FILE, "rb")
        data = read_file.read()
        if len(data) == 0:
            raise Exception('No downlinked data')
        data = msgpack.unpackb(data)

        Satellite.gyro = data['ω']
        Satellite.magnetic = data['b']

    except Exception as e:
        print(f'Error reading downlinked data:\n {e}')


# def bcross(b, ω, k=7e-4):
def bcross(b, ω, k=7e-4):
    b = (array(b))
    ω = (array(ω))
    b_hat = b / linalg.norm(b)
    bbt = matmul(array([b_hat]).transpose(), array([b_hat]))
    M = - k * matmul(identity(3) - bbt, ω)
    # control
    m = (1 / linalg.norm(b)) * cross(b_hat, M)
    if all(isfinite(m)):
        return m.tolist()
    return [0, 0, 0]


prev_time = time.time()

while True:
    # send control data to Julia sim
    try:
        control_data = {
            "m": bcross(Satellite.magnetic, Satellite.gyro),
            "dt": time.time() - prev_time
        }
        prev_time = time.time()
        payload = msgpack.packb(control_data, use_bin_type=True)
        # print(payload)
        if len(payload) != MAGIC_PACKET_SIZE:
            raise Exception('Invalid payload size')
            exit()
        uplink.write(payload)
        uplink.flush()
    except BrokenPipeError:
        print('Uplink broken')
        exit(BrokenPipeError)
    except Exception as e:
        print(f'Error sending uplinked data {e}')

    # read from the Julia data to update data
    t = threading.Thread(target=read_data(), name="read_thread")
    t.start()

    time.sleep(10)
