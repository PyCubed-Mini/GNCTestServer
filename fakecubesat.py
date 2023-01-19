import time
from multiprocessing import shared_memory
import random
import msgpack


control_data = {
    "m": [0.0, 0.0, 0.0],
    "time": time.time()
}
UPLINK_FILE = "/tmp/satuplink"
DOWNLINK_FILE = "/tmp/satdownlink"

print('Starting fake cubesat')
uplink = open(UPLINK_FILE, "wb")
print('Uplink established')

while True:
    # send control data to Julia sim
    try:
        payload = msgpack.packb(control_data, use_bin_type=True)
        print(payload)
        uplink.write(payload)
        uplink.flush()
    except Exception as e:
        print(f'Error reading uplinked data {e}')

    if random.uniform(0, 1) > 0.5:
        # read from the Julia data to update data
        try:
            read_file = open(DOWNLINK_FILE, "rb")
            data = msgpack.unpack(read_file)
            print(data)

        except Exception as e:
            print(f'Error reading downlinked data {e}')

    time.sleep(2.0)
