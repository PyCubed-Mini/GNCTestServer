import time
import os
from random import randint
from multiprocessing import shared_memory

control_data = 0
UP_FIFO = "/tmp/satuplink"
DOWN_FIFO = "/tmp/satdownlink"

while True:
    # send control data to Julia sim
    print("opening uplink...")
    try:
        write_file = open(UP_FIFO, "w")
        print("FIFO opened")
        write_file.write(time.time(), control_data)
        print("written to file")
    except Exception as e:
        print(e)

    if randint > 0.5:
        # read from the Julia data to update data
        try:
            read_file = open(DOWN_FIFO, "r")
            print(read_file.read())

        except Exception as e:
            print(e)

    time.sleep(0.1)
