import time
import threading
from random import randint
from multiprocessing import shared_memory

control_data = 0
UP_FIFO = "/tmp/satuplink"
DOWN_FIFO = "/tmp/satdownlink"

# Add threading to the random reads


def read_file():
    # read from the Julia data to update data
    # print("Starting data read...")
    try:
        read_file = open(DOWN_FIFO, "r")
        # print("file opened")
        print(read_file.read())

    except Exception as e:
        print(f"Failed to read data: {e}")


while True:
    # send control data to Julia sim
    # print("opening uplink...")
    try:
        write_file = open(UP_FIFO, "w")
        # print("FIFO opened")
        write_file.write(time.time(), control_data)
        # print("written to file")
    except Exception as e:
        print(f"Failed to Open satuplink: {e}")

    if randint > 0.5:
        # read from the Julia data to update data
        # print("Starting data read...")
        try:
            read_file = open(DOWN_FIFO, "r")
            # print("file opened")
            print(read_file.read())

        except Exception as e:
            print(f"Failed to read data: {e}")

    time.sleep(0.1)
