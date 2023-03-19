from threading import Thread, Lock
from multiprocessing import shared_memory
import copy
import msgpack
import time
import os
import signal
import sysv_ipc

UPLINK_SHARED_MEM = "gnc_uplink"
UPLINK_SEMAPHORE = 67
DOWNLINK_SHARED_MEM = "gnc_downlink"
DOWNLINK_SEMAPHORE = 68

TIME_INTERVAL = 10


class ThreadSafeState:

    def __init__(self, state):
        self.state = state
        self.lock = Lock()

    def get(self):
        with self.lock:
            return copy.copy(self.state)

    def set(self, state):
        with self.lock:
            self.state = state


class GNCTestClient:

    def __init__(self, uplink_interval=2.0, downlink_interval=1.0):
        try:
            self._state = {}
            self._prev_time = time.time()
            self._log_fd = open("/tmp/satlog.txt", "w")
            self._uplink_interval = uplink_interval
            self._downlink_interval = downlink_interval

            self._uplink_shm = shared_memory.SharedMemory(name=UPLINK_SHARED_MEM)
            self._uplink_buf = self._uplink_shm.buf
            self._uplink_sem = sysv_ipc.Semaphore(UPLINK_SEMAPHORE)
            self._uplink_id = 1

            self._downlink_shm = shared_memory.SharedMemory(
            name=DOWNLINK_SHARED_MEM)
            self._downlink_buf = self._downlink_shm.buf
            self._downlink_sem = sysv_ipc.Semaphore(DOWNLINK_SEMAPHORE)
        except Exception as e:
            self.log(f'Error initializing client:\n {e}')
            exit(1)

    def register_state(self, key, value):
        self._state[key] = ThreadSafeState(value)

    def __getitem__(self, key):
        return self._state[key].get()

    def __setitem__(self, key, value):
        self._state[key].set(value)

    def downlink(self):
        with self._downlink_sem:
            size = self._downlink_buf[0]
            if size != 0:
                data = self._downlink_buf[1:size+1]
                data = msgpack.unpackb(data)
                for key, value in data.items():
                    self[key] = value

    def uplink(self):
        control = {
            "m": self["control"],
            "dt": time.time() - self._prev_time
        }
        payload = msgpack.packb(control, use_bin_type=True)

        self._prev_time = time.time()

        with self._uplink_sem:
            self._uplink_buf[0:8] = self._uplink_id.to_bytes(8, 'little')
            self._uplink_buf[8:len(payload)+8] = payload
            self._uplink_id += 1

    def communication_thread(self):
        downlink_time = time.time()
        uplink_time = time.time()
        while True:
            if time.time() > downlink_time:
                self.downlink()
                downlink_time += self._downlink_interval
            if time.time() > uplink_time:
                self.uplink()
                uplink_time += self._uplink_interval

            next_uplink = uplink_time - time.time()
            next_downlink = downlink_time - time.time()

            sleep_time = min(next_uplink, next_downlink)

            if sleep_time < 0:
                self.log(f'Warning: sleep time is negative: {sleep_time}')
                os.kill(os.getpid(), signal.SIGINT)

            time.sleep(sleep_time)
    
    def com_thread_wrapper(self):
        try:
            self.communication_thread()
        except Exception as e:
            self.log(e)

    def launch(self):
        time.sleep(10.0)
        self.log('Launching communication thread...')
        t = Thread(target=self.com_thread_wrapper,
                   name="communication_thread")
        t.start()

    def log(self, text):
        self._log_fd.write(text)
        self._log_fd.write("\n")
        self._log_fd.flush()
