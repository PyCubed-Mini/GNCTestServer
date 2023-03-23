using InterProcessCommunication
using PyCall
sysv_ipc = PyNULL()

function __init__()
    copy!(sysv_ipc, pyimport("sysv_ipc")) # how to import this / make it a dependency
end

MAGIC_PACKET_SIZE = 43

function mk_shared(name, size)
    try
        shm = IPC.SharedMemory(name)
        rm(shm)
    catch
        # This should happen
    end
    shm = IPC.SharedMemory(name, size, perms=0o777)
    ptr = convert(Ptr{UInt8}, pointer(shm))
    buf = unsafe_wrap(Array, ptr, sizeof(shm))
    return buf, shm
end

function mk_semaphore(key)
    try
        sem = sysv_ipc.Semaphore(key)
        sem.remove()
        println("Semaphore should not have existed")
    catch
        # This should happen
    end
    sem = sysv_ipc.Semaphore(key, flags=sysv_ipc.IPC_CREAT)
    sem.release() # intialize it as available
    return sem
end


function downlink(buf, buf_sem, state, params)
    @pywith buf_sem as _ begin
        sensors = Dict(
            :ω => state.angular_velocity,
            :b => params.b,
        )
        payload = MsgPack.pack(sensors)
        if length(payload) + 1 > length(buf)
            psize = length(payload)
            throw(error("Payload of size $psize too large for downlink buffer"))
        end
        buf[1] = length(payload)
        buf[2:1+length(payload)] = payload
    end
end

function uplink(buf, buf_sem, itteration)
    start = time()
    while true
        id = reinterpret(Int64, buf[1:8])[1]
        if id == itteration
            break
        elseif id >= itteration
            throw(error("Simulation speed to high, satellite ahead of server (sat: $id, server: $itteration)"))
        end
        sleep(0.0001)
    end
    @pywith buf_sem as _ begin
        println(" $(time() - start)s")
        payload = buf[9:end]
        res = MsgPack.unpack(payload)
        return res
    end
end