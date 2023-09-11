try:
    from ulab.numpy import eye as identity, array, linalg, cross, dot as matmul, isfinite, all
except Exception:
    from numpy import identity, array, linalg, cross, matmul, isfinite, all

from GNCTestClient import GNCTestClient

PORT = 5555

def bcross(b, ω, k=7e-4):
    b = array(b)
    ω = array(ω)
    b_hat = b / linalg.norm(b)
    bbt = matmul(array([b_hat]).transpose(), array([b_hat]))
    M = - k * matmul(identity(3) - bbt, ω)
    # control
    m = (1 / linalg.norm(b)) * cross(b_hat, M)
    if all(isfinite(m)):
        return m.tolist()
    return [0, 0, 0]

def log(text):
    with open('log.txt', 'a') as f:
        f.write(text + '\n')

try:
    client = GNCTestClient(PORT)

    while True:
        sensors = client.uplink()
        angular_velocity = sensors['ω']
        magnetic_field = sensors['b']
        t = sensors['t']

        control = bcross(magnetic_field, angular_velocity)
        client.downlink(control)
except Exception as e:
    log(str(e))
