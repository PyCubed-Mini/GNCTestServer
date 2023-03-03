import time
try:
    from ulab.numpy import eye as identity, array, linalg, cross, dot as matmul, isfinite, all
except Exception:
    from numpy import identity, array, linalg, cross, matmul, isfinite, all

import GNCTestClient

client = GNCTestClient.GNCTestClient()
client.register_state("control", [0.0000001, 0.000000002, 0.00000003])
client.register_state("ω", [0.1, 0.2, 0.3])
client.register_state("b", [0.1, 1.1, -0.2])


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

client.launch()

while True:
    with open('./output.txt', 'w') as f:
        f.write(f'hi+{time.time()}\n')
        f.flush()
    client["control"] = bcross(client["b"], client["ω"])
    time.sleep(2)
