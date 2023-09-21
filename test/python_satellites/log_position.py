def log(text):
    with open('log.txt', 'a') as f:
        f.write(f'{text}' + '\n')


try:
    FILE = '/tmp/sat_log.json'

    SIMULATION_ITERATIONS = 10

    import sys
    sys.path.append('..')
    from GNCTestClient import GNCTestClient  # noqa: needs to see root of project directory for client
    import json

    PORT = 5555

    client = GNCTestClient(PORT)

    i = 0
    position_history = []

    while True:
        i += 1

        sensors = client.uplink()
        position = sensors['r']
        position_history.append(position)
        if i == SIMULATION_ITERATIONS:
            with open(FILE, "w") as f:
                f.write(json.dumps(position_history))

        client.downlink([0, 0, 0])

except Exception as e:
    log(str(e))
