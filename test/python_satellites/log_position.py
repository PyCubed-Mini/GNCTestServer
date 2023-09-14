def log(text):
    with open('log.txt', 'a') as f:
        f.write(f'{text}' + '\n')


log('hello world')

try:
    import sys
    sys.path.append('..')
    from GNCTestClient import GNCTestClient  # noqa: needs to see root of project directory for client

    PORT = 5555
    f = open('/tmp/position_logger.txt', 'w')

    f.write('hello world\n')

    client = GNCTestClient(PORT)

    while True:
        sensors = client.uplink()
        position = sensors['r']

        f.write(str(position) + '\n')
        log(str(position))

        client.downlink([0, 0, 0])

except Exception as e:
    log(str(e))
