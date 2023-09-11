default_data_state = SP.RBState([6.7751363e6
    0.0
    0.0
    0.0
    8155.162619089651
    0.0
    -0.9442165245970416
    -0.11229349782424022
    0.28649111413323164
    0.117337830843173
    0.013209720675656584
    0.36976717864576536
    -0.548490925366443])

default_data = (
    state = default_data_state,
    model = SP.default_model,
    env = copy(SP.state_view_environment(default_data_state, SP.default_environment))
)
