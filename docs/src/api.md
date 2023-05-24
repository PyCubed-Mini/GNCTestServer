```@meta
CurrentModule = SatellitePlayground
```

# SatellitePlayground.jl API

The API all revolves around the `simulate` function which has two modes of running:
- Software in the loop: A separate program pretends to be a satellite, and the simulator and program communicate sensor/control data between them.
    This is useful for testing your true satellite software to make sure the guidance, navigation, and controll code still works despite all the other things it is doing.
- Control function: A function is passed a measurement (by default the state, parameters, and time in body frame) and returns the magnetorquer control input.
    This is useful for directly testing attitude control schemes.
```@docs
simulate
```

## Other Simulation Modes
For simulating multiple satellites at once use:
```@docs
simulate_multiple
```

## Dynamics
There are two ways to change the dynamics of the simulation:
- Change the invidual satellite model (mass, inertia, control type, etc.)
- Change the environment model (e.g. earth gravity model, solar radiation pressure model, etc.)

## State
```@docs
RBState
```



## Defaults
```@docs
default_measure
default_terminate
```

## Useful Helper Functions
```@docs
initialize_orbit
world_to_body
vec_to_mat
```


## 