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

## State
```@docs
RBState
```

## Useful Helper Functions
```@docs
word_to_body
initialize_orbit
```


## 