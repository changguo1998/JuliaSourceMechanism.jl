# JuliaSourceMechanism

## Home

## Manual

### Data Type and Structure

1. `env`

```text
env
    *dataroot
    algorithm
    event
        depth
        latitude
        longitude
        magnitude
        origintime
    stations
        [1]
            network
            station
            component
            base_azimuth
            *base_begintime
            base_distance
            *base_record
            base_trim
            green_dt
            green_m
            green_tsource
            meta_lon
            meta_lat
            meta_el
            meta_dt
            meta_btime
            meta_file
            phases
                [1]
                    type
                    at
                    tt
                [2]
                ...
        [2]
        ...
```

## API
