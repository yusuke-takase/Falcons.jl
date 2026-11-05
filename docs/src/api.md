# API reference

```@docs
Falcons
```

## Scan description

```@docs
SpinningScanningStrategy
satellite_preset
satellite_presets
InstrumentInfo
bore2spin_quat
DetectorInfo
detector_quat
FreqChannelInfo
get_boresight_detector
IdealHWP
hwp_angle
Observation
split_observation
get_times
get_time_span
n_detectors
```

## Rotations

```@docs
Quaternion
quat_rotation
quat_rotation_x
quat_rotation_y
quat_rotation_z
rotate_vector
RotQuaternion
slerp
is_constant
time_range
spin_to_ecliptic
```

## Pointings and maps

```@docs
get_pointings
compute_pointing_and_orientation
orientation_angle
get_pointing_pixels
PointingPixels
get_hitmap
make_h_maps
HMapsResult
h_map
hitmap
nm_grid
save_h_maps
load_h_maps_from_file
```

## Coordinates

```@docs
CoordinateSystem
Ecliptic
Galactic
rotate_coordinates_e2g!
rotate_vector_e2g!
```

## Backends

```@docs
AbstractBackend
CPUBackend
GPUBackend
```

## IMO and parameter-file functions

```@docs
Imo
query
DataFile
detector_list_from_parameters
run_h_maps_pipeline
```

## Units

```@docs
arcmin2rad
rpm2angfreq
period2rpm
rpm2period
```
