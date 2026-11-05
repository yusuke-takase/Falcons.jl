"""
    IdealHWP(; ang_speed_radpsec, start_angle_rad = 0.0)

Ideal half-wave plate rotating at `ang_speed_radpsec` (rad/s), as `IdealHWP` in
litebird_sim. Its mechanical angle at the time `t` is

    χ(t) = start_angle_rad + ang_speed_radpsec · (t − t_sim)

wrapped to [−π, π], where `t_sim` is the `start_time_simulation_s` of the
[`Observation`](@ref). The angle is evaluated in closed form, so it does not
accumulate rounding errors over long observations.
"""
Base.@kwdef struct IdealHWP
	ang_speed_radpsec::Float64
	start_angle_rad::Float64 = 0.0
end

"""
    hwp_angle(hwp::IdealHWP, time_from_sim_start_s) -> χ

HWP angle (rad, in [−π, π]) at `time_from_sim_start_s` seconds after the start
of the simulation.
"""
@inline hwp_angle(hwp::IdealHWP, time_from_sim_start_s::Real) =
	rem2pi(hwp.start_angle_rad + hwp.ang_speed_radpsec * time_from_sim_start_s, RoundNearest)
