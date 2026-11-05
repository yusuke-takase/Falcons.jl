_bad_unit(unit) = throw(ArgumentError("unknown unit \"$unit\"; use \"min\", \"sec\" or \"hour\""))

"""
    arcmin2rad(x) -> Float64

Convert an angle from arcminutes to radians.
"""
arcmin2rad(x) = deg2rad(x / 60.0)

"""
    rpm2angfreq(rpm) -> Float64

Convert a rotation rate from revolutions per minute to an angular frequency
(rad/s).
"""
rpm2angfreq(rpm) = (2.0π / 60.0) * rpm

"""
    period2rpm(period; unit = "min") -> Float64

Rotation rate (rpm) of a rotation with the given `period`, expressed in
`unit` (`"min"`, `"sec"` or `"hour"`).
"""
function period2rpm(period; unit::AbstractString = "min")
	unit == "min" && return 1.0 / period
	unit == "sec" && return 1.0 / (period / 60.0)
	unit == "hour" && return 1.0 / (period * 60.0)
	_bad_unit(unit)
end

"""
    rpm2period(rpm; unit = "min") -> Float64

Period of a rotation at `rpm` revolutions per minute, expressed in `unit`
(`"min"`, `"sec"` or `"hour"`).
"""
function rpm2period(rpm; unit::AbstractString = "min")
	unit == "min" && return 1.0 / rpm
	unit == "sec" && return 1.0 / (rpm / 60.0)
	unit == "hour" && return 1.0 / (rpm * 60.0)
	_bad_unit(unit)
end
