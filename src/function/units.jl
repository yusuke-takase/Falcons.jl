function arcmin2rad(x)
    return deg2rad(x/60.)
end

function rpm2angfreq(rpm)
   return (2.0π / 60.0) * rpm
end

function period2rpm(period,; unit="min")
    if unit == "min"
        rpm = 1.0 / period
    end
    if unit == "sec"
        rpm = 1.0 / (period/60.0)
    end
    if unit == "hour"
        rpm = 1.0 / (period*60.0)
    end
    return rpm
end

function rpm2period(rpm,; unit="min")
    if unit == "min"
        period = 1.0 / rpm
    end
    if unit == "sec"
        period = 1.0 / (rpm/60.0)
    end
    if unit == "hour"
        period = 1.0 / (rpm/60.0/60.0)
    end
    return period
end
