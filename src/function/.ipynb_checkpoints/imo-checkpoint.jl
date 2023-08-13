struct Imo
    imo::Dict
end

function gen_imo(path)
    imo = JSON.parsefile(path)
    return Imo(imo)
end

function get_channel_info(imo::Imo)
    channel_info = []
    for i in eachindex(imo.imo["data_files"])
        if imo.imo["data_files"][i]["name"] == "channel_info"
            ch = imo.imo["data_files"][i]["metadata"]["channel"]
            push!(channel_info, ch)
        end
    end
    return reverse(channel_info)
end

function imo_channel!(ss::ScanningStrategy_imo, imo::Imo,; channel)
    switch = 0
    df = ""
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixtype") == true
            boloname = imo.imo["data_files"][i]["metadata"]["name"]
            teles = split.(imo.imo["data_files"][i]["metadata"]["channel"], "")[1]
            ch = imo.imo["data_files"][i]["metadata"]["channel"]
            if ch == channel
                metadata = imo.imo["data_files"][i]["metadata"]
                if switch == 0
                    df = DataFrame(permutedims(collect(values(metadata))), collect(keys(metadata)))
                    switch = 1
                else
                    push!(df, collect(values(metadata)))
                end
            end
        end
    end
    if df == ""
         @error "No such channel in the IMo."
    else
        ss.quat = Vector{Vector{Float64}}(df.quat)
        ss.name = df.name
        ss.info = df
        println("The channel `$(channel)` is set from IMo.")
    end
    return ss
end

function imo_telescope!(ss::ScanningStrategy_imo, imo::Imo,;telescope)
    switch = 0
    df = 0
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixtype") == true
            boloname = imo.imo["data_files"][i]["metadata"]["name"]
            teles = split.(imo.imo["data_files"][i]["metadata"]["channel"], "")[1]
            ch = imo.imo["data_files"][i]["metadata"]["channel"]
            if teles == telescope
                metadata = imo.imo["data_files"][i]["metadata"]
                if switch == 0
                    df = DataFrame(permutedims(collect(values(metadata))), collect(keys(metadata)))
                    switch = 1
                else
                    push!(df, collect(values(metadata)))
                end
            end
        end
    end
    if df == ""
         @error "No such channel in the IMo."
    else
        ss.quat = Vector{Vector{Float64}}(df.quat)
        ss.name = df.name
        ss.info = df
        println("The telescope `$(telescope)FT' is set from IMo.")
    end
    return ss
end

function imo_name!(ss::ScanningStrategy_imo, imo::Imo,;name::Vector)
    switch = 0
    df = 0
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixtype") == true
            boloname = imo.imo["data_files"][i]["metadata"]["name"]
            teles = split.(imo.imo["data_files"][i]["metadata"]["channel"], "")[1]
            ch = imo.imo["data_files"][i]["metadata"]["channel"]
            for j in eachindex(name)
                if boloname == name[j]
                    metadata = imo.imo["data_files"][i]["metadata"]
                    if switch == 0
                        df = DataFrame(permutedims(collect(values(metadata))), collect(keys(metadata)))
                        switch = 1
                    else
                        push!(df, collect(values(metadata)))
                    end
                end
            end
        end
    end
    if df == ""
         @error "No such detector in the IMo."
    else
        ss.quat = Vector{Vector{Float64}}(df.quat)
        ss.name = df.name
        ss.info = df
        println("The detector `$(name)` is set from IMo.")
    end
    return ss
end

function get_det(path)
    bolonames = []
    open(path) do file
        for (index, line) in enumerate(eachline(file))
            row = split(line, '\t')
            filter!(x -> x ≠ "", row)
            if index == 1
                header = split(line)[2:end]
            else
                push!(bolonames, row[6])
            end
        end
    end
    return bolonames
end

function get_instrument_info(imo::Imo, inst)
    instrument_info = ""
    for i in eachindex(imo.imo["data_files"])
        if imo.imo["data_files"][i]["name"] == "instrument_info"
            if imo.imo["data_files"][i]["metadata"]["name"] == inst
                instrument_info = imo.imo["data_files"][i]["metadata"]
            end
        end
    end
    return instrument_info
end