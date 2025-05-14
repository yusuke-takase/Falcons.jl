struct Imo
    imo::Dict
end

function gen_imo(path)
    imo = JSON.parsefile(path)
    return Imo(imo)
end

function get_channel_list(imo::Imo)
    channel_info = []
    for i in eachindex(imo.imo["data_files"])
        if imo.imo["data_files"][i]["name"] == "channel_info"
            ch = imo.imo["data_files"][i]["metadata"]["channel"]
            push!(channel_info, ch)
        end
    end
    return reverse(channel_info)
end

function get_channel_info(imo::Imo)
    df = DataFrame()  # 空のDataFrameを初期化

    for i in eachindex(imo.imo["data_files"])
        if imo.imo["data_files"][i]["name"] == "channel_info"
            metadata = imo.imo["data_files"][i]["metadata"]
            metadata_keys = keys(metadata)

            keys_to_exclude = Set(["detector_objs", "detector_names"])
            filtered_metadata = Dict(
                key => value for (key, value) in metadata if !(key in keys_to_exclude)
            )

            if isempty(df)
                df = DataFrame(permutedims(collect(values(filtered_metadata))), collect(keys(filtered_metadata)))
            else
                push!(df, collect(values(filtered_metadata)))
            end
        end
    end
    sort!(df, :bandcenter_ghz)
    return df
end

function get_channel_info(imo::Imo, channel_name::String)
    df = get_channel_info(imo)
    channel_info = filter(row -> row.channel == channel_name, df)
    return channel_info
end

function imo_channel!(ss::ScanningStrategy, imo::Imo,; channel)
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
        #println("The channel `$(channel)` is set from IMo.")
    end
    return ss
end


function imo_channel!(ss::ScanningStrategy, imo::Imo, reformation::Bool, ; channel,)
    switch = 0
    df = ""
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixel") == true
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
        #println("The channel `$(channel)` is set from IMo.")
    end
    return ss
end

function imo_telescope!(ss::ScanningStrategy, imo::Imo,;telescope)
    switch = 0
    df = 0
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixtype") == true
            boloname = imo.imo["data_files"][i]["metadata"]["name"]
            teles = split.(imo.imo["data_files"][i]["metadata"]["channel"], "")[1] * "FT"
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
        #println("The telescope `$(telescope)FT' is set from IMo.")
    end
    return ss
end

function imo_name!(ss::ScanningStrategy, imo::Imo,;name::Vector)
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
        #println("The detector")
        #for i in eachindex(ss.name)
        #    println("     `$(ss.name[i])`")
        #end
        #println("is set from IMo.")
    end
    return ss
end


function imo_name!(ss::ScanningStrategy, imo::Imo, reformation::Bool, ;name::Vector)
    switch = 0
    df = 0
    for i in eachindex(imo.imo["data_files"])
        if haskey(imo.imo["data_files"][i]["metadata"], "pixel") == true
            boloname = imo.imo["data_files"][i]["metadata"]["name"]
            #println(boloname)
            #teles = split.(imo.imo["data_files"][i]["metadata"]["channel"], "")[1]
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
        #println("The detector")
        #for i in eachindex(ss.name)
        #    println("     `$(ss.name[i])`")
        #end
        #println("is set from IMo.")
    end
    return ss
end

function get_detectors(path)
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

function get_pol_angle(ss::ScanningStrategy, i::Int)
    if size(ss.info) == (0,0)
        # boresight
        return 0.
    end
    if hasproperty(ss.info, :orient)
        try
            polang = deg2rad(parse(Float64, ss.info.orient[i]))
            return polang
        catch
            if ss.info.orient[i] == "Q"
                if ss.info.pol[i] == "T"
                    polang = 0
                elseif ss.info.pol[i] == "B"
                    polang = π/2
                end
            end
            if ss.info.orient[i] == "U"
                if ss.info.pol[i] == "T"
                    polang = π/4
                elseif ss.info.pol[i] == "B"
                    polang = 3π/4
                end
            end
            return polang
        end
    else
        return ss.info.pol_angle_rad[i]
    end
end
