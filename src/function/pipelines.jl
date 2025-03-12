function sim_det_scanfields(tomlfile_path::String)
    tomlfile         = TOML.parsefile(tomlfile_path)
    base_path        = tomlfile["simulation"]["base_path"]
    imo_path         = tomlfile["general"]["imo_path"]
    imo              = gen_imo(imo_path)
    telescope        = tomlfile["general"]["telescope"]
    channel_name     = tomlfile["general"]["channel"]
    det_name         = tomlfile["general"]["det_name"]
    inst_info        = get_instrument_info(imo, telescope)
    channel_list     = get_channel_list(imo)
    channel_info     = get_channel_info(imo, channel_name)
    bolonames        = [det_name]

    ss               = gen_ScanningStrategy()
    division         = parse(Int64, tomlfile["simulation"]["division"])
    spin_n           = tomlfile["simulation"]["spin_n"]
    spin_m           = tomlfile["simulation"]["spin_m"]
    coord            = tomlfile["simulation"]["coord"]
    hwp_rpm          = tomlfile["simulation"]["hwp_rpm"]

    ss.nside         = parse(Int64, tomlfile["general"]["nside"])
    ss.sampling_rate = parse(Float64, tomlfile["general"]["sampling_rate"])
    ss.alpha         = parse(Float64, tomlfile["simulation"]["alpha"])
    ss.beta          = parse(Float64, tomlfile["simulation"]["beta"])
    ss.spin_rpm      = parse(Float64, tomlfile["simulation"]["spin_rpm"])

    ss.prec_rpm      = parse(Float64, tomlfile["simulation"]["prec_rpm"])
    ss.duration      = parse(Float64, tomlfile["simulation"]["duration_s"])
    ss.gamma         = parse(Float64, tomlfile["simulation"]["gamma"])

    println("hwp_rpm: ", hwp_rpm)
    if hwp_rpm == "IMO"
        ss.hwp_rpm       = inst_info["hwp_rpm"]
    else
        ss.hwp_rpm       = parse(Float64, hwp_rpm)
    end
    ss.coord         = coord
    if det_name != "boresight"
        imo_name!(ss, imo, name=bolonames)
    end
    show_ss(ss)
    @info ss

    field = get_scanfield(ss, division=division, spin_n=spin_n, spin_m=spin_m)
    if !isdir(base_path)
        mkdir(base_path)
    end
    if !isdir(base_path * "/" * channel_name)
        mkdir(base_path * "/" * channel_name)
    end
    base_path         = base_path * "/" * channel_name
    output_filename   = "/$(det_name).h5"
    create_h5_file(base_path, output_filename, field, savemap=true)
    println("\nSimulation was done successfully.")
    println("The fits file saved: ",  base_path * output_filename)
end


function create_h5_file(base_path::AbstractString, filename::AbstractString, field::scanfield,; savemap=true)
    h5open(base_path * filename, "w") do file
        write(file, "toml/spin_n", Float64.(field.n))
        write(file, "toml/spin_m", Float64.(field.m))
        for f in propertynames(field.ss)
            if f == :quat
                dets_quat = getfield(field.ss, :quat)
                boloname  = getfield(field.ss, :name)
                for i in eachindex(dets_quat)
                    write(file, "ss/$f/$(boloname[i])", dets_quat[i])
                end
            elseif f == :name
                boloname = getfield(field.ss, :name)
                write(file, "ss/$f", boloname)
            elseif f == :info
                nothing
            else
                data = getfield(field.ss, f)
                write(file, "ss/$f", data)
            end
        end
        write(file, "quantify/hitmap_std",   std(field.hitmap))
        write(file, "quantify/n",            Float64.(field.quantify.n))
        write(file, "quantify/m",            Float64.(field.quantify.m))
        write(file, "quantify/mean",         field.quantify.mean)
        write(file, "quantify/std",          field.quantify.std)
        write(file, "quantify/nanmean",      field.quantify.nanmean)
        write(file, "quantify/nanstd",       field.quantify.nanstd)
        write(file, "quantify/nan2one_mean", field.quantify.nan2one_mean)
        write(file, "quantify/nan2one_std",  field.quantify.nan2one_std)
        if savemap == true
            write(file, "hitmap", field.hitmap)
            write(file, "h", field.h)
        end
    end
end
