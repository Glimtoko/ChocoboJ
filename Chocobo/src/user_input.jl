using Parameters
using IniFile

# Custom types
@enum ANTI_HG_TYPES begin
    ANTI_HG_NONE = 1
    ANTI_HG_DYNA_1 = 2
    ANTI_HG_DYNA_2 = 3
    ANTI_HG_BF = 4
end

@with_kw mutable struct Input
    material_list = Dict(1 => 1, 2 => 2)

    # Artificial Viscosity
    qform = 2
    cl = 0.5
    cq = 0.75

    # Robustness
    anti_hg = [ANTI_HG_NONE, ANTI_HG_NONE]
    kappa = 0.7

    # Mesh Extents
    x0 = 0.0
    x1 = 1.0
    y0 = 0.0
    y1 = 1.0
    nregions = 1
    meshx = 50
    meshy = 50

    # Timestep, etc.
    t0 = 0.0
    tend = 0.205
    dtinit = 0.0001
    dtmax = 0.01
    growth = 1.02
    dtmin_hg = 1.0e-6

    # Output frequency
    dtsilo = 0.01

    # Profiling/Debug settings
    debug_step_count = 0   # Set to zero to disable debugging
    profile = false
end

function read_input(input::Input, filename::String)
    core_sections = ["Q", "SPHSOD MESH", "CONTROL", "DEBUG", "ROBUST"]

    ini = read(Inifile(), filename)
    s = sections(ini)
    for k in keys(s)
        if uppercase(k) in core_sections
            section = s[k]
            for v in keys(section)
                vlc = lowercase(v)
                vs = Symbol(vlc)
                if hasproperty(input, vs)
                    d = getproperty(input, vs)
                    t = typeof(d)
                    value = get(ini, k, v, d)
                    value = try
                        parse(t, value)
                    catch
                        data = split(value, ",")
                        value_arr = []
                        for d in data
                            val = parse(Int32, d)
                            push!(value_arr, val)
                        end
                        value_arr
                    end
                    setproperty!(input, vs, value)
                else
                    error("Unexpected variable $v in section $k")
                end
            end
        elseif uppercase(k) == "MATERIAL"
            section = s[k]
            for v in keys(section)
                if uppercase(v) == "MATERIAL_NUMBERS"
                    data = split(get(ini, k, v), ",")
                    index = 1
                    input.material_list = Dict()
                    for d in data
                        input.material_list[index] = parse(Int32, d)
                        index += 1
                    end
                else
                    error("Unexpected variable $v in section $k")
                end
            end
        else
            error("Unexpected section $k")
        end
    end

    # Sort out the anti_hg enum
    input.anti_hg = try
        map(x -> ANTI_HG_TYPES(x), input.anti_hg)
    catch
        input.anti_hg
    end
end
