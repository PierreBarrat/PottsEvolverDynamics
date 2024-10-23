struct YamlWrapper
    d::Dict{Any,Any}
end

function DrWatson._wsave(filename, data::YamlWrapper, args...; kwargs...)
    return YAML.write_file(filename, data.d)
end

