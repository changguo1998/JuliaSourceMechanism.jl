"""
module ArgParser

a simple module to help parse commandline parameters

usage:

 1. use `newoptionlist()` to get an empty option list

```julia
opts = ArgParser.newoptionlist()
```

 2. add options using `addopt!()``

```julia
ArgParser.addopt!(opts, "fullname", "defaultvalue", "abbrevationname")
...
```

 3. parse parameters

```julia
(named_parameter, positional_parameter) = ArgParser.parse(ARGS, opts)
```
"""
module ArgParser
struct Option
    key::String
    value::String
    abbr::String
    function Option(key::String, val::String = "", abbr::String = "")
        return new(key, val, abbr)
    end
end

function addopt!(l::Vector{Option}, x::Option)
    push!(l, x)
    return nothing
end

function addopt!(l::Vector{Option}, key::String, val::String = "", abbr::String = "")
    push!(l, Option(key, val, abbr))
    return nothing
end

function newoptionlist()
    return Vector{Option}(undef, 0)
end

"""
"""
function parse(args::Vector{String}, opts::Vector{Option})
    pars = Dict{String,String}()
    par_processed = falses(length(args))
    for o in opts
        flag = args .== ("--" * o.key)
        if !isempty(o.abbr)
            flag .|= args .== ("-" * o.abbr)
        end
        idx = findall(flag)
        if isempty(idx)
            pars[o.key] = o.value
        elseif length(idx) == 1
            pars[o.key] = args[idx[1]+1]
            par_processed[idx[1]] = true
            par_processed[idx[1]+1] = true
        else
            @error "More than one arguments match the same option."
        end
    end
    return (pars, args[.!par_processed])
end
end
