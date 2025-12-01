#______________________________________________________
# ADSIM Version Module
# Provides centralized version management
# Version: 0.1.0
#______________________________________________________

module ADSIMVersion

"""
    get_version()

Read and return the version string from the VERSION file at project root.
Returns the version as a string (e.g., "0.1.0").
"""
function get_version()
    version_file = joinpath(@__DIR__, "..", "VERSION")
    if isfile(version_file)
        return strip(read(version_file, String))
    else
        @warn "VERSION file not found, returning default version"
        return "0.0.0-dev"
    end
end

"""
    get_version_string()

Return a formatted version string suitable for display.
Returns string in format "ADSIM v0.1.0".
"""
function get_version_string()
    return "ADSIM v$(get_version())"
end

export get_version, get_version_string

end # module
