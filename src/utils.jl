"""
    regular_time_sampling(t::AbstractVector{T}, atol::T=0) where T

    This function checks if the time sampling is regular or not.
    It first checks exactly if the time sampling is regular or not and then
    also checks if the time sampling is regular within a certain tolerance.

    # Arguments
    - `t::AbstractVector{T}` : time vector
    - `atol::T=0` : Absolute tolerance for the time sampling for `isapprox` function

    # Returns
    - `Bool` : True if the time sampling is regular, False otherwise
"""
function regular_time_sampling(t::AbstractVector{T}; atol = 0.0) where {T}
    dt = diff(t)
    if length(unique(dt)) == 1  || all(isapprox.(dt, dt[1], atol = atol))
        return true
    end
    return false
end
