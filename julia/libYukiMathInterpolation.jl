using Interpolations;

"""
    libYukiMathInterpolationMonotonic(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real}
    )
Perform monotonic interpolation on the given data points `(x, y)`. 
The function returns an interpolation object that can be used to evaluate the interpolated values at any point within the range of 
`x` by calling it with the desired `x` value.
# Arguments
- `x`: A vector of independent variable values.
- `y`: A vector of dependent variable values.
# Returns
- An interpolation object that can be used to evaluate the 
interpolated values.
"""
function libYukiMathInterpolationMonotonic(
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real}
)
    length(x) == length(y) || 
        throw(DimensionMismatch("y and x must have equal lengths."))
    idX = sortperm(x)

    return interpolate(
        x[idX],
        y[idX],
        SteffenMonotonicInterpolation()
    );
end
