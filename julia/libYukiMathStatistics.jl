using MultivariateStats, Statistics;

"""
    libYukiMathIsotonicRegression(
        x, 
        y;
        weights = nothing,
        isDecreasing = false
    )
Perform isotonic regression on the given data points `(x, y)` with optional weights. The function returns the estimated values of `y` that are monotonically increasing (or decreasing if specified) with respect to `x`.
# Arguments
- `x`: A vector of independent variable values.
- `y`: A vector of dependent variable values.
- `weights`: An optional vector of weights for the data points. If not provided, equal weights are assumed.
- `isDecreasing`: A boolean flag indicating whether to perform decreasing isotonic regression. Default is `false` (increasing).
# Returns
- A vector of estimated `y` values that are monotonically increasing (or decreasing) with respect to `x`.
"""
function libYukiMathIsotonicRegression(
	x::AbstractVector{<:Real}, 
	y::AbstractVector{<:Real};
	weights::Union{AbstractVector{<:Real}, Nothing} = nothing,
    isDecreasing::Bool = false
)
	length(x) == length(y) || 
		throw(DimensionMismatch("y and x must have equal lengths."))
	weights === nothing || 
		length(weights) == length(x) || 
		throw(DimensionMismatch("weights must have equal length as y and x."))
    sortedY = isDecreasing ? -y : y;
    idX = sortperm(x);
    sortedX = x[idX];
    sortedY = sortedY[idX];
    sortedW = isnothing(weights) ? nothing : weights[idX];

    estimateY = isnothing(sortedW) ? MultivariateStats.isotonic(sortedX, sortedY) : MultivariateStats.isotonic(sortedX, sortedY, sortedW);

	return isDecreasing ? -estimateY[invperm(idX)] : estimateY[invperm(idX)];
end
