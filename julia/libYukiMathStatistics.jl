using MultivariateStats, Statistics, SmoothingSplines;

"""
    libYukiMathStatisticsSmoothingSplines(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real};
        lambda::Real = 1e-4
    )
Perform smoothing splines on the given data points `(x, y)` with
a specified smoothing parameter `lambda`. The function returns the
estimated values of `y` based on the smoothing splines fit.
# Arguments
- `x`: A vector of independent variable values.
- `y`: A vector of dependent variable values.
- `lambda`: A smoothing parameter that controls the trade-off between
the smoothness of the fitted curve and the closeness to the data points.
# Returns
- A vector of estimated `y` values based on the smoothing splines fit.
"""
function libYukiMathStatisticsSmoothingSplines(
	x::AbstractVector{<:Real},
	y::AbstractVector{<:Real};
	lambda::Real = 1e-4,
)
	idx = sortperm(x);
	spline = SmoothingSplines.fit(
		SmoothingSpline,
		x[idx],
		y[idx],
		lambda,
	);
	return predict(spline, x[idx]);
end

"""
    libYukiMathStatisticsIsotonicRegression(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real};
        weights::Union{AbstractVector{<:Real}, Nothing} = nothing,
        isDecreasing::Bool = false
    )
Perform isotonic regression on the given data points `(x, y)` with 
optional weights. The function returns the estimated values of `y` 
that are monotonically increasing (or decreasing if specified) with 
respect to `x`.
# Arguments
- `x`: A vector of independent variable values.
- `y`: A vector of dependent variable values.
- `weights`: An optional vector of weights for the data points. 
If not provided, equal weights are assumed.
- `isDecreasing`: A boolean flag indicating whether to perform 
decreasing isotonic regression. Default is `false` (increasing).
# Returns
- A vector of estimated `y` values that are monotonically 
increasing (or decreasing) with respect to `x`.
"""
function libYukiMathStatisticsIsotonicRegression(
	x::AbstractVector{<:Real}, 
	y::AbstractVector{<:Real};
	weights::Union{AbstractVector{<:Real}, Nothing} = nothing,
    isDecreasing::Bool = false
)
	length(x) == length(y) || 
		throw(DimensionMismatch("y and x must have equal lengths."))
	weights === nothing || 
		length(weights) == length(x) || 
		    throw(
                DimensionMismatch(
                    "weights must have equal length as y and x."
            ))
    sortedY = 
    idX = sortperm(x);
    sortedX = x[idX];
    sortedY = isDecreasing ? -y[idX] : y[idX];
    sortedW = isnothing(weights) ? nothing : weights[idX];

    estimateY = isnothing(sortedW) ? 
        MultivariateStats.isotonic(sortedX, sortedY) : 
        MultivariateStats.isotonic(sortedX, sortedY, sortedW);

	return isDecreasing ?
         -estimateY[invperm(idX)] : 
         estimateY[invperm(idX)];
end
