using Measurements, StaticArrays;
using ForwardDiff, QuadGK, LinearAlgebra;

include("libYukiBasic.jl")
include("libYukiMathVector.jl")

"""
    libYukiMathAngleBetweenVectors(
		vectorA, 
		vectorB
	)
Return the angle between two real vectors in radians.
# Arguments
- `vectorA::AbstractVector{<:Real}`: The first vector.
- `vectorB::AbstractVector{<:Real}`: The second vector.
# Returns
- The angle between `vectorA` and `vectorB`, in radians 
`[0, π]`.
# Notes
The value passed to `acos` is clamped to `[-1, 1]` to 
avoid floating-point round-off errors.
# Example
```julia
	vecA = [
		1 ± 0.01, 
		2 ± 0.02, 
		3 ± 0.03
	];
	vecB = [
		4 ± 0.04, 
		5 ± 0.05, 
		6 ± 0.06
	];
	libYukiMathAngleBetweenVectors(
		vecA, 
		vecB
	)
	# Output: 0.2257 ± 0.007
```
"""
libYukiMathAngleBetweenVectors(
	vectorA::AbstractVector{<:Real},
	vectorB::AbstractVector{<:Real}
)::Real = 
	(
		libYukiMathVectorDotProduct(vectorA, vectorB) /
			((vectorA |> libYukiMathVectorQuantity) * 
				(vectorB |> libYukiMathVectorQuantity))
	) |> acos;

"""
	libYukiMathIntegration(
		f, 
		limits...; 
		kwargs...
	)
Calculate the integral of a function `f` over specified 
	limits using adaptive quadrature.
# Arguments
- `f`: The function to be integrated. It can be a function 
of one or more variables.
- `limits...`: A variable number of tuples, each 
containing the lower and upper limits for integration. 
	For example, `(a, b)` for a single variable, or 
	`(a1, b1), (a2, b2), ...` for multiple variables.
- `kwargs...`: Additional keyword arguments to be passed 
to the `QuadGK.quadgk` function, such as `rtol` for relative 
	tolerance or `maxevals` for maximum evaluations.
# Returns
- The integral of the function `f` over the specified 
limits, along with an estimate of the uncertainty in the 
result. The result is returned as a `Measurement` object, 
which includes both the value of the integral and its 
uncertainty.
# Notes
- The function uses adaptive quadrature to efficiently 
compute the integral, and it can handle both single and 
multiple integrals by recursively calling itself for nested 
	integrations.
- The uncertainty in the result is calculated by combining 
the uncertainty from the integral result and the estimated 
error from the quadrature method using the `hypot` function 
	to ensure proper propagation of uncertainties.
# Example
```julia
	f(x) = x ^ 2;
	f(x, y) = x ^ 2 + y ^ 2;
	xLimit = 
		0 ± 0.01, 
		1 ± 0.05;
	yLimit = 0, 1;
	libYukiMathIntegration(
		f, 
		xLimit
	)
	# Output: 0.333 ± 0.05
	libYukiMathIntegration(
		f, 
		xLimit, 
		yLimit
	)
	# Output: 0.667 ± 0.038
```
"""
function libYukiMathIntegration(
	f, 
	limits...; 
	kwargs...
)
    return _libYukiMathIntegration(
		f, 
		(), 
		limits...; 
		kwargs...
	)
end
function _libYukiMathIntegration(
	f, 
	outerVariables, 
	limit; 
	kwargs...
)
    limitA, limitB = 
		_libYukiMathIntegrationEvaluateLimit(
			limit, 
			outerVariables...
		);
    result, error = QuadGK.quadgk(
        x -> f(outerVariables..., x),
        limitA,
        limitB;
        kwargs...
    )
    return Measurements.value(result) ± hypot(
        Measurements.uncertainty(result),
        Measurements.value(error),
        Measurements.uncertainty(error),
    )
end
function _libYukiMathIntegration(
	f, 
	outerVariables, 
	limit, 
	restLimits...; 
	kwargs...
)
    limitA, limitB = 
		_libYukiMathIntegrationEvaluateLimit(
			limit, 
			outerVariables...
		);
    result, error = QuadGK.quadgk(
        x -> _libYukiMathIntegration(
            f,
            (outerVariables..., x),
            restLimits...;
            kwargs...
        ),
        limitA,
        limitB;
        kwargs...
    );
    return Measurements.value(result) ± 
		hypot(
			Measurements.uncertainty(result),
			Measurements.value(error),
			Measurements.uncertainty(error),
    	);
end
function _libYukiMathIntegrationEvaluateLimit(
	limit, 
	variables...
)
    limitA, limitB = limit;
    return _libYukiMathIntegrationEvaluateBound(
			limitA, 
			variables...
		),
		_libYukiMathIntegrationEvaluateBound(
			limitB, 
			variables...
	);
end
_libYukiMathIntegrationEvaluateBound(bound, variables...) =
    bound isa Function ? bound(variables...) : bound;

"""
	libYukiMathDerivative(
		f, 
		x
	)
	libYukiMathJacobian(
		f, 
		x
	)
	libYukiMathGradient(
		f, 
		x
	)
	
Calculate the first-order derivative of a function `f` with respect to 
variable `x`.
# Arguments
- `f`: The function for which to calculate the derivative.
- `x`: The point at which to evaluate the derivative.
# Returns
- The derivative of `f` with respect to `x` evaluated at the 
point `x`.
# Notes
- This function is efficient and accurate for computing 
derivatives of functions.
# Example
```julia
	f(x) = x ^ 2;
	libYukiMathDerivative(
		f, 
		2 ± 0.1
	)
	# Output: 4.0 ± 0.2
```
"""
libYukiMathDerivative(f, x::Real) =
    ForwardDiff.derivative(f, x);
libYukiMathJacobian(f, x::AbstractVector{<:Real}) =
    ForwardDiff.jacobian(f, x);
libYukiMathGradient(f, x::AbstractVector{<:Real}) =
    ForwardDiff.gradient(f, x);
function libYukiMathDerivative(f, x::AbstractVector{<:Real})
	y = f(x);
    return (y isa Real) ? libYukiMathGradient(f, x) : (
		(y isa AbstractVector) ? libYukiMathJacobian(f, x) : throw(ArgumentError("f(x) must return a Real scalar or an AbstractVector"))
	);
end

"""
	libYukiMathHessian(
		f,
		x
	)
Calculate the second-order derivative (Hessian) of a function `f` with respect to 
the vector of variables `x`.
# Arguments
- `f`: The function for which to calculate the Hessian.
- `x`: The point at which to evaluate the Hessian.
# Returns
- The Hessian matrix of `f` with respect to `x` evaluated at the 
point `x`.
# Notes
- This function is efficient and accurate for computing 
second-order derivatives of multivariate functions.
# Example
```julia
	f(x) = x[1]^2 + x[2]^2;
	libYukiMathHessian(
		f, 
		[1.0, 2.0]
	)
	# Output: [2.0 0.0; 0.0 2.0]
```
"""
libYukiMathHessian(f, x::AbstractVector{<:Real}) =
    ForwardDiff.hessian(f, x);
