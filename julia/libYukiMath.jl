using Measurements, ForwardDiff, QuadGK, LinearAlgebra

# Derive angle between two vectors, return [0, π].
# Example: True.
function libYukiMathAngleBetweenVector(vectorA::AbstractVector{<:Real}, vectorB::AbstractVector{<:Real})::Real 
	return acos(libYukiMathVectorDotProduct(vectorA, vectorB) / (libYukiMathVectorQuantity(vectorA) * libYukiMathVectorQuantity(vectorB)));
end

# Derive unit vector from vector.
# Example: True.
function libYukiMathUnitVector(vectorA::AbstractVector{<:Real})::AbstractVector{<:Real} 
	return vectorA ./ libYukiMathVectorQuantity(vectorA);
end

# Rotate 3-dimension vector.
# Example: True.
function libYukiMath3DVectorRotate(vectorA::AbstractVector{<:Real}, rotatePolarAngle::Real, rotateAzimuthalAngle::Real)::AbstractVector{<:Real} 
	vectorMod, vectorPolarAngle, vectorAzimuthalAngle = libYukiMath3DVectorModAndAngle(vectorA);
	return libYukiMath3DVectorFromISO31_11(vectorMod, vectorPolarAngle + rotatePolarAngle, vectorAzimuthalAngle + rotateAzimuthalAngle);
end

# Rotate 2-dimension vector.
# Example: True.
function libYukiMath2DVectorRotate(vectorA::AbstractVector{<:Real}, rotateAngle::Real)::AbstractVector{<:Real} 
	vectorMod, vectorAngle = libYukiMath2DVectorModAndAngle(vectorA);
	return libYukiMath2DVectorFromModAngle(vectorMod, vectorAngle + rotateAngle);
end

# Calculate 3-dimension vector mod(r), polar angle(θ), azimuthal angle(φ).
# Example: True.
function libYukiMath3DVectorModAndAngle(vectorA::AbstractVector{<:Real})::Tuple{Real, Real, Real} 
	return sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2) + (vectorA[3] ^ 2)), 
		atan(sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2)) / vectorA[3]),
		atan(vectorA[2] / vectorA[1]);
end

# Calculate 2-dimension vector mod and angle.
# Example: True.
function libYukiMath2DVectorModAndAngle(vectorA::AbstractVector{<:Real})::Tuple{Real, Real} 
	return sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2)), atan(vectorA[2] / vectorA[1]);
end

# Calculate 3-dimension coordinate of vector from (r, θ, φ). ISO 31-11.
# Example: True.
function libYukiMath3DVectorFromISO31_11(radialDistance::Real, polarAngle::Real, azimuthalAngle::Real)::AbstractVector{<:Real} 
	return [radialDistance * sin(polarAngle) * cos(azimuthalAngle), radialDistance * sin(polarAngle) * sin(azimuthalAngle), radialDistance * cos(azimuthalAngle)];
end

# Calculate 2-dimension coordinate of vector from mod and angle.
# Example: True.
function libYukiMath2DVectorFromModAngle(vectorMod::Real, vectorAngle::Real)::AbstractVector{<:Real} 
	return [vectorMod * cos(vectorAngle), vectorMod * sin(vectorAngle)];
end

# Calculate one dimensional integration of function f from limit A to B.
# Dependency: QuadGK.
# Example: True.
function libYukiMath1DIntegration(f::Function, limitA::Real, limitB::Real)::Real 
	result, _ = QuadGK.quadgk(x -> f(x), limitA, limitB);
	return result;
end

# Calculate one dimensional integration of function f from limit A to B (with uncertainty propagation).
# Dependency: Measurements, QuadGK.
# Example: True.
function libYukiMath1DIntegration(f::Function, limitA::Real, limitB::Real)::Real
	result, error = QuadGK.quadgk(x -> f(x), limitA, limitB);
	return Measurements.value(result) ± sqrt(Measurements.uncertainty(result) ^ 2 + Measurements.value(error) ^ 2 + Measurements.uncertainty(error) ^ 2);
end

# Calculate forward difference of function f by variable x.
# Dependency: ForwardDiff.
# Example: True.
function libYukiMathForwardDifference(f::Function, x::Real)::Real 
	return ForwardDiff.derivative(x -> f(x), x);
end

# Calculate dot product of two vectors.
# Dependency: LinearAlgebra.
# Example: True.
function libYukiMathVectorDotProduct(vectorA::AbstractVector{<:Real}, vectorB::AbstractVector{<:Real})::Real 
	return LinearAlgebra.dot(vectorA, vectorB);
end

# Calculate cross product of two vectors.
# Dependency: LinearAlgebra.
# Example: True.
function libYukiMathVectorCrossProduct(vectorA::AbstractVector{<:Real}, vectorB::AbstractVector{<:Real})::AbstractVector{<:Real} 
	return LinearAlgebra.cross(vectorA, vectorB);
end

# Calculate quantity of vector.
# Dependency: LinearAlgebra.
# Example: True.
function libYukiMathVectorQuantity(VectorA::AbstractVector{<:Real})::Real 
	return LinearAlgebra.norm(VectorA);
end
