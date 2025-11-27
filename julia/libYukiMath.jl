using Measurements, ForwardDiff, QuadGK, LinearAlgebra

# Derive angle between two vectors.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMathAngleBetweenVectors(vectorA::Vector{Measurement{Float64}}, vectorB::Vector{Measurement{Float64}})::Measurement{Float64}
	return acos(libYukiMathVectorDotProduct(vectorA, vectorB) / (libYukiMathVectorQuantity(vectorA) * libYukiMathVectorQuantity(vectorB)));
end
# function libYukiMathAngleBetweenVectors(vectorA::Vector{Measurement{Float64}}, vectorsB::Vector{Vector{Measurement{Float64}}})::Vector{Measurement{Float64}}
# 	return map(x -> libYukiMathAngleBetweenVectors(vectorA, x), vectorsB);
# end

# Derive unit vector from vector.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMathUnitVector(vectorA::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	return vectorA ./ libYukiMathVectorQuantity(vectorA);
end

# Rotate 3-dimension vector.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMath3DVectorRotate(vectorA::Vector{Measurement{Float64}}, rotatePolarAngle::Measurement{Float64}, rotateAzimuthalAngle::Measurement{Float64})::Vector{Measurement{Float64}}
	vectorMod::Measurement{Float64}, vectorPolarAngle::Measurement{Float64}, vectorAzimuthalAngle::Measurement{Float64} = libYukiMath3DVectorModAndAngle(vectorA);
	return libYukiMath3DVectorFromISO31_11(vectorMod, vectorPolarAngle + rotatePolarAngle, vectorAzimuthalAngle + rotateAzimuthalAngle);
end

# Rotate 2-dimension vector.
# Dependency: Measurement.
# TODO: Validate & Example.
function libYukiMath2DVectorRotate(vectorA::Vector{Measurement{Float64}}, rotateAngle::Measurement{Float64})::Vector{Measurement{Float64}}
	vectorMod::Measurement{Float64}, vectorAngle::Measurement{Float64} = libYukiMath2DVectorModAndAngle(vectorA);
	return libYukiMath2DVectorFromModAngle(vectorMod, vectorAngle * rotateAngle);
end

# Calculate 3-dimension vector mod(r), polar angle(θ), azimuthal angle(φ).
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMath3DVectorModAndAngle(vectorA::Measurement{Float64})::Tuple{Measurement{Float64}, Measurement{Float64}, Measurement{Float64}}
	return sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2) + (vectorA[3] ^ 2)), 
		atan(sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2)) / vectorA[3]),
		atan(vectorA[2] / vectorA[1]);
end

# Calculate 2-dimension vector mod and angle.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMath2DVectorModAndAngle(vectorA::Measurement{Float64})::Tuple{Measurement{Float64}, Measurement{Float64}}
	return sqrt((vectorA[1] ^ 2) + (vectorA[2] ^ 2)), atan(vectorA[2] / vectorA[1]);
end

# Calculate 3-dimension coordinate of vector from (r, θ, φ). ISO 31-11.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMath3DVectorFromISO31_11(radialDistance::Measurement{Float64}, polarAngle::Measurement{Float64}, azimuthalAngle::Measurement{Float64})::Vector{Measurement{Float64}}
	return [radialDistance * sin(polarAngle) * cos(azimuthalAngle), radialDistance * sin(polarAngle) * sin(azimuthalAngle), radialDistance * cos(azimuthalAngle)];
end

# Calculate 2-dimension coordinate of vector from mod and angle.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiMath2DVectorFromModAngle(vectorMod::Measurement{Float64}, vectorAngle::Measurement{Float64})::Vector{Measurement{Float64}}
	return [vectorMod * cos(vectorAngle), vectorMod * sin(vectorAngle)];
end

# Calculate one dimensional integration of function f from limit A to B.
# Dependency: Measurements, QuadGK.
# Example: True.
function libYukiMath1DIntegration(f::Function, limitA::Measurement{Float64}, limitB::Measurement{Float64})::Measurement{Float64}
	result::Measurement{Float64}, error::Measurement{Float64} = QuadGK.quadgk(x -> f(x), limitA, limitB);
	return Measurements.value(result) ± sqrt(Measurements.uncertainty(result) ^ 2 + Measurements.value(error) ^ 2 + Measurements.uncertainty(error) ^ 2);
end

# Calculate forward difference of function f by variable x.
# Dependency: Measurements, ForwardDiff.
# Example: True.
function libYukiMathForwardDifference(f::Function, x::Measurement{Float64})::Measurement{Float64}
	return ForwardDiff.derivative(x -> f(x), x);
end

# Calculate dot product of two vectors.
# Dependency: Measurements, LinearAlgebra.
# TODO: Validate & Example.
function libYukiMathVectorDotProduct(vectorA::Vector{Measurement{Float64}}, vectorB::Vector{Measurement{Float64}})::Measurement{Float64}
	return LinearAlgebra.dot(vectorA, vectorB);
end

# Calculate cross product of two vectors.
# Dependency: Measurements, LinearAlgebra.
# TODO: Validate & Example.
function libYukiMathVectorCrossProduct(vectorA::Vector{Measurement{Float64}}, vectorB::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	return LinearAlgebra.cross(vectorA, vectorB);
end

# Calculate quantity of vector.
# Dependency: Measurements, LinearAlgebra.
# Example: True.
function libYukiMathVectorQuantity(VectorA::Vector{Measurement{Float64}})::Measurement{Float64}
	return LinearAlgebra.norm(VectorA);
end
