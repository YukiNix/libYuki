using Measurements, ForwardDiff, QuadGK

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

# Calculate distance between position A to position B (vector A->B).
# Dependency: Measurements.
# Example: True.
function libYukiMathVectorToDistance(positionVectorAB::Vector{Measurement{Float64}})::Measurement{Float64}
    return sqrt(sum(positionVectorAB .^ 2));
end
