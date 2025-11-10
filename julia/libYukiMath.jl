using Measurements, ForwardDiff, QuadGK

# Calculate one dimensional integration of function f from limit A to B.
# Dependency: Measurements, QuadGK.
function libYukiMathOneDimensionalIntegration(f::Function, limitA::Measurement{Float64}, limitB::Measurement{Float64})::Measurement{Float64}
    result::Measurement{Float64}, error::Measurement{Float64} = QuadGK.quadgk(x -> f(x), limitA, limitB);
    return Measurements.value(result) ± sqrt(Measurements.uncertainty(result) ^ 2 + Measurements.value(error) ^ 2 + Measurements.uncertainty(error) ^ 2);
end

# Calculate forward difference of function f by variable x.
# Dependency: Measurements, ForwardDiff.
function libYukiMathForwardDifference(f::Function, x::Measurement{Float64})::Measurement{Float64}
    return ForwardDiff.derivative(x -> f(x), x);
end

# Calculate distance between position A to position B (vector A->B).
# Dependency: Measurements.
function libYukiMathVectorToDistance(positionVectorAB::Vector{Measurement{Float64}})::Measurement{Float64}
    return sqrt(sum(positionVectorAB .^ 2));
end
