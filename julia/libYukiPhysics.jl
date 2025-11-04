using Measurements, ForwardDiff

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

# TODO: Adfinitas

# Derive acceleration from force. The second Newton.
# Dependency: Measurements.
function libYukiPhysicsForceToAcceleration(mass::Measurement{Float64}, force::Measurement{Float64})
    return force ./ mass;
end

# Derive force from potential. vector(source -> object).
# Dependency: Measurements.
function libYukiPhysicsPotentialForce(potentialFunction::Function, objectPosition::Vector{Measurement{Float64}}, differenceStep::Measurement{Float64})::Vector{Measurement{Float64}}
    positionVectorMod::Measurement{Float64} = libYukiMathVectorToDistance(objectPosition);
    return libYukiMathForwardDifference(x -> potentialFunction(x), positionVectorMod) .* (objectPosition ./ positionVectorMod);
end

# Calculate gravitational potential between source and motion object. vector(source -> object).
# Dependency: Measurements.
function libYukiPhysicsGravitaionalPotential(gravitationalConstant::Measurement{Float64}, sourceMass::Measurement{Float64}, objectPosition::Vector{Measurement{Float64}})::Measurement{Float64}
    return (gravitationalConstant * sourceMass) / libYukiMathVectorToDistance(objectPosition);
end
