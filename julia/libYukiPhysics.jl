using Measurements, ForwardDiff, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

# Physical object data structure
# Dependency: Measurements.
# Example: True.
mutable struct libYukiPhysicsBody
	name::String
	position::Vector{Vector{Measurement{Float64}}}
	velocity::Vector{Vector{Measurement{Float64}}}
	mass::Measurement{Float64}
	charge::Measurement{Float64}
	radius::Measurement{Float64}
	libYukiPhysicsBody(name, position, velocity, mass, charge, radius) = new(name, [position], [velocity], mass, charge, radius)
end

# Derive gravitational circular motion's orbit radius from source mass and angular velocity quantity.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassAngularVelocity(sourceMass::Measurement{Float64}, angularVelocityQuantity::Measurement{Float64}, gravitationalConstant::Measurement{Float64})::Measurement{Float64}
	return (gravitationalConstant * sourceMass / (angularVelocityQuantity ^ 2)) ^ (1 / 3.);
end

# Derive gravitational circular motion's velocity quantity from source mass and orbit radius.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionVelocityQuantity(sourceMass::Measurement{Float64}, orbitRadius::Measurement{Float64}, gravitationalConstant::Measurement{Float64})::Measurement{Float64}
	return sqrt(gravitationalConstant * sourceMass / orbitRadius);
end

# Derive circular motion's velocity from acceleration and angular velocity.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAccelerationAngularVelocity(acceleration::Vector{Measurement{Float64}}, angularVelocity::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	return libYukiMathUnitVector(cross(acceleration, angularVelocity)) .* (libYukiMathVectorQuantity(acceleration) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's orbit radius from velocity and angular velocity.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionOrbitDisplacementFromVelocityAngularVelocity(velocity::Vector{Measurement{Float64}}, angularVelocity::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	return libYukiMathUnitVector(cross(velocity, angularVelocity)) .* (libYukiMathVectorQuantity(velocity) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's angular velocity to velocity.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAngularVelocity(angularVelocity::Vector{Measurement{Float64}}, orbitDisplacement::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	return cross(angularVelocity, orbitDisplacement);
end

# Derive acceleration from mass. 
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsForceAcceleration(force::Vector{Measurement{Float64}}, objectMass::Measurement{Float64})::Vector{Measurement{Float64}}
	return force ./ objectMass;
end

# Derive force from potential energy at specified displacement. Displacement(source -> object).
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsPotentialEnergyForce(potentialEnergyFunction::Function, objectDisplacement::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
	displacementVectorMod::Measurement{Float64} = libYukiMathVectorQuantity(objectDisplacement);
	return -libYukiMathForwardDifference(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Calculate gravitational potential between source and object. Displacement(source -> object).
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Measurement{Float64}, sourceMass::Measurement{Float64}, objectDistance::Measurement{Float64})::Measurement{Float64}
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
# Differentiable version
# Example: True.
function libYukiPhysicsGravitationalPotentialDifferentiable(gravitationalConstant::Measurement{Float64}, sourceMass::Measurement{Float64}, objectDistance)
	return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
