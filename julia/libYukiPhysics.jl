using Measurements, ForwardDiff, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

# Physical object data structure
# Example: True.
mutable struct libYukiPhysicsBody{T <: AbstractFloat}
	name::String
	position::Vector{Vector{T}}
	velocity::Vector{Vector{T}}
	mass::T
	charge::T
	radius::T
	libYukiPhysicsBody(name, position::Vector{T}, velocity::Vector{T}, mass::T, charge::T, radius::T) where {T <: AbstractFloat} = new{T}(name, [position], [velocity], mass, charge, radius)
end

# Kinetic energy of a body.
# TODO: Validate & Example.
function libYukiPhysicsKineticEnergy(bodyMass::T, bodyVelocity::Vector{T})::T where {T <: AbstractFloat}
	return 0.5 * bodyMass * libYukiMathVectorQuantity(bodyVelocity .* bodyVelocity);
end

# Derive gravitational circular motion's orbit radius from source mass and angular velocity quantity.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassAngularVelocity(sourceMass::T, angularVelocityQuantity::T, gravitationalConstant::T)::T where {T <: AbstractFloat}
	return (gravitationalConstant * sourceMass / (angularVelocityQuantity ^ 2)) ^ (1 / 3.);
end

# Derive gravitational circular motion's velocity quantity from source mass and orbit radius.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionVelocityQuantity(sourceMass::T, orbitRadius::T, gravitationalConstant::T)::T where {T <: AbstractFloat}
	return sqrt(gravitationalConstant * sourceMass / orbitRadius);
end

# Derive circular motion's velocity from acceleration and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAccelerationAngularVelocity(acceleration::Vector{T}, angularVelocity::Vector{T})::Vector{T} where {T <: AbstractFloat}
	return libYukiMathUnitVector(cross(acceleration, angularVelocity)) .* (libYukiMathVectorQuantity(acceleration) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's orbit radius from velocity and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionOrbitDisplacementFromVelocityAngularVelocity(velocity::Vector{T}, angularVelocity::Vector{T})::Vector{T} where {T <: AbstractFloat}
	return libYukiMathUnitVector(cross(velocity, angularVelocity)) .* (libYukiMathVectorQuantity(velocity) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's angular velocity to velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAngularVelocity(angularVelocity::Vector{T}, orbitDisplacement::Vector{T})::Vector{T} where {T <: AbstractFloat}
	return cross(angularVelocity, orbitDisplacement);
end

# Derive acceleration from mass. 
# Example: True.
function libYukiPhysicsForceAcceleration(force::Vector{T}, objectMass::T)::Vector{T} where {T <: AbstractFloat}
	return force ./ objectMass;
end

# Derive force from potential energy at specified displacement. Displacement(source -> object).
# Example: True.
function libYukiPhysicsPotentialEnergyForce(potentialEnergyFunction::Function, objectDisplacement::Vector{T})::Vector{T} where {T <: AbstractFloat}
	displacementVectorMod = libYukiMathVectorQuantity(objectDisplacement);
	return -libYukiMathForwardDifference(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Calculate gravitational potential between source and object. Displacement(source -> object).
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::T, sourceMass::T, objectDistance::T)::T where {T <: AbstractFloat}
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
# Differentiable version (strips Measurement wrapper for ForwardDiff compatibility).
# When objectDistance is also Measurement, use the standard version above via full uncertainty propagation.
# When objectDistance is a non-Measurement Real (e.g. ForwardDiff.Dual), strip Measurement from G and M.
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Measurement{T}, sourceMass::Measurement{T}, objectDistance::S) where {T <: AbstractFloat, S <: Real}
	return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
function libYukiPhysicsGravitationalPotential(gravitationalConstant::T, sourceMass::T, objectDistance::T) where {T <: Measurement}
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
