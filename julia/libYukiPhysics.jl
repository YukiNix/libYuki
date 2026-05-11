using Measurements, ForwardDiff, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

# Physical object data structure
# Example: True.
mutable struct libYukiPhysicsBody{T <: Real}
	name::String
	position::AbstractVector{AbstractVector{T}}
	velocity::AbstractVector{AbstractVector{T}}
	mass::T
	charge::T
	radius::T
	libYukiPhysicsBody(name, position::AbstractVector{T}, velocity::AbstractVector{T}, mass::T, charge::T, radius::T) where {T <: Real} = new{T}(name, [position], [velocity], mass, charge, radius)
end

# Kinetic energy of a body.
# TODO: Validate & Example.
function libYukiPhysicsKineticEnergy(bodyMass::T, bodyVelocity::AbstractVector{T})::T where {T <: Real}
	return 0.5 * bodyMass * libYukiMathVectorQuantity(bodyVelocity .* bodyVelocity);
end

# Derive gravitational circular motion's orbit radius from source mass and orbit period.
# TODO: Validate & Example.
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassOrbitPeriod(sourceMass::T, orbitPeriod::T, gravitationalConstant::T)::T where {T <: Real}
	return (gravitationalConstant * sourceMass * (orbitPeriod / (2 * π)) ^ 2) ^ (1 / 3.);
end

# Derive gravitational circular motion's orbit radius from source mass and angular velocity quantity.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassAngularVelocity(sourceMass::T, angularVelocityQuantity::T, gravitationalConstant::T)::T where {T <: Real}
	return (gravitationalConstant * sourceMass / (angularVelocityQuantity ^ 2)) ^ (1 / 3.);
end

# Derive gravitational circular motion's velocity quantity from source mass and orbit radius.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionVelocityQuantity(sourceMass::T, orbitRadius::T, gravitationalConstant::T)::T where {T <: Real}
	return sqrt(gravitationalConstant * sourceMass / orbitRadius);
end

# Derive circular motion's velocity from acceleration and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAccelerationAngularVelocity(acceleration::AbstractVector{T}, angularVelocity::AbstractVector{T})::AbstractVector{T} where {T <: Real}
	return libYukiMathUnitVector(cross(acceleration, angularVelocity)) .* (libYukiMathVectorQuantity(acceleration) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's orbit radius from velocity and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionOrbitDisplacementFromVelocityAngularVelocity(velocity::AbstractVector{T}, angularVelocity::AbstractVector{T})::AbstractVector{T} where {T <: Real}
	return libYukiMathUnitVector(cross(velocity, angularVelocity)) .* (libYukiMathVectorQuantity(velocity) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's angular velocity to velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAngularVelocity(angularVelocity::AbstractVector{T}, orbitDisplacement::AbstractVector{T})::AbstractVector{T} where {T <: Real}
	return cross(angularVelocity, orbitDisplacement);
end

# Derive acceleration from mass. 
# Example: True.
function libYukiPhysicsForceAcceleration(force::AbstractVector{T}, objectMass::T)::AbstractVector{T} where {T <: Real}
	return force ./ objectMass;
end

# Derive force from potential energy at specified displacement. Displacement(source -> object).
# Example: True.
function libYukiPhysicsPotentialEnergyForce(potentialEnergyFunction::Function, objectDisplacement::AbstractVector{T})::AbstractVector{T} where {T <: Real}
	displacementVectorMod = libYukiMathVectorQuantity(objectDisplacement);
	return -libYukiMathForwardDifference(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Calculate gravitational potential between source and object. Displacement(source -> object).
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::T, sourceMass::T, objectDistance::T)::T where {T <: Real}
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
# Differentiable version (strips Measurement wrapper for ForwardDiff compatibility).
# When objectDistance is also Measurement, use the standard version above via full uncertainty propagation.
# When objectDistance is a non-Measurement Real (e.g. ForwardDiff.Dual), strip Measurement from G and M.
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Measurement{T}, sourceMass::Measurement{T}, objectDistance::S) where {T <: Real, S <: Real}
	return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
function libYukiPhysicsGravitationalPotential(gravitationalConstant::T, sourceMass::T, objectDistance::T) where {T <: Measurement}
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
