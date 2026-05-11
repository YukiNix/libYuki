using Measurements, ForwardDiff, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

# Physical object data structure
# Example: True.
mutable struct libYukiPhysicsBody{T <: Real}
	name::String
	position::AbstractVector{AbstractVector{<:Real}}
	velocity::AbstractVector{AbstractVector{<:Real}}
	mass::Real
	charge::Real
	radius::Real
	libYukiPhysicsBody(name, position::AbstractVector{<:Real}, velocity::AbstractVector{<:Real}, mass::Real, charge::Real, radius::Real) = new{<:Real}(name, [position], [velocity], mass, charge, radius)
end

# Kinetic energy of a body.
# TODO: Validate & Example.
function libYukiPhysicsKineticEnergy(bodyMass::Real, bodyVelocity::AbstractVector{<:Real})::Real 
	return 0.5 * bodyMass * libYukiMathVectorQuantity(bodyVelocity .* bodyVelocity);
end

# Derive gravitational circular motion's orbit radius from source mass and orbit period.
# TODO: Validate & Example.
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassOrbitPeriod(sourceMass::Real, orbitPeriod::Real, gravitationalConstant::Real)::Real
	return (gravitationalConstant * sourceMass * (orbitPeriod / (2 * π)) ^ 2) ^ (1 / 3.);
end

# Derive gravitational circular motion's orbit radius from source mass and angular velocity quantity.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionOrbitRadiusFromSourceMassAngularVelocity(sourceMass::Real, angularVelocityQuantity::Real, gravitationalConstant::Real)::Real
	return (gravitationalConstant * sourceMass / (angularVelocityQuantity ^ 2)) ^ (1 / 3.);
end

# Derive gravitational circular motion's velocity quantity from source mass and orbit radius.
# TODO: Validate & Example. 
function libYukiPhysicsGravitationalCircularMotionVelocityQuantity(sourceMass::Real, orbitRadius::Real, gravitationalConstant::Real)::Real 
	return sqrt(gravitationalConstant * sourceMass / orbitRadius);
end

# Derive circular motion's velocity from acceleration and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAccelerationAngularVelocity(acceleration::AbstractVector{<:Real}, angularVelocity::AbstractVector{<:Real})::AbstractVector{<:Real} 
	return libYukiMathUnitVector(cross(acceleration, angularVelocity)) .* (libYukiMathVectorQuantity(acceleration) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's orbit radius from velocity and angular velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionOrbitDisplacementFromVelocityAngularVelocity(velocity::AbstractVector{<:Real}, angularVelocity::AbstractVector{<:Real})::AbstractVector{<:Real} 
	return libYukiMathUnitVector(cross(velocity, angularVelocity)) .* (libYukiMathVectorQuantity(velocity) / libYukiMathVectorQuantity(angularVelocity));
end

# Derive circular motion's angular velocity to velocity.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionVelocityFromAngularVelocity(angularVelocity::AbstractVector{<:Real}, orbitDisplacement::AbstractVector{<:Real})::AbstractVector{<:Real} 
	return cross(angularVelocity, orbitDisplacement);
end

# Derive acceleration from mass. 
# Example: True.
function libYukiPhysicsForceAcceleration(force::AbstractVector{<:Real}, objectMass::Real)::AbstractVector{<:Real} 
	return force ./ objectMass;
end

# Derive force from potential energy at specified displacement. Displacement(source -> object).
# Example: True.
function libYukiPhysicsPotentialEnergyForce(potentialEnergyFunction::Function, objectDisplacement::AbstractVector{<:Real})::AbstractVector{<:Real} 
	displacementVectorMod = libYukiMathVectorQuantity(objectDisplacement);
	return -libYukiMathForwardDifference(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Calculate gravitational potential between source and object. Displacement(source -> object).
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Real, sourceMass::Real, objectDistance::Real)::Real 
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
# Differentiable version (strips Measurement wrapper for ForwardDiff compatibility).
# When objectDistance is also Measurement, use the standard version above via full uncertainty propagation.
# When objectDistance is a non-Measurement Real (e.g. ForwardDiff.Dual), strip Measurement from G and M.
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Real, sourceMass::Real, objectDistance::Real) 
	return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Real, sourceMass::Real, objectDistance::Real)
	return -(gravitationalConstant * sourceMass) / objectDistance;
end
