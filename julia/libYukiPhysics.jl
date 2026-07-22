using Measurements, ForwardDiff, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysicsBody.jl")
include("libYukiPhysicsField.jl")

"""
	libYukiPhysicsKineticEnergy(
		body
	)
	libYukiPhysicsKineticEnergy(
		mass, 
		velocity
	)
Calculate the kinetic energy of a body given its mass 
and velocity.
# Arguments
- `body`: A `libYukiPhysicsBody` instance.
- `mass`: Mass of the body.
- `velocity`: Velocity vector of the body (3D).
# Returns
- Kinetic energy of the body.
# Example
```julia
	libYukiPhysicsKineticEnergy(
		100 ± 0.3,
		[
			0.1 ± 0, 
			1 ± 0.1, 
			0.2 ± 0.01
		]
	)
	# Output: 52.0 ± 10.0
```
"""
function libYukiPhysicsKineticEnergy(
	body::libYukiPhysicsBody
)
	return libYukiPhysicsKineticEnergy.(
		body.mass, 
		body.trajectory.states[end].velocity
	);
end
function libYukiPhysicsKineticEnergy(
	mass::Real, 
	velocity::AbstractVector{<:Real}
)
	length(velocity) == 3 || 
		throw(ArgumentError("velocity must have length 3"))
	return libYukiPhysicsKineticEnergy(
		mass, 
		SVector{3}(velocity)
	);
end
libYukiPhysicsKineticEnergy(
	mass::Real, 
	velocity::SVector{3, <:Real}
) = 
	0.5 * mass * libYukiMathVectorDotProduct(
		velocity, 
		velocity
	);

"""
	libYukiPhysicsMomentum(
		body
	)
	libYukiPhysicsMomentum(
		mass, 
		velocity
	)
Calculate the momentum of a body given its mass and velocity.
# Arguments
- `body`: A `libYukiPhysicsBody` instance.
- `mass`: Mass of the body.
- `velocity`: Velocity vector of the body (3D).
# Returns
- Momentum vector of the body.
# Example
```julia
	libYukiPhysicsMomentum(
		100 ± 0.3,
		[
			0.1 ± 0, 
			1 ± 0.1, 
			0.2 ± 0.01
		]
	)
	# Output [
		10.0 ± 0.03
		100.0 ± 10.0
  		20.0 ± 1.0
	]
```
"""
libYukiPhysicsMomentum(
	body::libYukiPhysicsBody
) = 
	libYukiPhysicsMomentum.(
		body.mass, 
		body.trajectory.states[end].velocity
	);
function libYukiPhysicsMomentum(
	mass::Real, 
	velocity::AbstractVector{<:Real}
) 
	length(velocity) == 3 || 
		throw(ArgumentError("velocity must have length 3"))
	return libYukiPhysicsMomentum(
		mass, 
		SVector{3}(velocity)
	);
end
libYukiPhysicsMomentum(
	mass::Real, 
	velocity::SVector{3, <:Real}
) = mass .* velocity;






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
	return -libYukiMathDerivative(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Differentiable version (strips Measurement wrapper for ForwardDiff compatibility).
# When objectDistance is also Measurement, use the standard version above via full uncertainty propagation.
# When objectDistance is a non-Measurement Real (e.g. ForwardDiff.Dual), strip Measurement from G and M.
# Dependency: Measurements.
# Example: True.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Real, sourceMass::Real, objectDistance::Real) 
	return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
# function libYukiPhysicsGravitationalPotential(gravitationalConstant::Real, sourceMass::Real, objectDistance::Real)
# 	return -(gravitationalConstant * sourceMass) / objectDistance;
# end
