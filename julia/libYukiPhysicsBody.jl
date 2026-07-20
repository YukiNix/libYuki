using StaticArrays

"""
	libYukiPhysicsBodyState(
		position, 
		velocity
	)
Define the state of a physics body with position and velocity.
# Arguments
- `position`: Position vector of the body (3D).
- `velocity`: Velocity vector of the body (3D).
# Returns
- `libYukiPhysicsBodyState` instance.
# Example
```julia
	state = libYukiPhysicsBodyState(
		[0.1, 0.2, 0.3], 
		[1.0, 0.0, 0.0]
	)
``` 
"""
struct libYukiPhysicsBodyState{
		TP<:Real,
		TV<:Real
}
	position::SVector{3, TP}
	velocity::SVector{3, TV}
end
function libYukiPhysicsBodyState(
	position::AbstractVector{TP}, 
	velocity::AbstractVector{TV}
) where {TP<:Real, TV<:Real}
	length(position) == 3 || 
		throw(ArgumentError("position must have length 3"))
	length(velocity) == 3 || 
		throw(ArgumentError("velocity must have length 3"))
	return libYukiPhysicsBodyState{TP, TV}(
		SVector{3}(position), 
		SVector{3}(velocity)
	);
end

"""
	libYukiPhysicsBodyTrajectory(
		states
	)
Define the trajectory of a physics body as a sequence of states.
# Arguments
- `states`: A vector of `libYukiPhysicsBodyState` instances representing the states of the body over time.
# Returns
- `libYukiPhysicsBodyTrajectory` instance.
# Example
```julia
	states = [
		libYukiPhysicsBodyState([0.1, 0.2, 0.3], [1.0, 0.0, 0.0]),
		libYukiPhysicsBodyState([0.2, 0.3, 0.4], [0.0, 1.0, 0.0])
	]
	trajectory = libYukiPhysicsBodyTrajectory(states)
``` 
"""
struct libYukiPhysicsBodyTrajectory
	states::Vector{libYukiPhysicsBodyState}
end
function libYukiPhysicsBodyTrajectory(
	states::Vector{libYukiPhysicsBodyState}
)
	return libYukiPhysicsBodyTrajectory(
		states
	);
end

"""
	libYukiPhysicsBody(
		name, 
		mass, 
		charge, 
		radius, 
		trajectory
	)
Define a physics body with mass, charge, 
radius, and trajectory.
# Arguments
- `name`: Name of the body.
- `mass`: Mass of the body.
- `charge`: Electric charge of the body.
- `radius`: Radius of the body.
- `trajectory`: Trajectory of the body as a `libYukiPhysicsBodyTrajectory` instance.
# Returns
- `libYukiPhysicsBody` instance.
# Example
```julia
	body = libYukiPhysicsBody(
		"bodyA", 
		100 ± 0.3, 
		12 ± 0.24, 
		50 ± 1,
		libYukiPhysicsBodyTrajectory([])
	)
```
"""
struct libYukiPhysicsBody{
		TM<:Real, 
		TC<:Real, 
		TR<:Real
}
	name::String
	mass::TM
	charge::TC
	radius::TR
end
function libYukiPhysicsBody(
	name::String, 
	mass::TM, 
	charge::TC, 
	radius::TR,
	trajectory::libYukiPhysicsBodyTrajectory
) where {TM<:Real, TC<:Real, TR<:Real}
	return libYukiPhysicsBody{TM, TC, TR}(
		name, 
		mass, 
		charge, 
		radius,
		trajectory
	);
end
