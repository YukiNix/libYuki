using LinearAlgebra

"""
	libYukiMathVectorUnit(vectorA)
Return the unit vector of a given vector.
# Arguments
- `vectorA::AbstractVector{<:Real}`: The input vector.
# Returns
- The unit vector of `vectorA`.
# Notes
- The unit vector is calculated by dividing the input 
vector by its quantity (magnitude).
# Example
```julia
	vecA = [
		1 ± 0.01, 
		2 ± 0.02, 
		3 ± 0.03
	];
	libYukiMathVectorUnit(vecA)
	# Output: [
		0.2673 ± 0.0031, 
		0.5345 ± 0.0052, 
		0.8018 ± 0.0037
	]
```
"""
libYukiMathVectorUnit(
	vectorA::AbstractVector{<:Real}
)::AbstractVector{<:Real} = 
	vectorA ./ (vectorA |> libYukiMathVectorQuantity);

"""
	libYukiMathVectorShiftAngularCoordinates(
		vectorA, 
		rotateAngles...
	)
	libYukiMathVectorShiftAngularCoordinates!(
		vectorOut, 
		vectorA, 
		rotateAngles...
	)
Shift the angular coordinates of a vector by specified 
angles.
# Arguments
- `vectorA`: The input vector whose angular coordinates 
are to be shifted.
- `rotateAngles...`: A variable number of angles to shift 
the angular coordinates. For a 2D vector, provide one 
angle to shift the angle in the plane. For a 3D vector, 
provide two angles to shift the polar and azimuthal angles.
# Returns
- The non-mutated version returns a new vector with the 
shifted angular coordinates.
- The mutated version updates the provided `vectorOut` with 
the shifted angular coordinates.
# Notes
- The function calculates the magnitude and angles of the 
input vector, shifts the angles by the specified amounts, 
and then converts back to Cartesian coordinates to return 
the shifted vector.
- The angles are expected to be in radians, and the function 
handles both 2D and 3D vectors appropriately based on the 
provided `Val` type.
# Example
```julia
	vec2D = [
		3.0, 
		4.0
	];
	vecOut2D = [
		NaN, 
		NaN
	];
	rotateAngle = pi / 4;
	libYukiMathVectorShiftAngularCoordinates!(
		vecOut2D, 
		vec2D, 
		rotateAngle
	)
	libYukiMathVectorShiftAngularCoordinates(
		vec2D, 
		rotateAngle
	)
	# Output: [
		-0.7071067811865477, 
		4.949747468305833
	]
	vec3D = [
		1.0, 
		1.0, 
		1.0
	];
	vecOut3D = [
		NaN, 
		NaN, 
		NaN
	];
	rotatePolarAngle = pi / 4;
	rotateAzimuthalAngle = pi / 6;
	libYukiMathVectorShiftAngularCoordinates!(
		vecOut3D, 
		vec3D, 
		rotatePolarAngle, 
		rotateAzimuthalAngle
	)
	libYukiMathVectorShiftAngularCoordinates(
		vec3D, 
		rotatePolarAngle, 
		rotateAzimuthalAngle
	)
	# Output: [
		0.4418317469947405, 
		1.6489385281812876, 
		-0.2928932188134523
	]
```
"""
function libYukiMathVectorShiftAngularCoordinates(
	vectorA::AbstractVector, 
	rotateAngles...
)
    return libYukiMathVectorShiftAngularCoordinates(
		Val(length(vectorA)), 
		vectorA, 
		rotateAngles...
	);
end
function libYukiMathVectorShiftAngularCoordinates(
	::Val{2}, 
	vectorA::AbstractVector, 
	rotateAngle
)
    vectorMod, vectorAngle = 
		libYukiMathVectorModAndAngle(vectorA);
    return libYukiMathVectorFromPolar(
		vectorMod, 
		vectorAngle + rotateAngle
	);
end
function libYukiMathVectorShiftAngularCoordinates(
	::Val{3}, 
	vectorA::AbstractVector, 
	rotatePolarAngle, 
	rotateAzimuthalAngle
)
    vectorMod, vectorPolarAngle, vectorAzimuthalAngle =
		libYukiMathVectorModAndAngle(vectorA);
    return libYukiMathVectorFromPolar(
        vectorMod,
        vectorPolarAngle + rotatePolarAngle,
        vectorAzimuthalAngle + rotateAzimuthalAngle,
    );
end
function libYukiMathVectorShiftAngularCoordinates!(
	vectorOut::AbstractVector, 
	vectorA::AbstractVector, 
	rotateAngles...
)
    return libYukiMathVectorShiftAngularCoordinates!(
		vectorOut, 
		Val(length(vectorA)), 
		vectorA, 
		rotateAngles...
	);
end
function libYukiMathVectorShiftAngularCoordinates!(
	vectorOut::AbstractVector, 
	::Val{2}, 
	vectorA::AbstractVector, 
	rotateAngle
)
    vectorMod, vectorAngle = 
		libYukiMathVectorModAndAngle(vectorA);
    return libYukiMathVectorFromPolar!(
		vectorOut, 
		vectorMod, 
		vectorAngle + rotateAngle
	);
end
function libYukiMathVectorShiftAngularCoordinates!(
	vectorOut::AbstractVector, 
	::Val{3}, 
	vectorA::AbstractVector, 
	rotatePolarAngle, 
	rotateAzimuthalAngle
)
    vectorMod, vectorPolarAngle, vectorAzimuthalAngle =
		libYukiMathVectorModAndAngle(vectorA);
    return libYukiMathVectorFromPolar!(
        vectorOut,
        vectorMod,
        vectorPolarAngle + rotatePolarAngle,
        vectorAzimuthalAngle + rotateAzimuthalAngle,
    );
end

"""
	libYukiMathVectorModAndAngle(vectorA)
	libYukiMathVectorModAndAngle(
		::Val{2}, 
		vectorA
	)
	libYukiMathVectorModAndAngle(
		::Val{3}, 
		vectorA
	)
Calculate the magnitude and angle(s) of a vector.
# Arguments
- `vectorA`: The input vector.
- `::Val{2}`: Indicates that the input vector is 
2-dimensional.
- `::Val{3}`: Indicates that the input vector is 
3-dimensional.
# Returns
- For a 2D vector, returns a tuple containing the 
magnitude and the angle (in radians) of the vector.
- For a 3D vector, returns a tuple containing the 
magnitude, the polar angle (in radians), and the 
azimuthal angle (in radians) of the vector.
# Notes
- The magnitude is calculated using the Euclidean 
norm, while the angles are calculated using the `atan` 
function to ensure the correct quadrant is determined.
- The polar angle is measured from the positive z-axis, 
while the azimuthal angle is measured in the xy-plane 
	from the positive x-axis.
# Example
```julia
	vec2D = [
		3.0, 
		4.0
	];
	mod2D, angle2D = 
		libYukiMathVectorModAndAngle(vec2D)
	# Output: 
		mod2D = 5.0, 
		angle2D = 0.9273 rad
	vec3D = [
		1.0, 
		1.0, 
		1.0
	];
	mod3D, polarAngle3D, azimuthalAngle3D = 
		libYukiMathVectorModAndAngle(vec3D)
	# Output: 
		mod3D = 1.732, 
		polarAngle3D = 0.9553 rad, 
		azimuthalAngle3D = 0.7854 rad
```
"""
function libYukiMathVectorModAndAngle(
	vectorA::AbstractVector
)
    return libYukiMathVectorModAndAngle(
		Val(length(vectorA)), 
		vectorA
	);
end
function libYukiMathVectorModAndAngle(
	::Val{2}, 
	vectorA::AbstractVector
)
    x, y = vectorA
    return LinearAlgebra.hypot(x, y), 
		atan(y, x);
end
function libYukiMathVectorModAndAngle(
	::Val{3}, 
	vectorA::AbstractVector
)
    x, y, z = vectorA;
    modXY = LinearAlgebra.hypot(x, y);
    return LinearAlgebra.hypot(modXY, z), 
		atan(modXY, z), 
		atan(y, x);
end

"""
	libYukiMathVectorFromPolar(
		vector2DMod, 
		polarAngle
	)
	libYukiMathVectorFromPolar(
		vector3DMod, 
		polarAngle, 
		azimuthalAngle
	)
	libYukiMathVectorFromPolar!(
		vector2D, 
		vector2DMod, 
		polarAngle
	)
	libYukiMathVectorFromPolar!(
		vector3D, 
		vector3DMod, 
		polarAngle, 
		azimuthalAngle
	)
Calculate the Cartesian coordinates of a vector from 
its polar coordinates.
# Arguments
- `vector2DMod`: The magnitude of the 2D vector.
- `polarAngle`: The polar angle (in radians) for the 
2D vector.
- `vector3DMod`: The magnitude of the 3D vector.
- `azimuthalAngle`: The azimuthal angle (in radians) 
for the 3D vector.
- `vector2D`: A vector to store the Cartesian 
coordinates of the 2D vector.
- `vector3D`: A vector to store the Cartesian 
coordinates of the 3D vector.
# Returns
- The Cartesian coordinates of the vector in the form 
of an `SVector` for the non-mutating versions, or updates 
	the provided vectors for the mutating versions.
# Notes
- The non-mutating versions return a new `SVector` 
containing the Cartesian coordinates, while the mutating 
versions update the provided vectors in place.
- The polar angle is measured from the positive x-axis 
in the 2D case, and from the positive z-axis in the 3D 
case, while the azimuthal angle is measured in the 
xy-plane from the positive x-axis.	
# Example
```julia
	vec2D = [
		NaN, 
		NaN
	];
	vec3D = [
		NaN, 
		NaN, 
		NaN
	];
	mod = 1.0 ± 0.1;
	polarAngle = pi / 4 ± 0.01;
	azimuthalAngle = pi / 6 ± 0.02;
	libYukiMathVectorFromPolar(
		mod, 
		polarAngle)
	libYukiMathVectorFromPolar!(
		vec2D, 
		mod, 
		polarAngle
	)
	# Output: [
		0.707 ± 0.071, 
		0.707 ± 0.071
	]
	libYukiMathVectorFromPolar(
		mod, 
		polarAngle, 
		azimuthalAngle
	)
	libYukiMathVectorFromPolar!(
		vec3D, 
		mod, 
		polarAngle, 
		azimuthalAngle
	)
	# Output: [
		0.612 ± 0.062, 
		0.354 ± 0.038, 
		0.707 ± 0.071
	]
```
"""
function libYukiMathVectorFromPolar(
	vector2DMod, 
	polarAngle
)
	sinAngle, cosAngle = sincos(polarAngle);
	return SVector(
		vector2DMod * cosAngle, 
		vector2DMod * sinAngle
	);
end
function libYukiMathVectorFromPolar(
	vector3DMod, 
	polarAngle, 
	azimuthalAngle
)
	sinPolarAngle, cosPolarAngle = sincos(polarAngle);
	sinAzimuthalAngle, cosAzimuthalAngle = 
		sincos(azimuthalAngle);
	return SVector(
		vector3DMod * 
			sinPolarAngle * cosAzimuthalAngle,
		vector3DMod * 
			sinPolarAngle * sinAzimuthalAngle,
		vector3DMod * cosPolarAngle
	);
end
function libYukiMathVectorFromPolar!(
	vector2D, 
	vector2DMod, 
	polarAngle
)
    sinAngle, cosAngle = sincos(polarAngle);
    vector2D[1] = vector2DMod * cosAngle;
    vector2D[2] = vector2DMod * sinAngle;
end
function libYukiMathVectorFromPolar!(
	vector3D, 
	vector3DMod, 
	polarAngle, 
	azimuthalAngle
)
	sinPolarAngle, cosPolarAngle = sincos(polarAngle);
	sinAzimuthalAngle, cosAzimuthalAngle = 
		sincos(azimuthalAngle);
	vector3D[1] = vector3DMod * 
		sinPolarAngle * cosAzimuthalAngle;
	vector3D[2] = vector3DMod * 
		sinPolarAngle * sinAzimuthalAngle;
	vector3D[3] = vector3DMod * cosPolarAngle;
end

"""
	libYukiVectorGradient(
		f, 
		vector)
Calculate the gradient of a scalar function `f` with respect to a vector variable.
# Arguments
- `f`: The scalar function for which to calculate the gradient. It should take a vector as input and return a scalar.
- `vector`: The vector variable with respect to which the gradient is calculated. It should be an `AbstractVector` of real numbers.
# Returns
- The gradient of the function `f` with respect to the input `vector`, returned as an `SVector` of the same length as `vector`.
# Notes
- The function uses forward-mode automatic differentiation provided by the `ForwardDiff` package to compute the gradient efficiently and accurately.
- The input function `f` should be compatible with the automatic differentiation, meaning it should be differentiable and should not contain operations that are not supported by `ForwardDiff`.
# Example
```julia
	f(v) = v[1]^2 + v[2]^2 + v[3]^2;
	vec = [1.0, 2.0, 3.0];
	libYukiMathVectorGradient(f, vec)
# Output: SVector{3}([2.0, 4.0, 6.0])
```
"""
libYukiMathVectorGradient(
	f, 
	vector
) = 
	ForwardDiff.gradient(f, vector);


"""
	libYukiMathVectorDotProduct(
		vectorA, 
		vectorB
	)
Return the dot product of two vectors.
# Arguments
- `vectorA`: The first vector.
- `vectorB`: The second vector.
# Returns
- The dot product of `vectorA` and `vectorB`.
# Notes
- The dot product is calculated as the sum of the products 
of corresponding elements of the two vectors.
# Example
```julia
	vecA = [
		1 ± 0.01, 
		2 ± 0.02, 
		3 ± 0.03
	];
	vecB = [
		4 ± 0.04, 
		5 ± 0.05, 
		6 ± 0.06
	];
	libYukiMathVectorDotProduct(
		vecA, 
		vecB
	)
	# Output: 32.0 ± 0.3
```
"""
libYukiMathVectorDotProduct(
	vectorA, 
	vectorB
) = 
	LinearAlgebra.dot(vectorA, vectorB);

"""
	libYukiMathVectorCrossProduct(
		vectorA, 
		vectorB
	)
Return the cross product of two vectors.
# Arguments
- `vectorA`: The first vector.
- `vectorB`: The second vector.
# Returns
- The cross product of `vectorA` and `vectorB`.
# Notes
- The cross product is calculated using the right-hand rule.
# Example
```julia
	vecA = [
		1 ± 0.01, 
		2 ± 0.02, 
		3 ± 0.03
	];
	vecB = [
		4 ± 0.04, 
		5 ± 0.05, 
		6 ± 0.06
	];
	libYukiMathVectorCrossProduct(
		vecA, 
		vecB
	)
	# Output: [
		-3.0 ± 0.27, 
		6.0 ± 0.19, 
		-3.0 ± 0.13
	]
```
"""
libYukiMathVectorCrossProduct(
	vectorA, 
	vectorB
) = 
	LinearAlgebra.cross(vectorA, vectorB);

"""
	libYukiMathVectorQuantity(vectorA)
Return the quantity (magnitude) of a vector.
# Arguments
- `vectorA`: The vector.
# Returns
- The quantity (magnitude) of `vectorA`.
# Notes
- The quantity is calculated using the Euclidean norm.
# Example
```julia
	vecA = [
		1 ± 0.01, 
		1 ± 0.1, 
		1 ± 0.2
	];
	libYukiMathVectorQuantity(vecA) 
	# Output: 1.732 ± 0.2
```
"""
libYukiMathVectorQuantity(
	vectorA
) = 
	LinearAlgebra.norm(vectorA);

"""
	libYukiMathVectorFieldDivergence(
		f,
		x
	)
Calculate the divergence of a vector field `f` with respect to 
the vector of variables `x`.
# Arguments
- `f`: The vector field function for which to calculate the divergence.
- `x`: The point at which to evaluate the divergence.
# Returns
- The divergence of `f` with respect to `x` evaluated at the 
point `x`.
# Notes
- This function is efficient and accurate for computing 
the divergence of multivariate vector fields.
# Example
```julia
	f(x) = [x[1], x[2], x[3]];
	libYukiMathVectorFieldDivergence(
		f, 
		[1.0, 2.0, 3.0]
	)
	# Output: 3.0
``` 
"""
libYukiMathVectorFieldDivergence(
	f, 
	x::AbstractVector{<:Real}
) =
    LinearAlgebra.tr(libYukiMathJacobian(f, x));

"""
	libYukiMathVectorFieldCurl(
		f,
		x
	)
Calculates the curl of a 3-dimensional vector field `f` with respect to 
the vector of variables `x`.
# Arguments
- `f`: The 3-dimensional vector field function for which to calculate the curl.
- `x`: The point at which to evaluate the curl.
# Returns
- The curl of `f` with respect to `x` evaluated at the point `x`.
# Notes
- This function is defined only for 3-dimensional vector fields.
# Example
```julia
	f(x) = [
		x[1], 
		x[2], 
		x[3]
	];
	libYukiMathVectorFieldCurl(
		f, 
		[
			1.0, 
			2.0, 
			3.0
		]
	)
	# Output: [
		0.0, 
		0.0, 
		0.0
	]
``` 
"""
function libYukiMathVectorFieldCurl(
    f,
    x::AbstractVector{<:Real},
)
    length(x) == 3 ||
        throw(ArgumentError(
			"curl is defined here only for 3-dimensional" * 
			" vector fields"))
    jacobian = libYukiMathJacobian(f, x);
    return SVector(
        jacobian[3, 2] - jacobian[2, 3],
        jacobian[1, 3] - jacobian[3, 1],
        jacobian[2, 1] - jacobian[1, 2],
    );
end
