abstract type libYukiPhysicsFieldAbstractField end
abstract type libYukiPhusicsFieldAbstractPotentialField <: libYukiPhysicsFieldAbstractField end 
abstract type libYukiPhusicsFieldAbstractVectorField <: libYukiPhysicsFieldAbstractField end 

"""
    libYukiPhysicsFieldPotentialField(
        potentialFunction::F
    )
Calculate a potential field given a potential function.
# Arguments
- `potentialFunction`: A function that takes a position and returns 
the potential at that position.
# Returns
- An instance of `libYukiPhysicsFieldPotentialField`.
# Example
```julia
    potentialField = libYukiPhysicsFieldPotentialField(
        position -> position[1] ^ 2 + position[2] ^ 2 + position[3] ^ 2
    )
    potentialField([
        1, 
        2, 
        3
    ])
    # Output: 14
```
"""
struct libYukiPhysicsFieldPotentialField{F}<:libYukiPhusicsFieldAbstractPotentialField
    potentialFunction::F
end
(fieldObject::libYukiPhysicsFieldPotentialField{F})(position) where F =
    fieldObject.potentialFunction(position);

"""
    libYukiPhysicsFieldVectorField(
        vectorFunction::F
    )
Calculate a vector(strengths) field given a vector function.
# Arguments
- `vectorFunction`: A function that takes a position and returns 
the vector of strengths at that position.
# Returns
- An instance of `libYukiPhysicsFieldVectorField`.
# Example
```julia
    vectorField = libYukiPhysicsFieldVectorField(
        position -> [position[1], position[2], position[3]]
    )
    vectorField([
        1, 
        2, 
        3
    ])
    # Output: [
        1, 
        2, 
        3
    ]
```
"""
struct libYukiPhysicsFieldVectorField{F}<:libYukiPhusicsFieldAbstractVectorField
    vectorFunction::F
end
(fieldObject::libYukiPhysicsFieldVectorField{F})(position) where F =
    fieldObject.vectorFunction(position);

"""
    libYukiPhysicsFieldStrength(
        potentialField, 
        x
    )
Calculate the strength of the field at a given position.
# Arguments
- `potentialField`: An callable instance (object) of `libYukiPhysicsFieldPotentialField`.
- `x`: The position at which to calculate the field strength.
# Returns
- The strength of the field at the given position.
# Example
```julia
    potentialField = libYukiPhysicsFieldPotentialField(
        position -> position[1] ^ 2 + position[2] ^ 2 + position[3] ^ 2
    )
    libYukiPhysicsFieldStrength(potentialField, [
        1, 
        2, 
        3
    ])
    # Output: [
        -2, 
        -4, 
        -6
    ]
``` 
"""
function libYukiPhysicsFieldStrength(potentialField::libYukiPhysicsFieldPotentialField, x)
    return -libYukiMathDerivative(potentialField, x);
end

# function libYukiPhysicsFieldGravitationalField()
#     libYukiPhysicsFieldVectorField(
#         position -> begin
#             G = 6.67430e-11
#             mass = 1.0
#             r = sqrt(position[1]^2 + position[2]^2 + position[3]^2)
#             return SVector(
#                 -G * mass * position[1] / r^3,
#                 -G * mass * position[2] / r^3,
#                 -G * mass * position[3] / r^3
#             )
#         end
#     )
# end