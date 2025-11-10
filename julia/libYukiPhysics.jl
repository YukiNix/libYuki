using Measurements, ForwardDiff, DifferentialEquations, LinearAlgebra

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")

mutable struct libYukiPhysicsBody
    name::String
    position::Vector{Measurement{Float64}}
    velocity::Vector{Measurement{Float64}}
    mass::Measurement{Float64}
    charge::Measurement{Float64}
    radius::Measurement{Float64}
    libYukiPhysicsBody(name, position, velocity, mass, charge, radius) = new(name, position, velocity, mass, charge, radius)
end

function libYukiPhysicsGravitationalNBodySimulation(timeStart::Measurement{Float64}, timeEnd::Measurement{Float64}, timeStep::Measurement{Float64}, bodies::Vector{libYukiPhysicsBody}, integrator, gravitationalConstant::Measurement{Float64})

    function libYukiPhysicsGravitationalNBodySimulationAcceleration!(∂velocity::Vector{Measurement{Float64}}, velocity::Vector{Measurement{Float64}}, position::Vector{Measurement{Float64}}, partical, time::Measurement{Float64})
        bodiesNumber::Int64 = length(partical.masses);
        dimension::Int64 = 3;
        separatedPositions = reshape(position, dimension, bodiesNumber);
        ∂velocity .= 0.0 ± 0.0;
        for bodyAIndex::Int64 in 1 : bodiesNumber
            bodyAAcceleration::Vector{Measurement{Float64}} = zeros(Measurement{Float64}, dimension);
            for bodyBIndex::Int64 in 1 : bodiesNumber
                if bodyAIndex != bodyBIndex
                    displacementBA::Vector{Measurement{Float64}} = separatedPositions[:, bodyAIndex] - separatedPositions[:, bodyBIndex];
                    bodyAAcceleration += libYukiPhysicsForceAcceleration(
                        libYukiPhysicsPotentialEnergyForce(
                            x -> Measurements.value(partical.masses[bodyAIndex]) * libYukiPhysicsGravitationalPotentialDifferentiable(
                                gravitationalConstant, 
                                partical.masses[bodyBIndex], 
                                x),
                            displacementBA), 
                        partical.masses[bodyAIndex]);
                end
            end
            ∂velocity[(bodyAIndex - 1) * dimension + 1 : bodyAIndex * dimension] = bodyAAcceleration;
        end
    end

    masses::Vector{Measurement{Float64}} = map(x -> x.mass, bodies);
    velocities::Vector{Measurement{Float64}} = [body.velocity[dimension] for body in bodies for dimension in 1 : 3];
    positions::Vector{Measurement{Float64}} = [body.position[dimension] for body in bodies for dimension in 1 : 3];

    NBodyProblem = SecondOrderODEProblem(
        libYukiPhysicsGravitationalNBodySimulationAcceleration!,
        velocities,
        positions,
        (timeStart, timeEnd),
        (masses = masses, )
    );
    NBodySolution = solve(NBodyProblem, integrator, dt = timeStep);

    return NBodySolution;
end

# Derive acceleration from mass. 
# Dependency: Measurements.
function libYukiPhysicsForceAcceleration(force::Vector{Measurement{Float64}}, objectMass::Measurement{Float64})::Vector{Measurement{Float64}}
    return force ./ objectMass;
end

# Derive force from potential energy at specified displacement. Displacement(source -> object).
# Dependency: Measurements.
function libYukiPhysicsPotentialEnergyForce(potentialEnergyFunction::Function, objectDisplacement::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
    displacementVectorMod::Measurement{Float64} = libYukiMathVectorToDistance(objectDisplacement);
    return -libYukiMathForwardDifference(x -> potentialEnergyFunction(x), displacementVectorMod) .* objectDisplacement ./ displacementVectorMod;
end

# Calculate gravitational potential between source and object. Displacement(source -> object).
# Dependency: Measurements.
function libYukiPhysicsGravitationalPotential(gravitationalConstant::Measurement{Float64}, sourceMass::Measurement{Float64}, objectDistance::Measurement{Float64})::Measurement{Float64}
    return -(gravitationalConstant * sourceMass) / objectDistance;
end
# Differentiable version
function libYukiPhysicsGravitationalPotentialDifferentiable(gravitationalConstant::Measurement{Float64}, sourceMass::Measurement{Float64}, objectDistance)
    return -Measurements.value(gravitationalConstant * sourceMass) / objectDistance;
end
