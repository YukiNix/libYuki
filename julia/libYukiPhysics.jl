using Measurements, ForwardDiff, DifferentialEquations, LinearAlgebra, JLD2

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

# Derive circular motion's angular velocity to velocity.
# Dependency: Measurements.
# TODO: Validate & Example. 
function libYukiPhysicsCircularMotionAngularVelocityToVelocity(angularVelocity::Measurement{Float64}, orbitDisplacement::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
    displacementQuantity::Measurement{Float64}, displacementAngle::Measurement{Float64} = libYukiMath2DVectorModAndAngle(orbitDisplacement);
    return [-sin(displacementAngle), cos(displacementAngle)] .* libYukiPhysicsCircularMotionAngularVelocityToVelocityQuantity(angularVelocity, displacementQuantity);
end

# Derive circular motion's velocity quantity to angular velocity.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiPhysicsCircularMotionVelocityQuantityToAngularVelocity(velocityQuantity::Measurement{Float64}, orbitRadius::Measurement{Float64})::Measurement{Float64}
    return velocityQuantity / orbitRadius;
end

# Derive circular motion's angular velocity to velocity quantity.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiPhysicsCircularMotionAngularVelocityToVelocityQuantity(angularVelocity::Measurement{Float64}, orbitRadius::Measurement{Float64})::Measurement{Float64}
    return angularVelocity * orbitRadius;
end

# Load streaming saved gravitational N-body simulation result. 
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiPhysicsGravitationalNBodySimulationStreamingLoad(savingDirectory::String)
    splitIndex::Int64 = 0;
    timesLoaded = nothing;
    bodiesLoaded = nothing;
    while ispath("$savingDirectory/libYukiPhysicsGravitationalNBodySimulationStreamingResult_$splitIndex.jld2")
        @load "$savingDirectory/libYukiPhysicsGravitationalNBodySimulationStreamingResult_$splitIndex.jld2" times bodies
        if !isnothing(timesLoaded)
            append!(timesLoaded, times);
            for bodiesIndex::Int64 in 1 : length(bodies)
                append!(bodiesLoaded[bodiesIndex].velocity, bodies[bodiesIndex].velocity);
                append!(bodiesLoaded[bodiesIndex].position, bodies[bodiesIndex].position);
            end
        else
            timesLoaded = times;
            bodiesLoaded = bodies;
        end
        splitIndex = splitIndex + 1;
    end
    return timesLoaded, bodiesLoaded;
end

# Gravitational N-body simulation (stream save to file). 
# Dependency: Measurements, libYukiMath.
# TODO: Validate & Example.
function libYukiPhysicsGravitationalNBodySimulationStreaming(timeStart::Measurement{Float64}, timeEnd::Measurement{Float64}, timeStep::Measurement{Float64}, bodiesSimulation::Vector{libYukiPhysicsBody}, integrator, gravitationalConstant::Measurement{Float64}, savingDirectory::String, splitSteps::Int64)
    splitIndex::Int64 = 0;
    timesSimulation::Vector{Measurement{Float64}} = [];

    timeSimulating::Measurement{Float64} = timeStart;
    timeNext::Measurement{Float64} = timeSimulating + splitSteps * timeStep;
    mkdir(savingDirectory);
    while timeSimulating < timeEnd
        timeNext = (timeNext < timeEnd) ? timeNext : timeEnd;
        bodiesSimulation = [ 
            libYukiPhysicsBody(bodiesSimulation[bodyIndex].name, bodiesSimulation[bodyIndex].position[end], bodiesSimulation[bodyIndex].velocity[end], bodiesSimulation[bodyIndex].mass, bodiesSimulation[bodyIndex].charge, bodiesSimulation[bodyIndex].radius)
            for bodyIndex::Int64 in 1 : length(bodiesSimulation)];
        timesSimulation = libYukiPhysicsGravitationalNBodySimulation(libYukiConstantZero, timeNext - timeSimulating, timeStep, bodiesSimulation, integrator, gravitationalConstant) .+ timeSimulating;
        timeSimulating = timeNext + timeStep;
        timeNext = timeSimulating + (splitSteps * timeStep);
        function libYukiPhysicsGravitationalNBodySimulationStreamingAsyncSave!(timesSaving::Vector{Measurement{Float64}}, bodiesSaving::Vector{libYukiPhysicsBody}, savingDirectory::String, splitIndex::Int64)
            times::Vector{Measurement{Float64}} = timesSaving;
            bodies::Vector{libYukiPhysicsBody} = bodiesSaving;
            @save "$savingDirectory/libYukiPhysicsGravitationalNBodySimulationStreamingResult_$splitIndex.jld2" times bodies
        end
        libYukiPhysicsGravitationalNBodySimulationStreamingAsyncSave!(timesSimulation, bodiesSimulation, savingDirectory, splitIndex);
        splitIndex = splitIndex + 1;
    end
end

# Gravitational N-body simulation. 
# Dependency: Measurements, libYukiMath.
# Example: True.
function libYukiPhysicsGravitationalNBodySimulation(timeStart::Measurement{Float64}, timeEnd::Measurement{Float64}, timeStep::Measurement{Float64}, bodies::Vector{libYukiPhysicsBody}, integrator, gravitationalConstant::Measurement{Float64})::Vector{Measurement{Float64}}

    dimension::Int64 = 3;
    bodiesNumber::Int64 = length(bodies);

    function libYukiPhysicsGravitationalNBodySimulationAcceleration!(∂velocity::Vector{Measurement{Float64}}, velocity::Vector{Measurement{Float64}}, position::Vector{Measurement{Float64}}, partical, time::Measurement{Float64})
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
    velocities::Vector{Measurement{Float64}} = [body.velocity[end][dimensionIndex] for body in bodies for dimensionIndex in 1 : dimension];
    positions::Vector{Measurement{Float64}} = [body.position[end][dimensionIndex] for body in bodies for dimensionIndex in 1 : dimension];

    NBodyProblem = SecondOrderODEProblem(
        libYukiPhysicsGravitationalNBodySimulationAcceleration!,
        velocities,
        positions,
        (timeStart, timeEnd),
        (masses = masses, )
    );
    NBodySolution = solve(NBodyProblem, integrator, dt = timeStep);

    simulatedSteps::Int64 = length(NBodySolution.u);
    simulatedVelocityMatrix = [reshape(NBodySolution.u[i][1 : bodiesNumber * dimension], dimension, bodiesNumber) for i in 1 : simulatedSteps];
    simulatedPositionMatrix = [reshape(NBodySolution.u[i][bodiesNumber * dimension + 1 : 2 * bodiesNumber * dimension], dimension, bodiesNumber) for i in 1 : simulatedSteps];

    for bodyIndex in 1 : bodiesNumber
        bodies[bodyIndex].position = [
            [simulatedPositionMatrix[tIndex][dimIndex, bodyIndex]
            for dimIndex in 1 : dimension] 
            for tIndex in 1 : simulatedSteps];
        bodies[bodyIndex].velocity = [
            [simulatedVelocityMatrix[tIndex][dimIndex, bodyIndex]
            for dimIndex in 1 : dimension] 
            for tIndex in 1 : simulatedSteps];
    end
 
    return NBodySolution.t;
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
    displacementVectorMod::Measurement{Float64} = libYukiMathVectorToDistance(objectDisplacement);
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
