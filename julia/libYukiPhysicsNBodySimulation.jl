using Measurements, DifferentialEquations, JLD2, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl");
include("libYukiPhysics.jl");

# # Derive angles between position and . 
# # Dependency: Measurements, JLD2.
# # TODO: Validate & Example.
# function libYukiPhysicsNBodySimulationAngleBetweenPositionVector(position::Vector{Vector{Measurement{Float64}}}, vectorA::Vector{Measurement{Float64}})::Vector{Measurement{Float64}}
# 	return map(x -> libYukiMathAngleBetweenVectors(x, vectorA), position);
# end

# Load streaming saved gravitational N-body simulation result. 
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingLoad(savingDirectory::String)
	splitIndex::Int64 = 1;
	timesLoaded = nothing;
	bodiesLoaded = nothing;
	while ispath("$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2")
		@load "$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2" times bodies
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
# Dependency: Measurements, libYukiMath, JLD2, Dates.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreaming(timeStart::Measurement{Float64}, timeEnd::Measurement{Float64}, timeStep::Measurement{Float64}, bodiesSimulation::Vector{libYukiPhysicsBody}, integrator, gravitationalConstant::Measurement{Float64}, savingDirectory::String, splitSteps::Int64)
	splitIndex::Int64 = 1;
	timesSimulation::Vector{Measurement{Float64}} = [];

	totalSteps::Int64 = Int64(Measurements.value(ceil((timeEnd - timeStart) / timeStep)));
	totalGroups::Int64 = ceil(totalSteps / splitSteps);
	println("#INFO:[" * string(now()) * "] Simulation Groups: $totalGroups");

	timeSimulating::Measurement{Float64} = timeStart;
	timeNext::Measurement{Float64} = timeSimulating + splitSteps * timeStep;
	mkdir(savingDirectory);
	while timeSimulating < timeEnd
		println("#INFO:[" * string(now()) * "] Simulating Group $splitIndex / $totalGroups");
		timeNext = (timeNext < timeEnd) ? timeNext : timeEnd;
		bodiesSimulation = [ 
			libYukiPhysicsBody(bodiesSimulation[bodyIndex].name, bodiesSimulation[bodyIndex].position[end], bodiesSimulation[bodyIndex].velocity[end], bodiesSimulation[bodyIndex].mass, bodiesSimulation[bodyIndex].charge, bodiesSimulation[bodyIndex].radius)
			for bodyIndex::Int64 in 1 : length(bodiesSimulation)];
		timesSimulation = libYukiPhysicsNBodySimulationGravitational(libYukiConstantZero, timeNext - timeSimulating, timeStep, bodiesSimulation, integrator, gravitationalConstant) .+ timeSimulating;
		timeSimulating = timeNext + timeStep;
		timeNext = timeSimulating + (splitSteps * timeStep);
		function libYukiPhysicsNBodySimulationGravitationalStreamingSave!(timesSaving::Vector{Measurement{Float64}}, bodiesSaving::Vector{libYukiPhysicsBody}, savingDirectory::String, splitIndex::Int64)
			times::Vector{Measurement{Float64}} = timesSaving;
			bodies::Vector{libYukiPhysicsBody} = bodiesSaving;
			@save "$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2" times bodies
		end
		println("#INFO:[" * string(now()) * "]\t Group $splitIndex finished, storaging...");
		libYukiPhysicsNBodySimulationGravitationalStreamingSave!(timesSimulation, bodiesSimulation, savingDirectory, splitIndex);
		println("#INFO:[" * string(now()) * "]\t Progress " * string(round((splitIndex / totalGroups) * 100, digits = 2)) * "%.")
		splitIndex = splitIndex + 1;
	end
	println("#INFO:[" * string(now()) * "] Simulation Finished.");
end

# Gravitational N-body simulation. 
# Dependency: Measurements, libYukiMath.
# Example: True.
function libYukiPhysicsNBodySimulationGravitational(timeStart::Measurement{Float64}, timeEnd::Measurement{Float64}, timeStep::Measurement{Float64}, bodies::Vector{libYukiPhysicsBody}, integrator, gravitationalConstant::Measurement{Float64})::Vector{Measurement{Float64}}

	dimension::Int64 = 3;
	bodiesNumber::Int64 = length(bodies);

	function libYukiPhysicsNBodySimulationGravitationalAcceleration!(∂velocity::Vector{Measurement{Float64}}, velocity::Vector{Measurement{Float64}}, position::Vector{Measurement{Float64}}, partical, time::Measurement{Float64})
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
		libYukiPhysicsNBodySimulationGravitationalAcceleration!,
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
