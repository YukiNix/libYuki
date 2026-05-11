using Measurements, DifferentialEquations, JLD2, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl");
include("libYukiPhysics.jl");

# Calculate N-Body hamiltonian. 
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalHamiltonian(bodiesSimulation::AbstractVector{libYukiPhysicsBody{<:Real}}, gravitationalConstant::Real)
	simulationStep = length(bodiesSimulation[1].velocity);
	for bodyLengthChecking in bodiesSimulation
		simulationVelocitySteps = length(bodyLengthChecking.velocity);
		(simulationVelocitySteps == length(bodyLengthChecking.position) && simulationStep == simulationVelocitySteps) || error("Bodies simulation steps are not equal.")
		simulationStep = simulationVelocitySteps;
	end

	bodiesNumber = length(bodiesSimulation);
	kineticEnergy = zeros(T, simulationStep);
	potentialEnergy = zeros(T, simulationStep);
	for stepIndex in 1 : simulationStep
		kineticEnergy[stepIndex] = sum(
			libYukiPhysicsKineticEnergy(body.mass, body.velocity[stepIndex]) for body in bodiesSimulation
		);
		for bodyAIndex in 1 : (bodiesNumber - 1)
			for bodyBIndex in (bodyAIndex + 1) : bodiesNumber
				potentialEnergy[stepIndex] = potentialEnergy[stepIndex] + bodiesSimulation[bodyBIndex].mass * libYukiPhysicsGravitationalPotential(
					gravitationalConstant,
					bodiesSimulation[bodyAIndex].mass,
					libYukiMathVectorQuantity(bodiesSimulation[bodyBIndex].position[stepIndex] .- bodiesSimulation[bodyAIndex].position[stepIndex])
				);
			end
		end
	end
	
    return kineticEnergy .+ potentialEnergy;
end

# # Find angular crossing time by specified vector. 
# # Dependency: Measurements.
# # TODO: Validate & Example.
# function libYukiPhysicsNBodySimulationTimeRelativeAngleZeroCrossing(vectorA::AbstractVector{Measurement{Float64}}, time::AbstractVector{Measurement{Float64}}, position::AbstractVector{AbstractVector{Measurement{Float64}}})::Union{AbstractVector{Measurement{Float64}}, Missing}
# 	sortedTime::AbstractVector{Measurement{Float64}}, sortedPosition::AbstractVector{AbstractVector{Measurement{Float64}}} = libYukiBasicSortElementsByOrderVector(time, position);
# 	angles::AbstractVector{Measurement{Float64}} = map(x -> libYukiMathAngleBetweenVector(vectorA, x), sortedPosition);

# 	crossIndices = findall(x -> angles[x - 1] > angles[x] && angles[x + 1] > angles[x], 2 : length(angles) - 2);
# 	if length(crossIndices) > 0
# 		k0_0::AbstractVector{Measurement{Float64}} = (angles[crossIndices] .- angles[crossIndices .- 1]) ./ (sortedTime[crossIndices] .- sortedTime[crossIndices .- 1]);
# 		t0_0::AbstractVector{Measurement{Float64}} = ((angles[crossIndices] .* sortedTime[crossIndices .- 1]) .- (angles[crossIndices .- 1] .* sortedTime[crossIndices])) ./ (angles[crossIndices] .- angles[crossIndices .- 1]);
# 		nk0_1::AbstractVector{Measurement{Float64}} = -angles[crossIndices .+ 1] ./ (sortedTime[crossIndices .+ 1] .- t0_0);
# 		dk0mul::AbstractVector{Measurement{Float64}} = k0_0 .* nk0_1;
# 		dk0plus::AbstractVector{Measurement{Float64}} = abs.(k0_0 .+ nk0_1);

# 		k1_0::AbstractVector{Measurement{Float64}} = (angles[crossIndices .+ 1] .- angles[crossIndices]) ./ (sortedTime[crossIndices .+ 1] .- sortedTime[crossIndices]);
# 		t1_0::AbstractVector{Measurement{Float64}} = ((angles[crossIndices .+ 1] .* sortedTime[crossIndices]) .- (angles[crossIndices] .* sortedTime[crossIndices .+ 1])) ./ (angles[crossIndices .+ 1] .- angles[crossIndices]);
# 		nk1_1::AbstractVector{Measurement{Float64}} = -angles[crossIndices .- 1] ./ (sortedTime[crossIndices .- 1] .- t1_0);
# 		dk1mul::AbstractVector{Measurement{Float64}} = k1_0 .* nk1_1;
# 		dk1plus::AbstractVector{Measurement{Float64}} = abs.(k1_0 .+ nk1_1);

# 		t0::AbstractVector{Measurement{Float64}} = map(x -> 
# 				dk0mul[x] > 0 ? t1_0[x] : (
# 				dk1mul[x] > 0 ? t0_0[x] : (
# 				dk0plus[x] < dk1plus[x] ? t0_0[x] : t1_0[x]
# 			)),
# 			1 : length(k0_0));
# 		return t0;
# 	end
# 	return missing;
# end

# Gravitational N-body simulation (stream save to file). 
# Dependency: JLD2, Dates.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreaming(timeStart::Real, timeEnd::Real, timeStep::Real, bodiesSimulation::AbstractVector{libYukiPhysicsBody{<:Real}}, integrator, gravitationalConstant::Real, savingDirectory::String, splitSteps::Integer)
	splitIndex = 1;
	timesSimulation = T[];

	totalSteps = Int(ceil(Measurements.value((timeEnd - timeStart) / timeStep)));
	totalGroups = ceil(Int, totalSteps / splitSteps);
	println("#INFO:[" * string(now()) * "] Simulation Groups: $totalGroups");

	timeSimulating = timeStart;
	timeNext = timeSimulating + splitSteps * timeStep;
	mkdir(savingDirectory);
	while timeSimulating < timeEnd
		println("#INFO:[" * string(now()) * "] Simulating Group $splitIndex / $totalGroups");
		timeNext = (timeNext < timeEnd) ? timeNext : timeEnd;
		bodiesSimulation = [ 
			libYukiPhysicsBody(bodiesSimulation[bodyIndex].name, bodiesSimulation[bodyIndex].position[end], bodiesSimulation[bodyIndex].velocity[end], bodiesSimulation[bodyIndex].mass, bodiesSimulation[bodyIndex].charge, bodiesSimulation[bodyIndex].radius)
			for bodyIndex in 1 : length(bodiesSimulation)];
		timesSimulation = libYukiPhysicsNBodySimulationGravitational(zero(T), timeNext - timeSimulating, timeStep, bodiesSimulation, integrator, gravitationalConstant) .+ timeSimulating;
		timeSimulating = timeNext + timeStep;
		timeNext = timeSimulating + (splitSteps * timeStep);
		function libYukiPhysicsNBodySimulationGravitationalStreamingSave!(timesSaving, bodiesSaving, savingDirectory::String, splitIndex)
			times = timesSaving;
			bodies = bodiesSaving;
			@save "$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2" times bodies
		end
		println("#INFO:[" * string(now()) * "]\t Group $splitIndex finished, storaging...");
		libYukiPhysicsNBodySimulationGravitationalStreamingSave!(timesSimulation, bodiesSimulation, savingDirectory, splitIndex);
		println("#INFO:[" * string(now()) * "]\t Progress " * string(round((splitIndex / totalGroups) * 100, digits = 2)) * "%.")
		splitIndex = splitIndex + 1;
	end
	println("#INFO:[" * string(now()) * "] Simulation Finished.");
end

# Gravitational N-body simulation (stream follow a file). 
# Dependency: JLD2, Dates.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingFollowSplit(timeSimulating::Real, timeStep::Real, integrator, gravitationalConstant::Real, savingDirectory::String, followingSplit::Integer)
	timesSimulation = nothing;
	bodiesSimulation = nothing;
	times, bodies = libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory, followingSplit);

	if !isnothing(times)
		bodiesSimulation = [
			libYukiPhysicsBody(bodies[bodyIndex].name, bodies[bodyIndex].position[end], bodies[bodyIndex].velocity[end], bodies[bodyIndex].mass, bodies[bodyIndex].charge, bodies[bodyIndex].radius)
			for bodyIndex in eachindex(bodies)];
		timesSimulation = libYukiPhysicsNBodySimulationGravitational(zero(T), timeSimulating, timeStep, bodiesSimulation, integrator, gravitationalConstant) .+ timeSimulating;
	else
		return timesSimulation, bodiesSimulation;
	end
	return timesSimulation, bodiesSimulation;
end

# Load streaming saved gravitational N-body simulation result. 
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingLoad(savingDirectory::String)
	splitIndex = 1;
	timesLoaded = nothing;
	bodiesLoaded = nothing;
	finishedFlag = false;
	while !finishedFlag
		times, bodies = libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory, splitIndex);
		if !isnothing(timesLoaded)
			append!(timesLoaded, times);
			for bodiesIndex in eachindex(bodies)
				append!(bodiesLoaded[bodiesIndex].velocity, bodies[bodiesIndex].velocity);
				append!(bodiesLoaded[bodiesIndex].position, bodies[bodiesIndex].position);
			end
		else
			if splitIndex == 1
				timesLoaded = times;
				bodiesLoaded = bodies;
			else
				finishedFlag = true;
			end
		end
		splitIndex += 1;
	end
	return timesLoaded, bodiesLoaded;
end

# Load streaming saved gravitational N-body simulation result partially. 
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory::String, splitIndex::Integer)
	times = nothing;
	bodies = nothing;
	if ispath("$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2")
		@load "$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2" times bodies
	end
	return times, bodies;
end

# Gravitational N-body simulation. 
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitational(timeStart::Real, timeEnd::Real, timeStep::Real, bodies::AbstractVector{libYukiPhysicsBody{<:Real}}, integrator, gravitationalConstant::Real)

	dimension = 3;
	bodiesNumber = length(bodies);

	function libYukiPhysicsNBodySimulationGravitationalAcceleration!(∂velocity, velocity, position, partical, time)
		separatedPositions = reshape(position, dimension, bodiesNumber);
		∂velocity .= zero(T);
		for bodyAIndex in 1 : bodiesNumber
			bodyAAcceleration = zeros(T, dimension);
			for bodyBIndex in 1 : bodiesNumber
				if bodyAIndex != bodyBIndex
					displacementBA = separatedPositions[:, bodyAIndex] - separatedPositions[:, bodyBIndex];
					bodyAAcceleration += libYukiPhysicsForceAcceleration(
						libYukiPhysicsPotentialEnergyForce(
							x -> Measurements.value(partical.masses[bodyAIndex]) * libYukiPhysicsGravitationalPotential(
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

	masses = map(x -> x.mass, bodies);
	velocities = [body.velocity[end][dimensionIndex] for body in bodies for dimensionIndex in 1 : dimension];
	positions = [body.position[end][dimensionIndex] for body in bodies for dimensionIndex in 1 : dimension];

	NBodyProblem = SecondOrderODEProblem(
		libYukiPhysicsNBodySimulationGravitationalAcceleration!,
		velocities,
		positions,
		(timeStart, timeEnd),
		(masses = masses, )
	);
	NBodySolution = solve(NBodyProblem, integrator, dt = timeStep);

	simulatedSteps = length(NBodySolution.u);
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
