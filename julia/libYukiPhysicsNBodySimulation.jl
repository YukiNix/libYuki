using Measurements, DifferentialEquations, JLD2, Dates

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl");
include("libYukiPhysics.jl");

# Calculate N-Body hamiltonian. 
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalHamiltonian(bodiesSimulation::Vector{libYukiPhysicsBody}, gravitationalConstant::Measurement{Float64})::Vector{Measurement{Float64}}
	simulationStep::Int64 = length(bodiesSimulation[1].velocity);
	for bodyLengthChecking::libYukiPhysicsBody in bodiesSimulation
		simulationVelocitySteps::Int64 = length(bodyLengthChecking.velocity);
		(simulationVelocitySteps == length(bodyLengthChecking.position) && simulationStep == simulationVelocitySteps) || error("Bodies simulation steps are not equal.")
		simulationStep = simulationVelocitySteps;
	end

	bodiesNumber::Int64 = length(bodiesSimulation);
	kineticEnergy::Vector{Measurement{Float64}} = zeros(Measurement{Float64}, simulationStep);
	potentialEnergy::Vector{Measurement{Float64}} = zeros(Measurement{Float64}, simulationStep);
	for stepIndex::Int64 in 1 : simulationStep
		kineticEnergy[stepIndex] = sum(
			libYukiPhysicsKineticEnergy(body.mass, body.velocity[stepIndex]) for body in bodiesSimulation
		);
		for bodyAIndex::Int64 in 1 : (bodiesNumber - 1)
			for bodyBIndex::Int64 in (bodyAIndex + 1) : bodiesNumber
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

# Find angular crossing time by specified vector. 
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationTimeRelativeAngleZeroCrossing(vectorA::Vector{Measurement{Float64}}, time::Vector{Measurement{Float64}}, position::Vector{Vector{Measurement{Float64}}})::Union{Vector{Measurement{Float64}}, Missing}
	sortedTime::Vector{Measurement{Float64}}, sortedPosition::Vector{Vector{Measurement{Float64}}} = libYukiBasicSortElementsByOrderVector(time, position);
	angles::Vector{Measurement{Float64}} = map(x -> libYukiMathAngleBetweenVector(vectorA, x), sortedPosition);

	crossIndices = findall(x -> angles[x - 1] > angles[x] && angles[x + 1] > angles[x], 2 : length(angles) - 2);
	if length(crossIndices) > 0
		k0_0::Vector{Measurement{Float64}} = (angles[crossIndices] .- angles[crossIndices .- 1]) ./ (sortedTime[crossIndices] .- sortedTime[crossIndices .- 1]);
		t0_0::Vector{Measurement{Float64}} = ((angles[crossIndices] .* sortedTime[crossIndices .- 1]) .- (angles[crossIndices .- 1] .* sortedTime[crossIndices])) ./ (angles[crossIndices] .- angles[crossIndices .- 1]);
		nk0_1::Vector{Measurement{Float64}} = -angles[crossIndices .+ 1] ./ (sortedTime[crossIndices .+ 1] .- t0_0);
		dk0mul::Vector{Measurement{Float64}} = k0_0 .* nk0_1;
		dk0plus::Vector{Measurement{Float64}} = abs.(k0_0 .+ nk0_1);

		k1_0::Vector{Measurement{Float64}} = (angles[crossIndices .+ 1] .- angles[crossIndices]) ./ (sortedTime[crossIndices .+ 1] .- sortedTime[crossIndices]);
		t1_0::Vector{Measurement{Float64}} = ((angles[crossIndices .+ 1] .* sortedTime[crossIndices]) .- (angles[crossIndices] .* sortedTime[crossIndices .+ 1])) ./ (angles[crossIndices .+ 1] .- angles[crossIndices]);
		nk1_1::Vector{Measurement{Float64}} = -angles[crossIndices .- 1] ./ (sortedTime[crossIndices .- 1] .- t1_0);
		dk1mul::Vector{Measurement{Float64}} = k1_0 .* nk1_1;
		dk1plus::Vector{Measurement{Float64}} = abs.(k1_0 .+ nk1_1);

		t0::Vector{Measurement{Float64}} = map(x -> 
				dk0mul[x] > 0 ? t1_0[x] : (
				dk1mul[x] > 0 ? t0_0[x] : (
				dk0plus[x] < dk1plus[x] ? t0_0[x] : t1_0[x]
			)),
			1 : length(k0_0));
		return t0;
	end
	return missing;
end

# Gravitational N-body simulation (stream save to file). 
# Dependency: Measurements, JLD2, Dates.
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

# Gravitational N-body simulation (stream follow a file). 
# Dependency: Measurements, JLD2, Dates.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingFollowSplit(timeSimulating::Measurement{Float64}, timeStep::Measurement{Float64}, integrator, gravitationalConstant::Measurement{Float64}, savingDirectory::String, followingSplit::Int64)
	timesSimulation = nothing;
	bodiesSimulation = nothing;
	times, bodies = libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory, followingSplit);

	if !isnothing(times)
		bodiesSimulation = [
			libYukiPhysicsBody(bodies[bodyIndex].name, bodies[bodyIndex].position[end], bodies[bodyIndex].velocity[end], bodies[bodyIndex].mass, bodies[bodyIndex].charge, bodies[bodyIndex].radius)
			for bodyIndex::Int64 in eachindex(bodies)];
		timesSimulation = libYukiPhysicsNBodySimulationGravitational(libYukiConstantZero, timeSimulating, timeStep, bodiesSimulation, integrator, gravitationalConstant) .+ timeSimulating;
	else
		return timesSimulation, bodiesSimulation;
	end
	return timesSimulation, bodiesSimulation;
end

# Load streaming saved gravitational N-body simulation result. 
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingLoad(savingDirectory::String)
	splitIndex::Int64 = 1;
	timesLoaded = nothing;
	bodiesLoaded = nothing;
	finishedFlag::Bool = false;
	while !finishedFlag
		times, bodies = libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory, splitIndex);
		if !isnothing(timesLoaded)
			append!(timesLoaded, times);
			for bodiesIndex::Int64 in eachindex(bodies)
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
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiPhysicsNBodySimulationGravitationalStreamingPartLoad(savingDirectory::String, splitIndex::Int64)
	times = nothing;
	bodies = nothing;
	if ispath("$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2")
		@load "$savingDirectory/libYukiPhysicsNBodySimulationGravitationalStreamingResult_$splitIndex.jld2" times bodies
	end
	return times, bodies;
end

# Gravitational N-body simulation. 
# Dependency: Measurements.
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
