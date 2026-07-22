using Optim

include("libYukiAstronomyTransitLimbDarkening.jl")

"""
    libYukiAstronomyTransitTTV(;
        oc::Vector{Float64} = Float64[], 
        ocErr::Vector{Float64} = Float64[], 
        isAcceptable::Vector{Bool} = Bool[], 
        predictedMidPoints::Vector{Float64} = Float64[]
    )
Define the transit timing variations (TTV) of a planet.
# Arguments
- `oc`: A vector of observed minus calculated (O-C) transit times.
- `ocErr`: A vector of errors associated with the O-C transit times.
- `isAcceptable`: A vector of boolean values indicating whether each
transit timing measurement is acceptable.
- `predictedMidPoints`: A vector of predicted mid-transit times based
on the planet's orbital period and initial transit time.
# Returns
- A `libYukiAstronomyTransitTTV` instance containing the TTV data.
"""
mutable struct libYukiAstronomyTransitTTV
    oc::AbstractVector{<:Real}
    ocErr::AbstractVector{<:Real}
    isAcceptable::AbstractVector{Bool}
    predictedMidPoints::AbstractVector{<:Real}
    libYukiAstronomyTransitTTV(;
        oc::AbstractVector{<:Real} = Float64[], 
		ocErr::AbstractVector{<:Real} = Float64[], 
		isAcceptable::AbstractVector{Bool} = Bool[], 
		predictedMidPoints::AbstractVector{<:Real} = Float64[]
	) = new(
		oc, 
		ocErr, 
		isAcceptable, 
		predictedMidPoints
	);
end

"""
	libYukiAstronomyTransitTTVCorrectLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve,
		ttv::libYukiAstronomyTransitTTV,
		transit::libYukiAstronomyTransit
	) 
    libYukiAstronomyTransitTTVCorrectLightCurve(
        lightCurve::libYukiAstronomyTransitLightCurve, 
        ttv::libYukiAstronomyTransitTTV,
        planetOrbitPeriod::Real
    )
	libYukiAstronomyTransitTTVCorrectLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve,
		ttv::libYukiAstronomyTransitTTV,
		transit::libYukiAstronomyTransit
	)
    libYukiAstronomyTransitTTVCorrectLightCurve!(
        lightCurve::libYukiAstronomyTransitLightCurve, 
        ttv::libYukiAstronomyTransitTTV,
        planetOrbitPeriod::Real
    )
Correct the timing of a light curve based on the transit timing 
variations (TTV) of a planet. The function can be called with either 
a `libYukiAstronomyTransitLightCurve` instance, in which case a new 
corrected light curve will be returned, or with a 
`libYukiAstronomyTransitLightCurve` instance and the `!` suffix, 
in which case the input light curve will be modified in place.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing time values and flux values to be corrected.
- `ttv`: An instance of `libYukiAstronomyTransitTTV` containing the 
transit timing variations to be applied.
- `transit`: An instance of `libYukiAstronomyTransit` containing the
planet's orbital period, which is used to determine the time window
for applying the TTV corrections.
- `planetOrbitPeriod`: The orbital period of the planet, used to 
determine the time window for applying the TTV corrections.
# Returns
- A new `libYukiAstronomyTransitLightCurve` instance with corrected 
time values,
- Or modifies the input `libYukiAstronomyTransitLightCurve` instance 
in place if the `!` suffix is used.
"""
libYukiAstronomyTransitTTVCorrectLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve,
	ttv::libYukiAstronomyTransitTTV,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitTTVCorrectLightCurve(
	lightCurve,
	ttv,
	transit.planet.orbit.period
);
function libYukiAstronomyTransitTTVCorrectLightCurve(
    lightCurve::libYukiAstronomyTransitLightCurve,
    ttv::libYukiAstronomyTransitTTV,
    planetOrbitPeriod::Real
)
    lightCurveCorrected = deepcopy(lightCurve);
    libYukiAstronomyTransitTTVCorrectLightCurve!(
        lightCurveCorrected,
        ttv,
        planetOrbitPeriod,
    );
    return lightCurveCorrected;
end
libYukiAstronomyTransitTTVCorrectLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	ttv::libYukiAstronomyTransitTTV,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitTTVCorrectLightCurve!(
	lightCurve,
	ttv,
	transit.planet.orbit.period
);
function libYukiAstronomyTransitTTVCorrectLightCurve!(
    lightCurve::libYukiAstronomyTransitLightCurve,
    ttv::libYukiAstronomyTransitTTV,
    planetOrbitPeriod::Real
)
    originalTime = copy(lightCurve.time);
    for (index, transitMidTime) in enumerate(ttv.predictedMidPoints)
        mask = abs.(originalTime .- transitMidTime) .< 
            0.5 * planetOrbitPeriod;
        lightCurve.time[mask] .-= ttv.oc[index];
    end
    libYukiAstronomyTransitExtractValidLightCurve!(lightCurve);
    return lightCurve;
end

"""
	libYukiAstronomyTransitTTVOCBrent(
		lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
		transit::libYukiAstronomyTransit,
		limbDarkeningFunc::AbstractLimbDark;
		ocMeasurementLengthDurations::Real = 1.0,
		ocSearchLengthDurations::Real = 0.1,
		minimumPointCount::Integer = 5,
		minimumPointsPerSide::Integer = 2,
		maximumTimingErrorDurations::Real = 0.1,
		profilePointCount::Integer = 101,
		searchCentreFitPointCount::Integer = 5,
		searchWidthMaximumMultiplier::Real = 8.0
	)
	libYukiAstronomyTransitTTVOCBrent(
		lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
		transit::libYukiAstronomyTransit,
		orbit::Orbits.AbstractOrbit,
		limbDarkeningFunc::AbstractLimbDark;
		ocMeasurementLengthDurations::Real = 1.0,
		ocSearchLengthDurations::Real = 0.1,
		minimumPointCount::Integer = 5,
		minimumPointsPerSide::Integer = 2,
		maximumTimingErrorDurations::Real = 0.1,
		profilePointCount::Integer = 101,
		searchCentreFitPointCount::Integer = 5,
		searchWidthMaximumMultiplier::Real = 8.0
	)
    libYukiAstronomyTransitTTVOCBrent(
        lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
		planetOrbitPeriod::Real,
		planetTransitDuration::Real,
		planetTransitCentreTime::Real,
		planetStellarRadiusRatio::Real,
		orbit::Orbits.AbstractOrbit,
		limbDarkeningFunc::AbstractLimbDark;
		ocMeasurementLengthDurations::Real = 1.0,
		ocSearchLengthDurations::Real = 0.1,
		minimumPointCount::Integer = 5,
		minimumPointsPerSide::Integer = 2,
		maximumTimingErrorDurations::Real = 0.1,
		profilePointCount::Integer = 101,
		searchCentreFitPointCount::Integer = 5,
		searchWidthMaximumMultiplier::Real = 8.0,
    )
Calculate the transit timing variations (TTV) using the Brent 
optimization method.
# Arguments
- `lightCurvesSplited`: A vector of `libYukiAstronomyTransitLightCurve` instances, each representing a light curve segment for a transit.
- `transit`: An instance of `libYukiAstronomyTransit` containing the planet's orbital parameters.
- `orbit`: An instance of `Orbits.AbstractOrbit` representing the planet's orbit.
- `limbDarkeningFunc`: An instance of `AbstractLimbDark` representing the limb darkening function.
- `ocMeasurementLengthDurations`: The length of the measurement window in units of the planet's transit duration.
- `ocSearchLengthDurations`: The length of the search window in units of the planet's transit duration.
- `minimumPointCount`: The minimum number of data points required in the measurement window to perform the TTV calculation.
- `minimumPointsPerSide`: The minimum number of data points required on each side of the transit center to perform the TTV calculation.
- `maximumTimingErrorDurations`: The maximum allowed timing error in units of the planet's transit duration.
- `profilePointCount`: The number of points to use in the profile likelihood calculation.
- `searchCentreFitPointCount`: The number of previous accepted TTV measurements to use for fitting the search center.
- `searchWidthMaximumMultiplier`: The maximum multiplier for the search width in case of boundary hits or incomplete profiles.
# Returns
- A vector of calculated transit timing variations (TTV) for each 
transit in the input light curves.
"""
libYukiAstronomyTransitTTVOCBrent(
	lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
	transit::libYukiAstronomyTransit,
	limbDarkeningFunc::AbstractLimbDark;
	ocMeasurementLengthDurations::Real = 1.0,
	ocSearchLengthDurations::Real = 0.1,
	minimumPointCount::Integer = 5,
	minimumPointsPerSide::Integer = 2,
	maximumTimingErrorDurations::Real = 0.1,
	profilePointCount::Integer = 101,
	searchCentreFitPointCount::Integer = 5,
	searchWidthMaximumMultiplier::Real = 8.0
) = libYukiAstronomyTransitTTVOCBrent(
	lightCurvesSplited,
	transit,
	SimpleOrbit(
        period = transit.planet.orbit.period, 
        duration = transit.planetTransitDuration
    ),
	limbDarkeningFunc;
	ocMeasurementLengthDurations = ocMeasurementLengthDurations,
	ocSearchLengthDurations = ocSearchLengthDurations,
	minimumPointCount = minimumPointCount,
	minimumPointsPerSide = minimumPointsPerSide,
	maximumTimingErrorDurations = maximumTimingErrorDurations,
	profilePointCount = profilePointCount,
	searchCentreFitPointCount = searchCentreFitPointCount,
	searchWidthMaximumMultiplier = searchWidthMaximumMultiplier
);
libYukiAstronomyTransitTTVOCBrent(
	lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
	transit::libYukiAstronomyTransit,
	orbit::Orbits.AbstractOrbit,
	limbDarkeningFunc::AbstractLimbDark;
	ocMeasurementLengthDurations::Real = 1.0,
	ocSearchLengthDurations::Real = 0.1,
	minimumPointCount::Integer = 5,
	minimumPointsPerSide::Integer = 2,
	maximumTimingErrorDurations::Real = 0.1,
	profilePointCount::Integer = 101,
	searchCentreFitPointCount::Integer = 5,
	searchWidthMaximumMultiplier::Real = 8.0
) = libYukiAstronomyTransitTTVOCBrent(
	lightCurvesSplited,
	transit.planet.orbit.period,
	transit.planetTransitDuration,
	transit.planetTransitCentreTime,
	transit.planetStellarRadiusRatio,
	orbit,
	limbDarkeningFunc;
	ocMeasurementLengthDurations = ocMeasurementLengthDurations,
	ocSearchLengthDurations = ocSearchLengthDurations,
	minimumPointCount = minimumPointCount,
	minimumPointsPerSide = minimumPointsPerSide,
	maximumTimingErrorDurations = maximumTimingErrorDurations,
	profilePointCount = profilePointCount,
	searchCentreFitPointCount = searchCentreFitPointCount,
	searchWidthMaximumMultiplier = searchWidthMaximumMultiplier
);
function libYukiAstronomyTransitTTVOCBrent(
	lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
	planetOrbitPeriod::Real,
	planetTransitDuration::Real,
	planetTransitCentreTime::Real,
	planetStellarRadiusRatio::Real,
	orbit::Orbits.AbstractOrbit,
	limbDarkeningFunc::AbstractLimbDark;
	ocMeasurementLengthDurations::Real = 1.0,
	ocSearchLengthDurations::Real = 0.1,
	minimumPointCount::Integer = 5,
	minimumPointsPerSide::Integer = 2,
	maximumTimingErrorDurations::Real = 0.1,
	profilePointCount::Integer = 101,
	searchCentreFitPointCount::Integer = 5,
	searchWidthMaximumMultiplier::Real = 8.0,
)
	nTransits = length(lightCurvesSplited);

	oc = fill(NaN, nTransits);
	ocErr = fill(NaN, nTransits);
	isAcceptable = fill(false, nTransits);
	epochs = collect(0 : nTransits - 1);

	predictedMidPoints = planetTransitCentreTime .+ 
		epochs .* planetOrbitPeriod;

	searchWidth = ocSearchLengthDurations * planetTransitDuration;

	for (index, t0) in enumerate(predictedMidPoints)
		lightCurve = lightCurvesSplited[index];

		searchCentre = 0.0;
		acceptedIndices = index > 1 ?
			findall(isAcceptable[1 : index - 1]) :
			Int[];
		
		if length(acceptedIndices) == 1
			searchCentre = oc[only(acceptedIndices)];

		elseif length(acceptedIndices) >= 2
			firstFitIndex = max(
				1,
				length(acceptedIndices) - 
					searchCentreFitPointCount + 1,
			);
			fitIndices = acceptedIndices[firstFitIndex : end];
			fitEpoch = Float64.(epochs[fitIndices]);
			fitOC = oc[fitIndices];

			designMatrix = hcat(
				ones(length(fitEpoch)),
				fitEpoch,
			);
			driftParameters = designMatrix \ fitOC;

			searchCentre = driftParameters[1] +
				driftParameters[2] * epochs[index];
		end

		measurementCentre = t0 + searchCentre;
		mask = abs.(lightCurve.time .- measurementCentre) .<
			ocMeasurementLengthDurations * planetTransitDuration;
		count(mask) >= minimumPointCount ||
			continue

		t = lightCurve.time[mask] .- t0;
		f = lightCurve.flux[mask];
		all(isfinite, t) || continue
		all(isfinite, f) || continue

		fluxErrProvided = lightCurve.fluxErr !== nothing &&
			!isempty(lightCurve.fluxErr);
		if fluxErrProvided
			length(lightCurve.fluxErr) ==
				length(lightCurve.time) ||
				continue
			fErr = lightCurve.fluxErr[mask];
			all(isfinite, fErr) || continue
			all(fErr .> 0) || continue
			sqrtWeight = 1.0 ./ fErr;
		else
			sqrtWeight = ones(length(t));
		end

		function loss(dt)
			templateFlux = 
				libYukiAstronomyTransitLimbDarkeningTransitFlux(
					t .- dt,
					planetStellarRadiusRatio,
					orbit,
					limbDarkeningFunc,
				);
			all(isfinite, templateFlux) ||
				return Inf
			
			designMatrix = hcat(
				ones(length(t)),
				templateFlux .- 1.0,
			);
			weightedDesignMatrix = designMatrix .*
				reshape(sqrtWeight, :, 1);
			parameters = weightedDesignMatrix \ (f .* sqrtWeight);
			residual = (f .- designMatrix * parameters) .*
				sqrtWeight;

			return sum(abs2, residual);
		end

		widthMultiplier = 1.0;
		while widthMultiplier <= searchWidthMaximumMultiplier
			currentSearchWidth =
				widthMultiplier * searchWidth;
			searchLower = searchCentre - currentSearchWidth;
			searchUpper = searchCentre + currentSearchWidth;

			result = optimize(
				loss,
				searchLower,
				searchUpper,
				Brent(),
			);
			if !Optim.converged(result)
				widthMultiplier *= 2.0;
				continue;
			end

			dtBest = Optim.minimizer(result);
			lossBest = Optim.minimum(result);
			if !(isfinite(dtBest) && isfinite(lossBest))
				widthMultiplier *= 2.0;
				continue;
			end
			if !(abs(dtBest - searchCentre) <
				0.9 * currentSearchWidth)
				widthMultiplier *= 2.0;
				continue;
			end

			shiftedTime = t .- dtBest;
			count(shiftedTime .< 0) >=
				minimumPointsPerSide ||
				break
			count(shiftedTime .> 0) >=
				minimumPointsPerSide ||
				break

			degreesOfFreedom = length(t) - 3;
			degreesOfFreedom > 0 ||
				break
			profileLossIncrease = fluxErrProvided ?
				1.0 :
				lossBest / degreesOfFreedom;
			(isfinite(profileLossIncrease) &&
				profileLossIncrease > 0) ||
				break

			dtGrid = collect(
				range(
					searchLower,
					searchUpper;
					length = profilePointCount,
				),
			);
			lossGrid = loss.(dtGrid);
			targetLoss = lossBest + profileLossIncrease;
			
			leftIndices = findall(
				(dtGrid .< dtBest) .&
					(lossGrid .>= targetLoss),
			);
			rightIndices = findall(
				(dtGrid .> dtBest) .&
					(lossGrid .>= targetLoss),
			);
			if isempty(leftIndices) || isempty(rightIndices)
				widthMultiplier *= 2.0;
				continue;
			end

			dtLeft = dtGrid[last(leftIndices)];
			dtRight = dtGrid[first(rightIndices)];
			timingError = (dtRight - dtLeft) / 2;
			(isfinite(timingError) &&
				timingError > 0) ||
				break
			timingError <=
				maximumTimingErrorDurations * planetTransitDuration ||
				break

			oc[index] = dtBest;
			ocErr[index] = timingError;
			isAcceptable[index] = true;
			break;
		end
	end

	return libYukiAstronomyTransitTTV(;
		oc = oc,
		ocErr = ocErr,
		isAcceptable = isAcceptable,
		predictedMidPoints = predictedMidPoints,
	);
end