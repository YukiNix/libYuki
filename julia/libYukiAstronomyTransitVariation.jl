using Optim, JLD2

include("libYukiAstronomyTransitLimbDarkening.jl")

"""
    libYukiAstronomyTransitVariation(;
        ttvOC::Vector{Float64} = Float64[], 
        ttvOCErr::Vector{Float64} = Float64[], 
        ttvOCIsAcceptable::Vector{Bool} = Bool[], 
		tdvOR::Vector{Float64} = Float64[],
		tdvORErr::Vector{Float64} = Float64[],
		tdvORIsAcceptable::Vector{Bool} = Bool[],
		transitDurationRef::Float64 = NaN,
		transitMidPoints::Vector{Float64} = Float64[],
        predictedMidPoints::Vector{Float64} = Float64[]
    )
Define the transit timing variations (TTV) of a planet.
# Arguments
- `ttvOC`: A vector of observed minus calculated (Observed minus 
Calculated) transit times.
- `ttvOCErr`: A vector of errors associated with the O-C transit times.
- `ttvOCIsAcceptable`: A vector of boolean values indicating whether 
each transit timing measurement is acceptable.
- `tdvOR`: A vector of observed minus calculated (Observed minus 
Reference) transit durations.
- `tdvORErr`: A vector of errors associated with the O-R transit 
durations.
- `tdvORIsAcceptable`: A vector of boolean values indicating whether 
each transit duration measurement is acceptable.
- `transitDurationRef`: A reference transit duration used for 
calculating the O-R values.
- `transitMidPoints`: A vector of observed mid-transit times.
- `predictedMidPoints`: A vector of predicted mid-transit times based
on the planet's orbital period and initial transit time.
# Returns
- A `libYukiAstronomyTransitVariation` instance containing the TTV data.
"""
mutable struct libYukiAstronomyTransitVariation
    ttvOC::AbstractVector{<:Real}
    ttvOCErr::AbstractVector{<:Real}
    ttvOCIsAcceptable::AbstractVector{Bool}
	tdvOR::AbstractVector{<:Real}
	tdvORErr::AbstractVector{<:Real}
	tdvORIsAcceptable::AbstractVector{Bool}
	transitDurationRef::Real
	transitMidPoints::AbstractVector{<:Real}
    predictedMidPoints::AbstractVector{<:Real}
    libYukiAstronomyTransitVariation(;
        ttvOC::AbstractVector{<:Real} = Float64[], 
		ttvOCErr::AbstractVector{<:Real} = Float64[], 
		ttvOCIsAcceptable::AbstractVector{Bool} = Bool[],
		tdvOR::AbstractVector{<:Real} = Float64[],
		tdvORErr::AbstractVector{<:Real} = Float64[],
		tdvORIsAcceptable::AbstractVector{Bool} = Bool[],
		transitDurationRef::Real = NaN,
		transitMidPoints::AbstractVector{<:Real} = Float64[],
		predictedMidPoints::AbstractVector{<:Real} = Float64[]
	) = new(
		ttvOC, 
		ttvOCErr, 
		ttvOCIsAcceptable, 
		tdvOR,
		tdvORErr,
		tdvORIsAcceptable,
		transitDurationRef,
		transitMidPoints,
		predictedMidPoints
	);
end

""" 
	libYukiAstronomyTransitVariationTDVCorrectLightCurve(
		lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit,
		variation::libYukiAstronomyTransitVariation
	)
	libYukiAstronomyTransitVariationTDVCorrectLightCurve!(
		lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit,
		variation::libYukiAstronomyTransitVariation
	)
Correct the timing of a light curve based on the transit duration 
variations (TDV).
# Arguments
- `lightCurveTTVCorrected`: An instance of 
`libYukiAstronomyTransitLightCurve` containing the light curve data 
that has been corrected for transit timing variations (TTV).
- `transit`: An instance of `libYukiAstronomyTransit` containing 
the transit parameters.
- `variation`: An instance of `libYukiAstronomyTransitVariation` 
containing the transit duration variations to be applied.
# Returns
- A new `libYukiAstronomyTransitLightCurve` instance with corrected 
time values, or modifies the input `libYukiAstronomyTransitLightCurve` 
instance in place if the `!` suffix is used.
"""
function libYukiAstronomyTransitVariationTDVCorrectLightCurve(
	lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit,
	variation::libYukiAstronomyTransitVariation
)
	lightCurveTDVCorrected = deepcopy(lightCurveTTVCorrected);
	libYukiAstronomyTransitVariationTDVCorrectLightCurve!(
		lightCurveTDVCorrected,
		transit,
		variation
	);
	return lightCurveTDVCorrected;
end
function libYukiAstronomyTransitVariationTDVCorrectLightCurve!(
	lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit,
	variation::libYukiAstronomyTransitVariation
)
	originalTime = copy(lightCurveTTVCorrected.time);
    for (index, transitMidTime) in enumerate(variation.predictedMidPoints)
		mask = abs.(originalTime .- transitMidTime) .<
			0.5 * transit.planet.orbit.period;
		durationScale = 1 +
			(variation.tdvOR[index] / transit.planetTransitDuration);
		lightCurveTTVCorrected.time[mask] .= transitMidTime .+
			(originalTime[mask] .- transitMidTime) ./ durationScale;
	end
    libYukiAstronomyTransitExtractValidLightCurve!(lightCurveTTVCorrected);

    return lightCurveTTVCorrected;
end

"""
	libYukiAstronomyTransitVariationTDVORBrent(
		lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit,
		variation::libYukiAstronomyTransitVariation;
		searchFraction::Real = 0.1,
		profilePointCount::Integer = 101
	)
Measure the transit duration variations (TDV) of a planet using the
Brent optimization method. The function requires a light curve that has
been corrected for transit timing variations (TTV) and a reference
transit duration. If the reference transit duration is not provided, 
it will be measured from the folded light curve.
# Arguments
- `lightCurveTTVCorrected`: An instance of 
`libYukiAstronomyTransitLightCurve` containing the light curve data
that has been corrected for transit timing variations (TTV).
- `transit`: An instance of `libYukiAstronomyTransit` containing 
the transit parameters.
- `variation`: An instance of `libYukiAstronomyTransitVariation` to 
store the measured transit duration variations (TDV).
- `searchFraction`: A fraction of the initial transit duration to 
define the search range for optimization.
- `profilePointCount`: The number of points used to profile the 
loss function for error estimation.
# Returns
- The updated `libYukiAstronomyTransitVariation` instance 
containing the measured transit duration variations (TDV) and 
their associated errors.
"""
function libYukiAstronomyTransitVariationTDVORBrent(
	lightCurveTTVCorrected::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit,
	variation::libYukiAstronomyTransitVariation;
	searchFraction::Real = 0.1,
	profilePointCount::Integer = 101
)
	variation.tdvOR = fill(NaN, length(variation.ttvOC));
	variation.tdvORErr = fill(NaN, length(variation.ttvOC));
	variation.tdvORIsAcceptable = fill(false, length(variation.ttvOC));
	if isnan(variation.transitDurationRef)
		libYukiAstronomyTransitVariationDurationBrent(
			libYukiAstronomyTransitFoldLightCurve(
				lightCurveTTVCorrected,
				transit
			),
			transit,
			variation;
			searchFraction = searchFraction,
			updateTransitDuration = true,
			updateDurationRef = true,
			profilePointCount = profilePointCount
		);
	end
	transit.planetTransitDuration = variation.transitDurationRef;

	lightCurvesSplited = libYukiAstronomyTransitSplitLightCurveByPeriod(
		lightCurveTTVCorrected,
		transit
	);
	for (index, lightCurve) in enumerate(lightCurvesSplited)
		if !variation.ttvOCIsAcceptable[index]
			continue
		end
		durationMeasured, durationMeasuredErr = 
			libYukiAstronomyTransitVariationDurationBrent(
				lightCurve,
				transit;
				searchFraction = searchFraction,
				profilePointCount = profilePointCount,
				updateTransitDuration = false
			);
		variation.tdvOR[index] = durationMeasured - variation.transitDurationRef;
		variation.tdvORErr[index] = durationMeasuredErr;
		variation.tdvORIsAcceptable[index] = isfinite(durationMeasured) && isfinite(durationMeasuredErr);
	end

	return variation;
end

"""
	libYukiAstronomyTransitVariationDurationBrent(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit,
		variation::libYukiAstronomyTransitVariation;
		searchFraction::Real = 0.1,
		profilePointCount::Integer = 101,
		updateDurationRef::Bool = false,
		updateTransitDuration::Bool = false
	)
	libYukiAstronomyTransitVariationDurationBrent(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit;
		searchFraction::Real = 0.1,
		profilePointCount::Integer = 101,
		updateTransitDuration::Bool = false
	)
Measure the transit duration of a planet using the Brent 
optimization method. The light curve must be phase-folded with 
the transit centre at time zero, since the duration is probed by 
rescaling the time axis about zero.
# Arguments
- `lightCurve`: An instance of 
`libYukiAstronomyTransitLightCurve` containing the light 
curve data.
- `transit`: An instance of `libYukiAstronomyTransit` containing 
the transit parameters.
- `variation`: An instance of `libYukiAstronomyTransitVariation` 
to store the measured transit duration.
- `searchFraction`: A fraction of the initial transit duration to 
define the search range for optimization.
- `profilePointCount`: The number of points used to profile the loss 
function.
- `updateDurationRef`: A boolean indicating whether to update the 
reference transit duration in the `variation` object with the 
measured value.
- `updateTransitDuration`: A boolean indicating whether to update 
the transit duration in the `transit` object with the measured value.
# Returns
- A tuple `(durationRef, durationErr)` containing the measured 
transit duration and its estimated uncertainty (`NaN` if the 
uncertainty cannot be determined within the search range).
"""
function libYukiAstronomyTransitVariationDurationBrent(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit,
	variation::libYukiAstronomyTransitVariation;
	searchFraction::Real = 0.1,
	profilePointCount::Integer = 101,
	updateDurationRef::Bool = false,
	updateTransitDuration::Bool = false
)
	durationRef, durationErr = libYukiAstronomyTransitVariationDurationBrent(
		lightCurve,
		transit;
		searchFraction = searchFraction,
		profilePointCount = profilePointCount,
		updateTransitDuration = updateTransitDuration
	);
	if updateDurationRef 
		variation.transitDurationRef = durationRef;
	end
	return durationRef, durationErr;
end
function libYukiAstronomyTransitVariationDurationBrent(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit;
	searchFraction::Real = 0.1,
	profilePointCount::Integer = 101,
	updateTransitDuration::Bool = false
)
	durationInitial = transit.planetTransitDuration;
	hasFluxErr = lightCurve.fluxErr !== nothing &&
		!isempty(lightCurve.fluxErr);
	if hasFluxErr
		length(lightCurve.fluxErr) ==
			length(lightCurve.flux) ||
				throw(
					DimensionMismatch(
						"fluxErr and flux must have the same length.",
					),
				)
		all(isfinite, lightCurve.fluxErr) ||
			throw(
				ArgumentError(
					"fluxErr must contain only finite values.",
				),
			)
		all(lightCurve.fluxErr .> 0) ||
			throw(
				ArgumentError(
					"fluxErr must contain only positive values.",
				),
			)
	end

	function loss(duration)
		scaledTime = lightCurve.time .* 
			(durationInitial / duration);
		lightCurveScaled = libYukiAstronomyTransitLightCurve(
			time = scaledTime,
		);
		libYukiAstronomyTransitLimbDarkeningTransitFlux!(
			lightCurveScaled,
			transit,
		);
		designMatrix = hcat(
			ones(length(lightCurve.flux)),
			lightCurveScaled.flux .- 1,
		);
		if hasFluxErr
			weightedDesignMatrix =
				designMatrix ./ lightCurve.fluxErr;
			parameters = weightedDesignMatrix \
				(lightCurve.flux ./ lightCurve.fluxErr);
			residual = (
					designMatrix * parameters .- lightCurve.flux
				) ./ lightCurve.fluxErr;
		else
			parameters = designMatrix \ lightCurve.flux;
			residual = designMatrix * parameters .-
				lightCurve.flux;
		end
		return sum(abs2, residual);
	end

	durationLower = durationInitial * (1 - searchFraction);
	durationUpper = durationInitial * (1 + searchFraction);
	result = optimize(
		loss,
		durationLower,
		durationUpper,
		Brent(),
	);
	durationRef = Optim.minimizer(result);
	lossRef = Optim.minimum(result);
	degreesOfFreedom = length(lightCurve.flux) - 1;
	lossIncrease = hasFluxErr ?
		1.0 :
		lossRef / degreesOfFreedom;
	targetLoss = lossRef + lossIncrease;
	durationGrid = range(
		durationLower,
		durationUpper;
		length = profilePointCount
	);
	lossGrid = loss.(durationGrid);
	leftIndices = findall(
		(durationGrid .< durationRef) .&
		(lossGrid .>= targetLoss),
	);
	rightIndices = findall(
		(durationGrid .> durationRef) .&
		(lossGrid .>= targetLoss),
	);
	durationErr = (isempty(leftIndices) || isempty(rightIndices)) ?
		NaN :
		(
			durationGrid[first(rightIndices)] - 
				durationGrid[last(leftIndices)]
		) / 2;
	if updateTransitDuration
		transit.planetTransitDuration = durationRef;
	end

	return durationRef, durationErr;
end

"""
	libYukiAstronomyTransitVariationLoadVariation(
		filename::String
	)
Load the transit variations data from a file.
# Arguments
- `filename`: The name of the file to load the data from.
# Returns
- A `libYukiAstronomyTransitVariation` instance containing the 
transit variation data.
"""
function libYukiAstronomyTransitVariationLoadVariation(
	filename::String
)
	@load filename variation
	return variation;
end

"""
	libYukiAstronomyTransitVariationSaveVariation(
		variation::libYukiAstronomyTransitVariation, 
		filename::String
	)
Save the transit variations data to a file.
# Arguments
- `variation`: An instance of `libYukiAstronomyTransitVariation` containing the transit variation data.
- `filename`: The name of the file to save the data to.
"""
function libYukiAstronomyTransitVariationSaveVariation(
	variation::libYukiAstronomyTransitVariation, 
	filename::String
)
	@save filename variation
end

"""
	libYukiAstronomyTransitVariationTTVCorrectLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve,
		ttv::libYukiAstronomyTransitVariation,
		transit::libYukiAstronomyTransit
	) 
    libYukiAstronomyTransitVariationTTVCorrectLightCurve(
        lightCurve::libYukiAstronomyTransitLightCurve, 
        ttv::libYukiAstronomyTransitVariation,
        planetOrbitPeriod::Real
    )
	libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve,
		ttv::libYukiAstronomyTransitVariation,
		transit::libYukiAstronomyTransit
	)
    libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
        lightCurve::libYukiAstronomyTransitLightCurve, 
        ttv::libYukiAstronomyTransitVariation,
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
- `ttv`: An instance of `libYukiAstronomyTransitVariation` containing the 
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
libYukiAstronomyTransitVariationTTVCorrectLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve,
	variation::libYukiAstronomyTransitVariation,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitVariationTTVCorrectLightCurve(
	lightCurve,
	variation,
	transit.planet.orbit.period
);
function libYukiAstronomyTransitVariationTTVCorrectLightCurve(
    lightCurve::libYukiAstronomyTransitLightCurve,
    variation::libYukiAstronomyTransitVariation,
    planetOrbitPeriod::Real
)
    lightCurveCorrected = deepcopy(lightCurve);
    libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
        lightCurveCorrected,
        variation,
        planetOrbitPeriod,
    );
    return lightCurveCorrected;
end
libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	variation::libYukiAstronomyTransitVariation,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
	lightCurve,
	variation,
	transit.planet.orbit.period
);
function libYukiAstronomyTransitVariationTTVCorrectLightCurve!(
    lightCurve::libYukiAstronomyTransitLightCurve,
    variation::libYukiAstronomyTransitVariation,
    planetOrbitPeriod::Real
)
    originalTime = copy(lightCurve.time);
    for (index, transitMidTime) in enumerate(variation.predictedMidPoints)
        mask = abs.(originalTime .- transitMidTime) .< 
            0.5 * planetOrbitPeriod;
        lightCurve.time[mask] .-= variation.ttvOC[index];
    end
    libYukiAstronomyTransitExtractValidLightCurve!(lightCurve);
    return lightCurve;
end

"""
	libYukiAstronomyTransitVariationTTVOCBrent(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit;
		ocMeasurementLengthDurations::Real = 1.0,
		ocSearchLengthDurations::Real = 0.1,
		minimumPointCount::Integer = 5,
		minimumPointsPerSide::Integer = 2,
		maximumTimingErrorDurations::Real = 0.1,
		profilePointCount::Integer = 101,
		searchCentreFitPointCount::Integer = 5,
		searchWidthMaximumMultiplier::Real = 8.0
	)
	libYukiAstronomyTransitVariationTTVOCBrent(
		lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
		transit::libYukiAstronomyTransit;
		ocMeasurementLengthDurations::Real = 1.0,
		ocSearchLengthDurations::Real = 0.1,
		minimumPointCount::Integer = 5,
		minimumPointsPerSide::Integer = 2,
		maximumTimingErrorDurations::Real = 0.1,
		profilePointCount::Integer = 101,
		searchCentreFitPointCount::Integer = 5,
		searchWidthMaximumMultiplier::Real = 8.0
	)
	libYukiAstronomyTransitVariationTTVOCBrent(
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
		searchWidthMaximumMultiplier::Real = 8.0
	)
Transit timing variation (TTV) analysis using the Brent optimization 
method to find the best-fit transit times for a series of light curves.
The function can be called with either a single light curve and a 
transit object, or with a vector of split light curves and the 
necessary transit parameters.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing the light curve data for a single transit.
- `transit`: An instance of `libYukiAstronomyTransit` containing the 
transit parameters.
- `lightCurvesSplited`: A vector of `libYukiAstronomyTransitLightCurve` 
instances, each representing a split light curve for individual 
transits.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `planetTransitCentreTime`: The predicted mid-transit time based on 
the planet's orbital parameters.
- `planetStellarRadiusRatio`: The ratio of the planet's radius to the 
host star's radius.
- `orbit`: An instance of `Orbits.AbstractOrbit` representing the 
planet's orbit.
- `limbDarkeningFunc`: An instance of `AbstractLimbDark` representing 
the limb darkening function.
- `ocMeasurementLengthDurations`: The length of the time window (in 
units of transit duration) used to measure the O-C values.
- `ocSearchLengthDurations`: The length of the time window (in units 
of transit duration) used to search for the best-fit transit time.
- `minimumPointCount`: The minimum number of data points required 
within the measurement window to perform the TTV analysis.
- `minimumPointsPerSide`: The minimum number of data points required 
on each side of the transit to consider the measurement valid.
- `maximumTimingErrorDurations`: The maximum allowed timing error (in 
units of transit duration) for a valid measurement.
- `profilePointCount`: The number of points used to profile the loss 
function for error estimation.
- `searchCentreFitPointCount`: The number of previous accepted 
measurements used to fit a linear trend for the search centre.
- `searchWidthMaximumMultiplier`: The maximum multiplier for the 
search width to allow for a broader search if the initial search fails.
# Returns
- A `libYukiAstronomyTransitVariation` instance containing the TTV 
data, including the O-C values, their errors, and flags indicating 
whether each measurement is acceptable.
"""
function libYukiAstronomyTransitVariationTTVOCBrent(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit;
	ocMeasurementLengthDurations::Real = 1.0,
	ocSearchLengthDurations::Real = 0.1,
	minimumPointCount::Integer = 5,
	minimumPointsPerSide::Integer = 2,
	maximumTimingErrorDurations::Real = 0.1,
	profilePointCount::Integer = 101,
	searchCentreFitPointCount::Integer = 5,
	searchWidthMaximumMultiplier::Real = 8.0
)
	lightCurvesSplited = libYukiAstronomyTransitSplitLightCurveByPeriod(
		lightCurve,
		transit
	);
	return libYukiAstronomyTransitVariationTTVOCBrent(
		lightCurvesSplited,
		transit;
		ocMeasurementLengthDurations = ocMeasurementLengthDurations,
		ocSearchLengthDurations = ocSearchLengthDurations,
		minimumPointCount = minimumPointCount,
		minimumPointsPerSide = minimumPointsPerSide,
		maximumTimingErrorDurations = maximumTimingErrorDurations,
		profilePointCount = profilePointCount,
		searchCentreFitPointCount = searchCentreFitPointCount,
		searchWidthMaximumMultiplier = searchWidthMaximumMultiplier
	);
end
libYukiAstronomyTransitVariationTTVOCBrent(
	lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
	transit::libYukiAstronomyTransit;
	ocMeasurementLengthDurations::Real = 1.0,
	ocSearchLengthDurations::Real = 0.1,
	minimumPointCount::Integer = 5,
	minimumPointsPerSide::Integer = 2,
	maximumTimingErrorDurations::Real = 0.1,
	profilePointCount::Integer = 101,
	searchCentreFitPointCount::Integer = 5,
	searchWidthMaximumMultiplier::Real = 8.0
) = libYukiAstronomyTransitVariationTTVOCBrent(
	lightCurvesSplited,
	transit.planet.orbit.period,
	transit.planetTransitDuration,
	transit.planetTransitCentreTime,
	transit.planetStellarRadiusRatio,
	SimpleOrbit(
		period = transit.planet.orbit.period, 
		duration = transit.planetTransitDuration
	),
	transit.limbDarkeningFunc;
	ocMeasurementLengthDurations = ocMeasurementLengthDurations,
	ocSearchLengthDurations = ocSearchLengthDurations,
	minimumPointCount = minimumPointCount,
	minimumPointsPerSide = minimumPointsPerSide,
	maximumTimingErrorDurations = maximumTimingErrorDurations,
	profilePointCount = profilePointCount,
	searchCentreFitPointCount = searchCentreFitPointCount,
	searchWidthMaximumMultiplier = searchWidthMaximumMultiplier,
);
function libYukiAstronomyTransitVariationTTVOCBrent(
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
			templateLightCurve = 
				libYukiAstronomyTransitLimbDarkeningTransitFlux(
					t .- dt,
					planetStellarRadiusRatio,
					orbit,
					limbDarkeningFunc,
				);
			templateFlux = templateLightCurve.flux;
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

	return libYukiAstronomyTransitVariation(;
		ttvOC = oc,
		ttvOCErr = ocErr,
		ttvOCIsAcceptable = isAcceptable,
		transitMidPoints = predictedMidPoints .+ oc,
		predictedMidPoints = predictedMidPoints,
		transitDurationRef = planetTransitDuration
	);
end