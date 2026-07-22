using Optim

include("libYukiAstronomyTransitLimbDarkening.jl")

mutable struct libYukiAstronomyTransitTTV
    OC::Vector{Float64}
    OCErr::Vector{Float64}
    isAccepted::Vector{Bool}
    predictedMidPoints::Vector{Float64}
    libYukiAstronomyTransitTTV(
        OC, OCErr, isAccepted, predictedMidPoints) = new(
            OC, OCErr, isAccepted, predictedMidPoints);
end

"""
    libYukiAstronomyTransitTTVOCBrent(
        lightCurvesSplited, 
        planetOrbitPeriod, 
        planetTransitDuration, 
        planetTransitCentreTime, 
        planetStellarRadiusRatio, 
        orbit, 
        limbDarkeningFunc, 
        ocMeasurementLengthDurations = 1.0, 
        ocSearchLengthDurations = 0.1
    )
Calculate the transit timing variations (TTV) using the Brent 
optimization method.
# Arguments
- `lightCurvesSplited`: A vector of `libYukiAstronomyTransitLightCurve`
instances, each representing a light curve for a single transit.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `planetTransitCentreTime`: The expected center time of the first 
transit.
- `planetStellarRadiusRatio`: The ratio of the planet's radius to 
the star's radius.
- `orbit`: An instance of a subtype of `AbstractOrbit` representing 
the orbital parameters of the planet.
- `limbDarkeningFunc`: An instance of a subtype of `AbstractLimbDark` 
representing the limb darkening model to be used for the transit flux 
calculation.
- `ocMeasurementLengthDurations`: The length of time (in units of 
transit duration) to consider for measuring the transit timing 
variations. Default is 1.0.
- `ocSearchLengthDurations`: The length of time (in units of transit 
duration) to search for the optimal transit timing variations. 
Default is 0.1.
# Returns
- A vector of calculated transit timing variations (TTV) for each 
transit in the input light curves.
"""
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
)
	nTransits = length(lightCurvesSplited);
	oc = fill(NaN, nTransits);
	ocErr = fill(NaN, nTransits);
	isAccepted = falses(nTransits);

	predictedMidPoints = planetTransitCentreTime .+ 
        (0 : nTransits - 1) .* planetOrbitPeriod;

	searchWidth = ocSearchLengthDurations * planetTransitDuration;

	for (index, t0) in enumerate(predictedMidPoints)
		lightCurve = lightCurvesSplited[index];

		mask = abs.(lightCurve.time .- t0) .<
			ocMeasurementLengthDurations * planetTransitDuration;

		count(mask) >= minimumPointCount || continue

		t = lightCurve.time[mask] .- t0;
		f = lightCurve.flux[mask];

		all(isfinite, t) || continue
		all(isfinite, f) || continue

		fluxErrProvided =
			lightCurve.fluxErr !== nothing &&
			!isempty(lightCurve.fluxErr);
		if fluxErrProvided
			length(lightCurve.fluxErr) ==
				length(lightCurve.time) || continue
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

			all(isfinite, templateFlux) || return Inf
			designMatrix = hcat(
				ones(length(t)),
                templateFlux .- 1.0,
            );

			weightedDesignMatrix = designMatrix .* 
                reshape(sqrtWeight, :, 1);

			parameters = weightedDesignMatrix \
				(f .* sqrtWeight);

			residual = (f .- designMatrix * parameters) .*
				sqrtWeight;

			return sum(abs2, residual);
		end

		result = optimize(
            loss,
            -searchWidth,
            searchWidth,
            Brent(),
        );
		Optim.converged(result) || continue

		dtBest = Optim.minimizer(result);
		lossBest = Optim.minimum(result);

		isfinite(dtBest) || continue
		isfinite(lossBest) || continue
		abs(dtBest) < 0.9 * searchWidth || continue
		shiftedTime = t .- dtBest
		count(shiftedTime .< 0) >= minimumPointsPerSide ||
			continue
		count(shiftedTime .> 0) >= minimumPointsPerSide ||
			continue

		degreesOfFreedom = length(t) - 3;
		degreesOfFreedom > 0 || continue
		profileLossIncrease = fluxErrProvided ?
			1.0 : lossBest / degreesOfFreedom;

		isfinite(profileLossIncrease) || continue
		profileLossIncrease > 0 || continue

		dtGrid = collect(
            range(
                -searchWidth,
                searchWidth;
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

		isempty(leftIndices) && continue
		isempty(rightIndices) && continue
		dtLeft = dtGrid[last(leftIndices)];
		dtRight = dtGrid[first(rightIndices)];

		timingError = (dtRight - dtLeft) / 2;
		isfinite(timingError) || continue
		timingError > 0 || continue
		timingError <= maximumTimingErrorDurations * 
            planetTransitDuration || continue;

		oc[index] = dtBest;
		ocErr[index] = timingError;
		isAccepted[index] = true;
	end

    return libYukiAstronomyTransitTTV(
        oc,
        ocErr,
        isAccepted,
        predictedMidPoints,
    );
end
