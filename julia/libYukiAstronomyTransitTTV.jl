using Optim

include("libYukiAstronomyTransitLimbDarkening.jl")

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
Calculate the transit timing variations (TTV) using the Brent optimization method.
# Arguments
- `lightCurvesSplited`: A vector of `libYukiAstronomyTransitLightCurve` instances, each representing a light curve for a single transit.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `planetTransitCentreTime`: The expected center time of the first transit.
- `planetStellarRadiusRatio`: The ratio of the planet's radius to the star's radius.
- `orbit`: An instance of a subtype of `AbstractOrbit` representing the orbital parameters of the planet.
- `limbDarkeningFunc`: An instance of a subtype of `AbstractLimbDark` representing the limb darkening model to be used for the transit flux calculation.
- `ocMeasurementLengthDurations`: The length of time (in units of transit duration) to consider for measuring the transit timing variations. Default is 1.0.
- `ocSearchLengthDurations`: The length of time (in units of transit duration) to search for the optimal transit timing variations. Default is 0.1.
# Returns
- A vector of calculated transit timing variations (TTV) for each transit in the input light curves.
"""
function libYukiAstronomyTransitTTVOCBrent(
    lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve},
    planetOrbitPeriod::Real,
    planetTransitDuration::Real,
    planetTransitCentreTime::Real,
    planetStellarRadiusRatio::Real,
    orbit::Orbits.AbstractOrbit,
    limbDarkeningFunc::AbstractLimbDark,
    ocMeasurementLengthDurations::Real = 1.0,
    ocSearchLengthDurations::Real = 0.1,
)
    ocSerialBrentDay =
        fill(NaN, length(lightCurvesSplited))

    predictedMidPoints = planetTransitCentreTime .+
        (0:length(lightCurvesSplited)-1) .* planetOrbitPeriod

    for (index, t0) in enumerate(predictedMidPoints)
        lightCurve = lightCurvesSplited[index]

        mask = abs.(lightCurve.time .- t0) .<
            ocMeasurementLengthDurations * planetTransitDuration;

        count(mask) >= 5 || continue

        t = lightCurve.time[mask] .- t0;
        f = lightCurve.flux[mask];
        fErr = lightCurve.fluxErr[mask];

        all(isfinite, t) || continue
        all(isfinite, f) || continue
        all(isfinite, fErr) || continue
        all(fErr .> 0) || continue

        weight = 1.0 ./ fErr.^2;
        sqrtWeight = sqrt.(weight);

        function optLoss(dt)
            templateFlux = libYukiAstronomyTransitLimbDarkeningTransitFlux(
                t .- dt,
                planetStellarRadiusRatio,
                orbit,
                limbDarkeningFunc,
            );
            transitShape = templateFlux .- 1.0;
            designMatrix = hcat(
                ones(length(t)),
                transitShape,
            );
            parameters = (designMatrix .* sqrtWeight) \
                (f .* sqrtWeight);
            fittedFlux = designMatrix * parameters;
            return sum(weight .* (f .- fittedFlux).^2);
        end
        searchWidth = ocSearchLengthDurations * planetTransitDuration;
        result = optimize(
            optLoss,
            -searchWidth,
            searchWidth,
            Brent(),
        );

        ocSerialBrentDay[index] = Optim.minimizer(result);
    end
    return ocSerialBrentDay;
end