using JLD2, Measurements, Dates, Tables, DataFrames, ProgressMeter;
using Turing, Transits, Random, Orbits, StableRNGs;
using PyCall;
using FITSIO;

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysics.jl")

"""
	libYukiTransitLightCurve(
		time, 
		flux, 
		fluxErr
	)
Define a light curve for transit analysis with time, flux, and flux error.
# Arguments
- `time`: A vector of time values.
- `flux`: A vector of flux values corresponding to the time values.
- `fluxErr`: A vector of flux error values corresponding to the time values.
# Returns
- An instance of `libYukiTransitLightCurve`.
"""
mutable struct libYukiTransitLightCurve
	time::AbstractVector
	flux::AbstractVector
	fluxErr::AbstractVector
	libYukiTransitLightCurve(time, flux, fluxErr) = new(time, flux, fluxErr);
end

# Quadratic limb-darkening transit flux model.
# Dependency: Turing, Transits, Orbits.
# TODO: Validate & Example.
@model function libYukiTransitQuadraticLimbDarkeningTransitFluxModel(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, f)
	u ~ Kipping13()
	σ ~ Truncated(Cauchy(0, 0.01), 0, Inf)
	fModel = libYukiTransitQuadraticLimbDarkeningTransitFlux(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, collect(u))
	for i in eachindex(t)
		f[i] ~ Normal(fModel[i], σ)
	end
end

"""
	libYukiTransitLoadLightCurveFromFITS(
		fileName
	)
Load a light curve from a FITS file.
# Arguments
- `fileName`: The path to the FITS file containing the light curve data.
# Returns
- An instance of `libYukiTransitLightCurve` containing the time, flux, and flux error data.
"""
function libYukiTransitLoadLightCurveFromFITS(fileName::String)
	time, flux, fluxErr = FITS(fileName, "r") do file
		table = file[2];
		return (
			read(table, "TIME"),
			read(table, "FLUX"),
			read(table, "FLUX_ERR")
		);
	end
	return libYukiTransitLightCurve(
		time,
		flux,
		fluxErr
	);
end

"""
	libYukiTransitSaveLightCurveToFITS(
		lightCurve,
		filePath
	)
Save a light curve to a FITS file.
# Arguments
- `lightCurve`: An instance of `libYukiTransitLightCurve`.
- `filePath`: The path to the FITS file where the light curve will be saved.
"""
function libYukiTransitSaveLightCurveToFITS(lightCurve::libYukiTransitLightCurve, filePath::String)
	length(lightCurve.time) == length(lightCurve.flux) == length(lightCurve.fluxErr) ||
        throw(DimensionMismatch("The three vectors must have equal lengths."))
    FITS(filePath, "w") do file
        write(file, Dict(
            "TIME"     => Float64.(lightCurve.time),
            "FLUX"     => Float64.(lightCurve.flux),
            "FLUX_ERR" => Float64.(lightCurve.fluxErr),
        ));
    end
end

"""
	libYukiTransitNormalizeLightCurve(
		lightCurve
	)
	libYukiTransitNormalizeLightCurve!(
		lightCurve
	)
Normalize the flux of a light curve by dividing it by its median value. The flux error is also adjusted accordingly. 
# Arguments
- `lightCurve`: An instance of `libYukiTransitLightCurve`.
# Returns
- A new instance of `libYukiTransitLightCurve` with normalized flux and flux error.
# Notes
- The function `libYukiTransitNormalizeLightCurve!` modifies the input light curve in place, while `libYukiTransitNormalizeLightCurve` returns a new normalized light curve.
"""
function libYukiTransitNormalizeLightCurve(lightCurve::libYukiTransitLightCurve)
	lightCurveNormalized = deepcopy(lightCurve);
	libYukiTransitNormalizeLightCurve!(lightCurveNormalized);
	return lightCurveNormalized;
end
function libYukiTransitNormalizeLightCurve!(lightCurve::libYukiTransitLightCurve)
	medianFlux = median(lightCurve.flux);
	if !isfinite(medianFlux) || medianFlux == 0.0
		error("Normalization failed: median flux is not finite or zero.")
	end
	lightCurve.flux ./= medianFlux;
	lightCurve.fluxErr ./= abs(medianFlux);
end

"""
	libYukiTransitGetValidLightCurve(
		lightCurve
	)
	libYukiTransitGetValidLightCurve!(
		lightCurve
	)
Filter out invalid data points from a light curve, ensuring that only finite and positive flux error values are retained. The filtered light curve is sorted by time.
# Arguments
- `lightCurve`: An instance of `libYukiTransitLightCurve`.
# Returns
- A new instance of `libYukiTransitLightCurve` containing only valid data points.
- Notes
- The function `libYukiTransitGetValidLightCurve!` modifies the input light curve in place, while `libYukiTransitGetValidLightCurve` returns a new filtered light curve.
"""
function libYukiTransitGetValidLightCurve!(lightCurve::libYukiTransitLightCurve)
	validMask = isfinite.(lightCurve.time) .& isfinite.(lightCurve.flux) .& isfinite.(lightCurve.fluxErr) .& (lightCurve.fluxErr .> 0);
	lightCurve.time = lightCurve.time[validMask];
	lightCurve.flux = lightCurve.flux[validMask];
	lightCurve.fluxErr = lightCurve.fluxErr[validMask];

	sortIdx = sortperm(lightCurve.time);
	lightCurve.time = map(x -> Float64(x), lightCurve.time[sortIdx]);
	lightCurve.flux = map(x -> Float64(x), lightCurve.flux[sortIdx]);
	lightCurve.fluxErr = map(x -> Float64(x), lightCurve.fluxErr[sortIdx]);
end
function libYukiTransitGetValidLightCurve(lightCurve::libYukiTransitLightCurve)
	lightCurveValid = deepcopy(lightCurve);
	libYukiTransitGetValidLightCurve!(lightCurveValid);
	return lightCurveValid;
end

# Transit flux with limb-darkening.
# Dependency: Transits, Orbits.
# TODO: Validate & Example.
function libYukiTransitQuadraticLimbDarkeningTransitFlux(times::AbstractVector{<:Real}, orbitPeriod::Real, transitDuration::Real, planetStellarRadiusRatio::Real, limbQuadraticDarkeningParameters::AbstractVector{<:Real}) 
	orbit = SimpleOrbit(period = orbitPeriod, duration = transitDuration);
    ld = QuadLimbDark(limbQuadraticDarkeningParameters);
	return ld.(Ref(orbit), times, planetStellarRadiusRatio)
end

# Transit flux with polynomial limb-darkening.
# Dependency: Transits, Orbits.
# TODO: Validate & Example.
function libYukiTransitPolynomialLimbDarkeningTransitFlux(times::AbstractVector{<:Real}, orbitPeriod::Real, transitDuration::Real, planetStellarRadiusRatio::Real, limbDarkeningParameters::AbstractVector{<:Real}) 
	orbit = SimpleOrbit(period = orbitPeriod, duration = transitDuration);
	ld = PolynomialLimbDark(limbDarkeningParameters);
	return ld.(Ref(orbit), times, planetStellarRadiusRatio)
end

# Polynomial limb-darkening transit flux model.
# Dependency: Turing, Transits, Orbits.
# TODO: Validate & Example.
@model function libYukiTransitPolynomialLimbDarkeningTransitFluxModel(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, f, nLimbDarkeningParameters)
	u ~ filldist(Uniform(-1, 1), nLimbDarkeningParameters)
	σ ~ Truncated(Cauchy(0, 0.01), 0, Inf)
	fModel = libYukiTransitPolynomialLimbDarkeningTransitFlux(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, collect(u))
	for i in eachindex(t)
		f[i] ~ Normal(fModel[i], σ)
	end
end

# Bin lightcurve.
# TODO: Validate & Example.
function libYukiTransitBinLightCurve(times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}, binBoxNumber::Integer) 
	binEdges = range(minimum(times), maximum(times); length = binBoxNumber + 1)
    binnedTimes = @. 0.5 * (binEdges[1 : end - 1] + binEdges[2 : end])
    binnedFluxes = zeros(T, length(binnedTimes))
    @showprogress for index in eachindex(binnedTimes)
        mask = (times .>= binEdges[index]) .& (times .< binEdges[index + 1]);
        if any(mask)
            binnedFluxes[index] = mean(fluxes[mask]);
        else
            binnedFluxes[index] = NaN;
        end
    end
    nNaNIndices = findall(x -> !isnan(x), binnedFluxes);
    return binnedTimes[nNaNIndices], binnedFluxes[nNaNIndices]
end


"""
	libYukiTransitFoldLightCurve(
		lightCurve, 
		planetOrbitPeriod, 
		transitCentreTime
	)
	libYukiTransitFoldLightCurve!(
		lightCurve, 
		planetOrbitPeriod, 
		transitCentreTime
	)
Fold a light curve based on the planet's orbital period and the transit center time. The time values are adjusted to be within the range of -P/2 to P/2, where P is the planet's orbital period. The light curve is then sorted by the adjusted time values.
# Arguments
- `lightCurve`: An instance of `libYukiTransitLightCurve`.
- `planetOrbitPeriod`: The orbital period of the planet.
- `transitCentreTime`: The time of the transit center.
# Returns
- A new instance of `libYukiTransitLightCurve` with folded time values.
- Notes
- The function `libYukiTransitFoldLightCurve!` modifies the input light curve in place, while `libYukiTransitFoldLightCurve` returns a new folded light curve.
"""
function libYukiTransitFoldLightCurve(
	lightCurve::libYukiTransitLightCurve,
	planetOrbitPeriod::Real, 
	transitCentreTime::Real
)
	lightCurveFolded = deepcopy(lightCurve);
	libYukiTransitFoldLightCurve!(lightCurveFolded, planetOrbitPeriod, transitCentreTime);
	return lightCurveFolded;
end
function libYukiTransitFoldLightCurve!(
	lightCurve::libYukiTransitLightCurve,
	planetOrbitPeriod::Real, 
	transitCentreTime::Real
)
	lightCurve.time = ((lightCurve.time .- transitCentreTime) .% planetOrbitPeriod);
	lightCurve.time[lightCurve.time .> planetOrbitPeriod / 2] .-= planetOrbitPeriod;

	sortedIndex = sortperm(lightCurve.time);
	lightCurve.time = lightCurve.time[sortedIndex];
	lightCurve.flux = lightCurve.flux[sortedIndex];
	lightCurve.fluxErr = lightCurve.fluxErr[sortedIndex];
end

# Detrending lightcurve with Wotan.
# Dependency: PyCall, wotan(python package).
# TODO: Validate & Example.
function libYukiTransitDetrendByWotan(times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}, stellarRadius::Real, stellarMass::Real, planetOrbitPeriod::Real) 
	pyWotan = pyimport("wotan");
	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
		window = 0.75;
	else
		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
	end
	sortedTimes, sortedFluxes = libYukiBasicSortElementsByOrderVector(times, fluxes);
	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
	return sortedTimes, AbstractVector{<:Real}(detrendedFluxes);
end
function libYukiTransitDetrendByWotan(times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}, stellarRadius::Real, stellarMass::Real, planetOrbitPeriod::Real) 
	pyWotan = pyimport("wotan");
	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
		window = 0.75;
	else
		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
	end
	sortedTimes, sortedFluxes = libYukiBasicSortElementsByOrderVector(times, fluxes);
	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
	return sortedTimes, detrendedFluxes .± (Measurements.uncertainty.(sortedFluxes) ./ Measurements.value.(sortedFluxes));
end

"""
	libYukiTransitTransitDurationEvaluate(
		planetOrbitPeriod, 
		stellarRadius, 
		planetRadius, 
		planetOrbitSemiMajorAxis, 
		planetOrbitInclination
	)
Evaluate the transit duration based on the planet's orbital parameters and stellar properties.
# Arguments
- `planetOrbitPeriod`: The orbital period of the planet.
- `stellarRadius`: The radius of the host star.
- `planetRadius`: The radius of the planet.
- `planetOrbitSemiMajorAxis`: The semi-major axis of the planet's orbit.
- `planetOrbitInclination`: The inclination of the planet's orbit in degrees.
# Returns
- The transit duration.
"""
function libYukiTransitTransitDurationEvaluate(
    planetOrbitPeriod::Real,
    stellarRadius::Real,
    planetRadius::Real,
    planetOrbitSemiMajorAxis::Real,
    planetOrbitInclination::Real,
)
    radiusRatio = planetRadius / stellarRadius;
    scaledSemiMajorAxis = planetOrbitSemiMajorAxis / stellarRadius;
    inclinationSin = sind(planetOrbitInclination);
    impactParameter = scaledSemiMajorAxis * cosd(planetOrbitInclination);
    contactRadius = 1 + radiusRatio;
    impactParameter < contactRadius || return zero(float(planetOrbitPeriod))

    argument = sqrt(contactRadius ^ 2 - impactParameter ^ 2) / (scaledSemiMajorAxis * inclinationSin);

    return (planetOrbitPeriod / π * asin(clamp(argument, -1, 1)));
end
