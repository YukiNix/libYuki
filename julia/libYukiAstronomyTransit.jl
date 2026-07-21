using JLD2, Measurements, Dates, Tables, DataFrames, ProgressMeter;
using Turing, Transits, Random, Orbits, StableRNGs;
using PyCall;
using FITSIO;

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysics.jl")

"""
	libYukiAstronomyTransitLightCurve(
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
- An instance of `libYukiAstronomyTransitLightCurve`.
"""
mutable struct libYukiAstronomyTransitLightCurve
	time::AbstractVector
	flux::AbstractVector
	fluxErr::AbstractVector
	libYukiAstronomyTransitLightCurve(time, flux, fluxErr) = new(time, flux, fluxErr);
end

include("libYukiAstronomyTransitLimbDarkening.jl")

function libYukiAstronomyTransitSplitLightCurveByPeriod(lightCurve::libYukiAstronomyTransitLightCurve, transitCentreTime::Real, planetOrbitPeriod::Real)
	nPeriods = ceil(Int, (maximum(lightCurve.time) - minimum(lightCurve.time)) / planetOrbitPeriod);
	splitLightCurves = Vector{libYukiAstronomyTransitLightCurve}(undef, nPeriods);
	firstPeriodTime = transitCentreTime - planetOrbitPeriod / 2;
	for index in 1 : nPeriods
		startTime = firstPeriodTime + (index - 1) * planetOrbitPeriod;
		endTime = startTime + planetOrbitPeriod;
		mask = (lightCurve.time .>= startTime) .& (lightCurve.time .< endTime);
		if any(mask)
			splitLightCurves[index] = libYukiAstronomyTransitLightCurve(
				lightCurve.time[mask],
				lightCurve.flux[mask],
				lightCurve.fluxErr[mask]
			);
		else
			splitLightCurves[index] = libYukiAstronomyTransitLightCurve(
				Float64[],
				Float64[],
				Float64[]
			);
		end
	end
	return splitLightCurves;
end

"""
	libYukiAstronomyTransitLoadLightCurveFromFITS(
		fileName
	)
Load a light curve from a FITS file.
# Arguments
- `fileName`: The path to the FITS file containing the light curve data.
# Returns
- An instance of `libYukiAstronomyTransitLightCurve` containing the time, flux, and flux error data.
"""
function libYukiAstronomyTransitLoadLightCurveFromFITS(fileName::String)
	time, flux, fluxErr = FITS(fileName, "r") do file
		table = file[2];
		return (
			read(table, "TIME"),
			read(table, "FLUX"),
			read(table, "FLUX_ERR")
		);
	end
	return libYukiAstronomyTransitLightCurve(
		time,
		flux,
		fluxErr
	);
end

"""
	libYukiAstronomyTransitSaveLightCurveToFITS(
		lightCurve,
		filePath
	)
Save a light curve to a FITS file.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `filePath`: The path to the FITS file where the light curve will be saved.
"""
function libYukiAstronomyTransitSaveLightCurveToFITS(lightCurve::libYukiAstronomyTransitLightCurve, filePath::String)
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
	libYukiAstronomyTransitNormalizeLightCurve(
		lightCurve
	)
	libYukiAstronomyTransitNormalizeLightCurve!(
		lightCurve
	)
Normalize the flux of a light curve by dividing it by its median value. The flux error is also adjusted accordingly. 
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` with normalized flux and flux error.
# Notes
- The function `libYukiAstronomyTransitNormalizeLightCurve!` modifies the input light curve in place, while `libYukiAstronomyTransitNormalizeLightCurve` returns a new normalized light curve.
"""
function libYukiAstronomyTransitNormalizeLightCurve(lightCurve::libYukiAstronomyTransitLightCurve)
	lightCurveNormalized = deepcopy(lightCurve);
	libYukiAstronomyTransitNormalizeLightCurve!(lightCurveNormalized);
	return lightCurveNormalized;
end
function libYukiAstronomyTransitNormalizeLightCurve!(lightCurve::libYukiAstronomyTransitLightCurve)
	medianFlux = median(lightCurve.flux);
	if !isfinite(medianFlux) || medianFlux == 0.0
		error("Normalization failed: median flux is not finite or zero.")
	end
	lightCurve.flux ./= medianFlux;
	lightCurve.fluxErr ./= abs(medianFlux);
end

"""
	libYukiAstronomyTransitGetValidLightCurve(
		lightCurve
	)
	libYukiAstronomyTransitGetValidLightCurve!(
		lightCurve
	)
Filter out invalid data points from a light curve, ensuring that only finite and positive flux error values are retained. The filtered light curve is sorted by time.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` containing only valid data points.
- Notes
- The function `libYukiAstronomyTransitGetValidLightCurve!` modifies the input light curve in place, while `libYukiAstronomyTransitGetValidLightCurve` returns a new filtered light curve.
"""
function libYukiAstronomyTransitGetValidLightCurve!(lightCurve::libYukiAstronomyTransitLightCurve)
	validMask = isfinite.(lightCurve.time) .& isfinite.(lightCurve.flux) .& isfinite.(lightCurve.fluxErr) .& (lightCurve.fluxErr .> 0);
	lightCurve.time = lightCurve.time[validMask];
	lightCurve.flux = lightCurve.flux[validMask];
	lightCurve.fluxErr = lightCurve.fluxErr[validMask];

	sortIdx = sortperm(lightCurve.time);
	lightCurve.time = map(x -> Float64(x), lightCurve.time[sortIdx]);
	lightCurve.flux = map(x -> Float64(x), lightCurve.flux[sortIdx]);
	lightCurve.fluxErr = map(x -> Float64(x), lightCurve.fluxErr[sortIdx]);
end
function libYukiAstronomyTransitGetValidLightCurve(lightCurve::libYukiAstronomyTransitLightCurve)
	lightCurveValid = deepcopy(lightCurve);
	libYukiAstronomyTransitGetValidLightCurve!(lightCurveValid);
	return lightCurveValid;
end

"""
	libYukiAstronomyTransitBinLightCurve(
		lightCurve, 
		binBoxNumber
	)
Bin a light curve into a specified number of bins, averaging the flux and flux error within each bin. The binned light curve is returned with time values corresponding to the centre of each bin.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `binBoxNumber`: The number of bins to divide the light curve into.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` containing the binned time, flux, and flux error values.
"""
function libYukiAstronomyTransitBinLightCurve(lightCurve::libYukiAstronomyTransitLightCurve, binBoxNumber::Integer) 
	binEdge = range(minimum(lightCurve.time), maximum(lightCurve.time); length = binBoxNumber + 1)
    binnedTime = @. 0.5 * (binEdge[1 : end - 1] + binEdge[2 : end])
    binnedFlux = zeros(length(binnedTime))
	binnedFluxErr = zeros(length(binnedTime))
   for index in eachindex(binnedTime)
        mask = (lightCurve.time .>= binEdge[index]) .& (lightCurve.time .< binEdge[index + 1]);
        if any(mask)
            binnedFlux[index] = mean(lightCurve.flux[mask]);
			binnedFluxErr[index] = sqrt(sum(lightCurve.fluxErr[mask] .^ 2)) / sum(mask);
        else
            binnedFlux[index] = NaN;
			binnedFluxErr[index] = NaN;
        end
    end
    nNaNIndices = findall(x -> !isnan(x), binnedFlux);

    return libYukiAstronomyTransitLightCurve(binnedTime[nNaNIndices], binnedFlux[nNaNIndices], binnedFluxErr[nNaNIndices]);
end

"""
	libYukiAstronomyTransitFoldLightCurve(
		lightCurve, 
		planetOrbitPeriod, 
		transitCentreTime
	)
	libYukiAstronomyTransitFoldLightCurve!(
		lightCurve, 
		planetOrbitPeriod, 
		transitCentreTime
	)
Fold a light curve based on the planet's orbital period and the transit centre time. The time values are adjusted to be within the range of -P/2 to P/2, where P is the planet's orbital period. The light curve is then sorted by the adjusted time values.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `planetOrbitPeriod`: The orbital period of the planet.
- `transitCentreTime`: The time of the transit centre.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` with folded time values.
- Notes
- The function `libYukiAstronomyTransitFoldLightCurve!` modifies the input light curve in place, while `libYukiAstronomyTransitFoldLightCurve` returns a new folded light curve.
"""
function libYukiAstronomyTransitFoldLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve,
	planetOrbitPeriod::Real, 
	transitCentreTime::Real
)
	lightCurveFolded = deepcopy(lightCurve);
	libYukiAstronomyTransitFoldLightCurve!(lightCurveFolded, planetOrbitPeriod, transitCentreTime);
	return lightCurveFolded;
end
function libYukiAstronomyTransitFoldLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve,
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
# function libYukiAstronomyTransitDetrendByWotan(times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}, stellarRadius::Real, stellarMass::Real, planetOrbitPeriod::Real) 
# 	pyWotan = pyimport("wotan");
# 	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
# 		window = 0.75;
# 	else
# 		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
# 	end
# 	sortedTimes, sortedFluxes = libYukiBasicSortElementsByOrderVector(times, fluxes);
# 	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
# 	return sortedTimes, AbstractVector{<:Real}(detrendedFluxes);
# end
# function libYukiAstronomyTransitDetrendByWotan(times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}, stellarRadius::Real, stellarMass::Real, planetOrbitPeriod::Real) 
# 	pyWotan = pyimport("wotan");
# 	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
# 		window = 0.75;
# 	else
# 		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
# 	end
# 	sortedTimes, sortedFluxes = libYukiBasicSortElementsByOrderVector(times, fluxes);
# 	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
# 	return sortedTimes, detrendedFluxes .± (Measurements.uncertainty.(sortedFluxes) ./ Measurements.value.(sortedFluxes));
# end

"""
	libYukiAstronomyTransitTransitDurationEvaluate(
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
function libYukiAstronomyTransitTransitDurationEvaluate(
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
