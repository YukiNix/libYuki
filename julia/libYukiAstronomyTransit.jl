using JLD2, Measurements, Dates, Tables, DataFrames, ProgressMeter;
using Turing, Transits, Random, Orbits, StableRNGs;
using PyCall;
using FITSIO;

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysics.jl")

"""
	libYukiAstronomyTransit(;
		host::libYukiAstronomyBody = libYukiAstronomyBody(), 
		planet::libYukiAstronomyBody = libYukiAstronomyBody(), 
		planetTransitCentreTime::Real = NaN, 
		planetTransitDuration::Real = NaN, 
		planetStellarRadiusRatio::Real = NaN,
		limbDarkeningFunc::AbstractLimbDark = QuadLimbDark([0.4804, 0.1867])
	)
Create a new `libYukiAstronomyTransit` instance with the specified 
parameters.
# Arguments
- `host`: An instance of `libYukiAstronomyBody` representing the 
host star.
- `planet`: An instance of `libYukiAstronomyBody` representing 
the transiting planet.
- `planetTransitCentreTime`: The time of the transit centre.
- `planetTransitDuration`: The duration of the transit.
- `planetStellarRadiusRatio`: The ratio of the planet's radius 
to the host star's radius.
- `limbDarkeningFunc`: An instance of `AbstractLimbDark` representing the limb darkening function. Default is `QuadLimbDark([0.4804, 0.1867])`, from (Hippke, 2019).
# Returns
- An instance of `libYukiAstronomyTransit`.
"""
mutable struct libYukiAstronomyTransit
	host::libYukiAstronomyBody
	planet::libYukiAstronomyBody
	planetTransitCentreTime::Real
	planetTransitDuration::Real
	planetStellarRadiusRatio::Real
	systemDistance::Real
	limbDarkeningFunc::AbstractLimbDark
	libYukiAstronomyTransit(;
		host::libYukiAstronomyBody = libYukiAstronomyBody(),
		planet::libYukiAstronomyBody = libYukiAstronomyBody(),
		planetTransitCentreTime::Real = NaN,
		planetTransitDuration::Real = NaN,
		planetStellarRadiusRatio::Real = NaN,
		systemDistance::Real = NaN,
		limbDarkeningFunc::AbstractLimbDark = QuadLimbDark(
			[0.4804, 0.1867]
		)
	) = new(
		host, 
		planet, 
		planetTransitCentreTime, 
		planetTransitDuration, 
		planetStellarRadiusRatio,
		systemDistance,
		limbDarkeningFunc
	);
end

"""
	libYukiAstronomyTransitLightCurve(;
		time::AbstractVector{<:Real} = Float64[], 
		flux::AbstractVector{<:Real} = Float64[], 
		fluxErr::AbstractVector{<:Real} = Float64[]
	)
Define a light curve for transit analysis with time, flux, and 
flux error.
# Arguments
- `time`: A vector of time values.
- `flux`: A vector of flux values corresponding to the time values.
- `fluxErr`: A vector of flux error values corresponding to the 
time values.
# Returns
- An instance of `libYukiAstronomyTransitLightCurve`.
"""
mutable struct libYukiAstronomyTransitLightCurve
	time::AbstractVector{<:Real}
	flux::AbstractVector{<:Real}
	fluxErr::AbstractVector{<:Real}
	libYukiAstronomyTransitLightCurve(;
		time::AbstractVector{<:Real} = Float64[], 
		flux::AbstractVector{<:Real} = Float64[], 
		fluxErr::AbstractVector{<:Real} = Float64[]
	) = new(
		time, 
		flux, 
		fluxErr
	);
end

include("libYukiAstronomyTransitLimbDarkening.jl")
include("libYukiAstronomyTransitVariation.jl")

"""
	libYukiAstronomyTransitLoadTransit(
		filePath::String
	)
Load a `libYukiAstronomyTransit` instance from a file using JLD2 format.
# Arguments
- `filePath`: The path to the file from which the transit data will 
be loaded.
# Returns
- A `libYukiAstronomyTransit` instance containing the loaded 
transit data.
"""
function libYukiAstronomyTransitLoadTransit(
	filePath::String
)
	@load filePath transit;
	return transit;
end

"""
	libYukiAstronomyTransitSaveTransit(
		transit::libYukiAstronomyTransit,
		filePath::String
	)
Save a `libYukiAstronomyTransit` instance to a file using JLD2 format.
# Arguments
- `transit`: An instance of `libYukiAstronomyTransit` to be saved.
- `filePath`: The path to the file where the transit data will be saved.
"""
function libYukiAstronomyTransitSaveTransit(
	transit::libYukiAstronomyTransit,
	filePath::String
)
	@save filePath transit;
end

"""
	libYukiAstronomyTransitExtractTransitPart(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitExtractTransitPart(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		planetOrbitPeriod::Real,
		planetTransitCentreTime::Real,
		planetTransitDuration::Real
	)
	libYukiAstronomyTransitExtractTransitPart(
		lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve}, 
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitExtractTransitPart(
		lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve}, 
		planetOrbitPeriod::Real,
		planetTransitCentreTime::Real,
		planetTransitDuration::Real
	)
	libYukiAstronomyTransitExtractTransitPart!(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitExtractTransitPart!(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		planetOrbitPeriod::Real,
		planetTransitCentreTime::Real,
		planetTransitDuration::Real
	)
Extract the transit part of a light curve based on the planet's 
orbital period, transit centre time, and transit duration. The 
function can be called with either a single 
`libYukiAstronomyTransitLightCurve` instance or a vector of such 
instances. In the latter case, the function will extract the transit 
part from each light curve segment. The function with the `!` 
suffix modifies the input light curve in place.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing time values and flux values.
- `transit`: An instance of `libYukiAstronomyTransit` containing 
the planet's orbital period, transit centre time, and transit duration.
- `lightCurvesSplited`: A vector of 
`libYukiAstronomyTransitLightCurve` instances, each representing 
a segment of the light curve.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitCentreTime`: The time of the transit centre.
- `planetTransitDuration`: The duration of the transit.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` containing 
the extracted transit part of the light curve.
"""
libYukiAstronomyTransitExtractTransitPart(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitExtractTransitPart(
	lightCurve,
	transit.planet.orbit.period,
	transit.planetTransitCentreTime,
	transit.planetTransitDuration
);
function libYukiAstronomyTransitExtractTransitPart(
	lightCurve::libYukiAstronomyTransitLightCurve,
	planetOrbitPeriod::Real,
	planetTransitCentreTime::Real,
	planetTransitDuration::Real
)
	lightCurveExtracted = deepcopy(lightCurve);
	libYukiAstronomyTransitExtractTransitPart!(
		lightCurveExtracted,
		planetOrbitPeriod,
		planetTransitCentreTime,
		planetTransitDuration
	);
	return lightCurveExtracted;
end
libYukiAstronomyTransitExtractTransitPart(
	lightCurvesSplited::Vector{libYukiAstronomyTransitLightCurve}, 
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitExtractTransitPart(
	lightCurvesSplited,
	transit.planet.orbit.period,
	transit.planetTransitCentreTime,
	transit.planetTransitDuration
);
function libYukiAstronomyTransitExtractTransitPart(
	lightCurveSplited::Vector{libYukiAstronomyTransitLightCurve}, 
	planetOrbitPeriod::Real,
	planetTransitCentreTime::Real,
	planetTransitDuration::Real
)
	lightCurveTransit = libYukiAstronomyTransitLightCurve();
	for (index, lightCurve) in enumerate(lightCurveSplited)
		startTime = planetTransitCentreTime + 
			(index - 1) * planetOrbitPeriod - 
			planetTransitDuration / 2;
		endTime = startTime + planetTransitDuration;
		mask = (lightCurve.time .>= startTime) .& 
			(lightCurve.time .< endTime);
		if any(mask)
			append!(
				lightCurveTransit.time, 
				lightCurve.time[mask]
			);
			append!(
				lightCurveTransit.flux, 
				lightCurve.flux[mask]
			);
			append!(
				lightCurveTransit.fluxErr, 
				lightCurve.fluxErr[mask]
			);
		end
	end
	return lightCurveTransit;
end
libYukiAstronomyTransitExtractTransitPart!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitExtractTransitPart!(
	lightCurve,
	transit.planet.orbit.period,
	transit.planetTransitCentreTime,
	transit.planetTransitDuration
);
function libYukiAstronomyTransitExtractTransitPart!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	planetOrbitPeriod::Real,
	planetTransitCentreTime::Real,
	planetTransitDuration::Real
)
	lightCurveExtracted = libYukiAstronomyTransitExtractTransitPart(
		libYukiAstronomyTransitSplitLightCurveByPeriod(
			lightCurve, 
			planetOrbitPeriod, 
			planetTransitCentreTime
		),
		planetOrbitPeriod,
		planetTransitCentreTime,
		planetTransitDuration
	);
	lightCurve.time = lightCurveExtracted.time;
	lightCurve.flux = lightCurveExtracted.flux;
	lightCurve.fluxErr = lightCurveExtracted.fluxErr;
	return lightCurve;
end

"""
	libYukiAstronomyTransitSplitLightCurveByPeriod(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitSplitLightCurveByPeriod(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		planetOrbitPeriod::Real,
		planetTransitCentreTime::Real
	)
Split a light curve into segments based on the planet's orbital 
period and the transit centre time. Each segment corresponds to 
one orbital period of the planet.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `transit`: An instance of `libYukiAstronomyTransit` containing 
the planet's orbital period and transit centre time.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitCentreTime`: The time of the transit centre.
# Returns
- A vector of `libYukiAstronomyTransitLightCurve` instances, 
each representing a segment of the original light curve 
corresponding to one orbital period.
"""
libYukiAstronomyTransitSplitLightCurveByPeriod(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitSplitLightCurveByPeriod(
	lightCurve, 
	transit.planet.orbit.period,
	transit.planetTransitCentreTime
);
function libYukiAstronomyTransitSplitLightCurveByPeriod(
	lightCurve::libYukiAstronomyTransitLightCurve, 
	planetOrbitPeriod::Real,
	planetTransitCentreTime::Real
)
	nPeriods = ceil(Int, (maximum(lightCurve.time) - 
		minimum(lightCurve.time)) / planetOrbitPeriod);
	splitLightCurves = Vector{libYukiAstronomyTransitLightCurve}(
		undef, 
		nPeriods
	);
	firstPeriodTime = planetTransitCentreTime - planetOrbitPeriod / 2;
	for index in 1 : nPeriods
		startTime = firstPeriodTime + 
			(index - 1) * planetOrbitPeriod;
		endTime = startTime + planetOrbitPeriod;
		mask = (lightCurve.time .>= startTime) .& 
			(lightCurve.time .< endTime);
		if any(mask)
			splitLightCurves[index] = 
				libYukiAstronomyTransitLightCurve(;
					time = lightCurve.time[mask],
					flux = lightCurve.flux[mask],
					fluxErr = lightCurve.fluxErr[mask]
				);
		else
			splitLightCurves[index] = 
				libYukiAstronomyTransitLightCurve();
		end
	end
	return splitLightCurves;
end

"""
	libYukiAstronomyTransitLoadLightCurveFromFITS(
		fileName::String
	)
Load a light curve from a FITS file.
# Arguments
- `fileName`: The path to the FITS file containing the light curve data.
# Returns
- An instance of `libYukiAstronomyTransitLightCurve` containing the 
time, flux, and flux error data.
"""
function libYukiAstronomyTransitLoadLightCurveFromFITS(
	fileName::String
)
	time, flux, fluxErr = FITS(fileName, "r") do file
		table = file[2];
		return (
			read(table, "TIME"),
			read(table, "FLUX"),
			read(table, "FLUX_ERR")
		);
	end
	return libYukiAstronomyTransitLightCurve(;
		time = time,
		flux = flux,
		fluxErr = fluxErr
	);
end

"""
	libYukiAstronomyTransitSaveLightCurveToFITS(
		lightCurve::libYukiAstronomyTransitLightCurve,
		filePath::String
	)
Save a light curve to a FITS file.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `filePath`: The path to the FITS file where the light curve will 
be saved.
"""
function libYukiAstronomyTransitSaveLightCurveToFITS(
	lightCurve::libYukiAstronomyTransitLightCurve, 
	filePath::String
)
	length(lightCurve.time) == 
		length(lightCurve.flux) == 
		length(lightCurve.fluxErr) ||
        	throw(
				DimensionMismatch(
					"The three vectors must have equal lengths."))
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
		lightCurve::libYukiAstronomyTransitLightCurve
	)
	libYukiAstronomyTransitNormalizeLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve
	)
Normalize the flux of a light curve by dividing it by its median 
value. The flux error is also adjusted accordingly. 
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` with 
normalized flux and flux error.
# Notes
- The function `libYukiAstronomyTransitNormalizeLightCurve!` 
modifies the input light curve in place, while 
`libYukiAstronomyTransitNormalizeLightCurve` returns a new 
normalized light curve.
"""
function libYukiAstronomyTransitNormalizeLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve
)
	lightCurveNormalized = deepcopy(lightCurve);
	libYukiAstronomyTransitNormalizeLightCurve!(lightCurveNormalized);
	return lightCurveNormalized;
end
function libYukiAstronomyTransitNormalizeLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve
)
	medianFlux = median(lightCurve.flux);
	if !isfinite(medianFlux) || medianFlux == 0.0
		error(
			"Normalization failed: median flux is not finite or zero."
		);
	end
	lightCurve.flux ./= medianFlux;
	lightCurve.fluxErr ./= abs(medianFlux);
	return lightCurve;
end

"""
	libYukiAstronomyTransitExtractValidLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve
	)
	libYukiAstronomyTransitExtractValidLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve
	)
Filter out invalid data points from a light curve, ensuring that 
only finite and positive flux error values are retained. The 
filtered light curve is sorted by time.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` containing 
only valid data points.
- Notes
- The function `libYukiAstronomyTransitExtractValidLightCurve!` 
modifies the input light curve in place, while 
`libYukiAstronomyTransitExtractValidLightCurve` returns a new 
filtered light curve.
"""
function libYukiAstronomyTransitExtractValidLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve
)
	lightCurveValid = deepcopy(lightCurve);
	libYukiAstronomyTransitExtractValidLightCurve!(lightCurveValid);
	return lightCurveValid;
end
function libYukiAstronomyTransitExtractValidLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve
)
	validMask = isfinite.(lightCurve.time) .& 
		isfinite.(lightCurve.flux) .& 
		isfinite.(lightCurve.fluxErr) .& 
		(lightCurve.fluxErr .> 0);
	lightCurve.time = lightCurve.time[validMask];
	lightCurve.flux = lightCurve.flux[validMask];
	lightCurve.fluxErr = lightCurve.fluxErr[validMask];

	sortIdx = sortperm(lightCurve.time);
	lightCurve.time = map(
		x -> Float64(x), 
		lightCurve.time[sortIdx]
	);
	lightCurve.flux = map(
		x -> Float64(x), 
		lightCurve.flux[sortIdx]
	);
	lightCurve.fluxErr = map(
		x -> Float64(x), 
		lightCurve.fluxErr[sortIdx]
	);
	return lightCurve;
end

"""
	libYukiAstronomyTransitBinLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		binBoxNumber::Integer
	)
Bin a light curve into a specified number of bins, averaging the 
flux and flux error within each bin. The binned light curve is 
returned with time values corresponding to the centre of each bin.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `binBoxNumber`: The number of bins to divide the light curve into.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` 
containing the binned time, flux, and flux error values.
"""
function libYukiAstronomyTransitBinLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve, 
	binBoxNumber::Integer
) 
	binEdge = range(
		minimum(lightCurve.time), 
		maximum(lightCurve.time); 
		length = binBoxNumber + 1
	);
    binnedTime = @. 0.5 * (binEdge[1 : end - 1] + binEdge[2 : end]);
    binnedFlux = zeros(length(binnedTime));
	binnedFluxErr = zeros(length(binnedTime));
   for index in eachindex(binnedTime)
        mask = (lightCurve.time .>= binEdge[index]) .& 
			(lightCurve.time .< binEdge[index + 1]);
        if any(mask)
            binnedFlux[index] = mean(lightCurve.flux[mask]);
			binnedFluxErr[index] = sqrt(
				sum(lightCurve.fluxErr[mask] .^ 2)
				) / sum(mask);
        else
            binnedFlux[index] = NaN;
			binnedFluxErr[index] = NaN;
        end
    end
    nNaNIndices = findall(x -> !isnan(x), binnedFlux);

    return libYukiAstronomyTransitLightCurve(;
		time = binnedTime[nNaNIndices], 
		flux = binnedFlux[nNaNIndices], 
		fluxErr = binnedFluxErr[nNaNIndices]
	);
end

"""
	libYukiAstronomyTransitFoldLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve,
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitFoldLightCurve(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		planetOrbitPeriod::Real, 
		planetTransitCentreTime::Real
	)
	libYukiAstronomyTransitFoldLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitFoldLightCurve!(
		lightCurve::libYukiAstronomyTransitLightCurve, 
		planetOrbitPeriod::Real, 
		planetTransitCentreTime::Real
	)
Fold a light curve based on the planet's orbital period and the 
transit centre time. The time values are adjusted to be within the 
range of -P/2 to P/2, where P is the planet's orbital period. 
The light curve is then sorted by the adjusted time values.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve`.
- `transit`: An instance of `libYukiAstronomyTransit` containing
the planet's orbital period and transit centre time.
- `planetOrbitPeriod`: The orbital period of the planet.
- `planetTransitCentreTime`: The time of the transit centre.
# Returns
- A new instance of `libYukiAstronomyTransitLightCurve` with 
folded time values.
- Notes
- The function `libYukiAstronomyTransitFoldLightCurve!` modifies 
the input light curve in place, while 
`libYukiAstronomyTransitFoldLightCurve` returns a new folded 
light curve.
"""
libYukiAstronomyTransitFoldLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitFoldLightCurve(
	lightCurve,
	transit.planet.orbit.period,
	transit.planetTransitCentreTime
);
function libYukiAstronomyTransitFoldLightCurve(
	lightCurve::libYukiAstronomyTransitLightCurve,
	planetOrbitPeriod::Real, 
	planetTransitCentreTime::Real
)
	lightCurveFolded = deepcopy(lightCurve);
	libYukiAstronomyTransitFoldLightCurve!(
		lightCurveFolded, 
		planetOrbitPeriod, 
		planetTransitCentreTime
	);
	return lightCurveFolded;
end
function libYukiAstronomyTransitFoldLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	planetOrbitPeriod::Real, 
	planetTransitCentreTime::Real
)
	lightCurve.time = (
		(lightCurve.time .- planetTransitCentreTime) .% 
		planetOrbitPeriod
	);
	lightCurve.time[lightCurve.time .> planetOrbitPeriod / 2] .-= 
		planetOrbitPeriod;

	sortedIndex = sortperm(lightCurve.time);
	lightCurve.time = lightCurve.time[sortedIndex];
	lightCurve.flux = lightCurve.flux[sortedIndex];
	lightCurve.fluxErr = lightCurve.fluxErr[sortedIndex];
	return lightCurve;
end
libYukiAstronomyTransitFoldLightCurve!(
	lightCurve::libYukiAstronomyTransitLightCurve,
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitFoldLightCurve!(
	lightCurve,
	transit.planet.orbit.period,
	transit.planetTransitCentreTime
);

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
		transit::libYukiAstronomyTransit
	)
	libYukiAstronomyTransitTransitDurationEvaluate(
		planetOrbitPeriod::Real, 
		stellarRadius::Real, 
		planetRadius::Real, 
		planetOrbitSemiMajorAxis::Real, 
		planetOrbitInclination::Real
	)
Evaluate the transit duration based on the planet's orbital 
parameters and stellar properties.
# Arguments
- `transit`: An instance of `libYukiAstronomyTransit` containing the
necessary parameters for the evaluation.
- `planetOrbitPeriod`: The orbital period of the planet.
- `stellarRadius`: The radius of the host star (unit consistent 
with `planetRadius`).
- `planetRadius`: The radius of the planet (unit consistent 
with `stellarRadius`).
- `planetOrbitSemiMajorAxis`: The semi-major axis of the planet's 
orbit (unit consistent with `stellarRadius` and `planetRadius`).
- `planetOrbitInclination`: The inclination of the planet's orbit 
in degrees.
# Returns
- The transit duration.
"""
libYukiAstronomyTransitTransitDurationEvaluate(
	transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitTransitDurationEvaluate(
	transit.planet.orbit.period,
	transit.host.radius * 
		libYukiConstantValue(libYukiConstantSolarRadius),
	transit.planet.radius *
		libYukiConstantValue(libYukiConstantEarthRadius),
	transit.planet.orbit.semiMajorAxis *
		libYukiConstantValue(libYukiConstantAstronomicalUnit),
	transit.planet.orbit.inclination
);
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
    impactParameter = scaledSemiMajorAxis * cosd(
		planetOrbitInclination);
    contactRadius = 1 + radiusRatio;
    impactParameter < contactRadius || 
		return zero(Float64);

    argument = sqrt(contactRadius ^ 2 - impactParameter ^ 2) / 
		(scaledSemiMajorAxis * inclinationSin);

    return (planetOrbitPeriod / π * asin(clamp(argument, -1, 1)));
end
