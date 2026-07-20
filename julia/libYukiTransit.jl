using JLD2, Measurements, Dates, Tables, DataFrames, ProgressMeter;
using Turing, Transits, Random, Orbits, StableRNGs;
using PyCall;
using FITSIO;

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysics.jl")

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

function libYukiTransitNormalizeLightCurve!(lightCurve::libYukiTransitLightCurve)
	medianFlux = median(lightCurve.flux);
	if !isfinite(medianFlux) || medianFlux == 0.0
		error("Normalization failed: median flux is not finite or zero.")
	end
	lightCurve.flux ./= medianFlux;
	lightCurve.fluxErr ./= abs(medianFlux);
end

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

# Folding lightcurve.
# TODO: Example.
function libYukiTransitFoldLightCurve(time::AbstractVector{<:Real}, flux::AbstractVector{<:Real}, fluxErr::AbstractVector{<:Real}, planetOrbitPeriod::Real, transitCentreTime::Real)
    foldedTime = ((time .- transitCentreTime) .% planetOrbitPeriod);
    foldedTime[foldedTime .> planetOrbitPeriod / 2] .-= planetOrbitPeriod;

	sortedIndex = sortperm(foldedTime);
    foldedTime = foldedTime[sortedIndex];
	foldedFlux = flux[sortedIndex];
	foldedFluxErr = fluxErr[sortedIndex];

	return libYukiTransitLightCurve(foldedTime, foldedFlux, foldedFluxErr);
end
function libYukiTransitFoldLightCurve(lightCurve::libYukiTransitLightCurve, planetOrbitPeriod::Real, transitCentreTime::Real)
	foldedTime = ((lightCurve.time .- transitCentreTime) .% planetOrbitPeriod);
	foldedTime[foldedTime .> planetOrbitPeriod / 2] .-= planetOrbitPeriod;

	sortedIndex = sortperm(foldedTime);
	foldedTime = foldedTime[sortedIndex];
	foldedFlux = lightCurve.flux[sortedIndex];
	foldedFluxErr = lightCurve.fluxErr[sortedIndex];

	return libYukiTransitLightCurve(foldedTime, foldedFlux, foldedFluxErr);
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

# Detrending lightcurve with Wotan (with uncertainty propagation).
# Dependency: Measurements, PyCall, wotan(python package).
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
	return sortedTimes, detrendedFluxes .± (Measurements.uncertainty.(sortedFluxes) ./ Measurements.value.(sortedFluxes));
end

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

# Load lightcurve from file. 
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiTransitLoadLightCurveFromJLD2(stellarName)
	@load "$stellarName.jld2" times fluxes
	return libYukiBasicSortElementsByOrderVector(times, fluxes);
end

# Save lightcurve to file. 
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiTransitSaveLightCurveToJLD2(stellarName::String, times::AbstractVector{<:Real}, fluxes::AbstractVector{<:Real}) 
	@save "$stellarName.jld2" times fluxes
end
