using JLD2, Measurements, Dates, Tables, DataFrames, ProgressMeter;
using Turing, Transits, Random, Orbits;
using PyCall;

include("libYukiBasic.jl")

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

# Transit flux with limb-darkening.
# Dependency: Transits, Orbits.
# TODO: Validate & Example.
function libYukiTransitQuadraticLimbDarkeningTransitFlux(times::Vector{T}, orbitPeriod::T, transitDuration::T, planetStellarRadiusRatio::T, limbQuadraticDarkeningParameters::Vector{T}) where {T <: AbstractFloat}
	orbit = SimpleOrbit(period = orbitPeriod, duration = transitDuration);
    ld = QuadLimbDark(limbQuadraticDarkeningParameters);
	return ld.(Ref(orbit), times, planetStellarRadiusRatio)
end

# Transit flux with polynomial limb-darkening.
# Dependency: Transits, Orbits.
# TODO: Validate & Example.
function libYukiTransitPolynomialLimbDarkeningTransitFlux(times::Vector{T}, orbitPeriod::T, transitDuration::T, planetStellarRadiusRatio::T, limbDarkeningParameters::Vector{T}) where {T <: AbstractFloat}
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
function libYukiTransitBinLightCurve(times::Vector{T}, fluxes::Vector{T}, binBoxNumber::Integer) where {T <: AbstractFloat}
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
function libYukiTransitFoldLightCurve(times::Vector{T}, fluxes::Vector{T}, planetOrbitPeriod::T, transitCentreTime::T) where {T <: AbstractFloat}
    foldedTimes = ((times .- transitCentreTime) .% planetOrbitPeriod);
    foldedTimes[foldedTimes .> planetOrbitPeriod / 2] .-= planetOrbitPeriod;

	sortedIndex = sortperm(foldedTimes);
    foldedTimes = foldedTimes[sortedIndex];
	foldedFluxes = fluxes[sortedIndex];
	return foldedTimes, foldedFluxes;
end

# Detrending lightcurve with Wotan.
# Dependency: PyCall, wotan(python package).
# TODO: Validate & Example.
function libYukiTransitDetrendByWotan(times::Vector{T}, fluxes::Vector{T}, stellarRadius::T, stellarMass::T, planetOrbitPeriod::T) where {T <: AbstractFloat}
	pyWotan = pyimport("wotan");
	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
		window = 0.75;
	else
		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
	end
	sortedTimes, sortedFluxes = libYukiBasicSortElementsByOrderVector(times, fluxes);
	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
	return sortedTimes, Vector{T}(detrendedFluxes);
end

# Detrending lightcurve with Wotan (with uncertainty propagation).
# Dependency: Measurements, PyCall, wotan(python package).
# TODO: Validate & Example.
function libYukiTransitDetrendByWotan(times::Vector{Measurement{T}}, fluxes::Vector{Measurement{T}}, stellarRadius::Measurement{T}, stellarMass::Measurement{T}, planetOrbitPeriod::Measurement{T}) where {T <: AbstractFloat}
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
function libYukiTransitSaveLightCurveToJLD2(stellarName::String, times::Vector{T}, fluxes::Vector{T}) where {T <: AbstractFloat}
	@save "$stellarName.jld2" times fluxes
end
