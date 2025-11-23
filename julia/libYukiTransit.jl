using JLD2, Measurements, Dates, Tables, DataFrames;
using Turing, Transits, Random;
using PyCall;

include("libYukiBasic.jl")

# Detrending lightcurve with Wotan\.
# Dependency: Measurements, PyCall, wotan(python package).
# TODO: Validate & Example.
function libYukiTransitDetrendByWotan(times::Vector{Measurement{Float64}}, fluxes::Vector{Measurement{Float64}}, stellarRadius::Measurement{Float64}, stellarMass::Measurement{Float64}, planetOrbitPeriod::Measurement{Float64})::Tuple{Vector{Measurement{Float64}}, Vector{Measurement{Float64}}}
	pyWotan = pyimport("wotan");
	if stellarRadius == 0. || stellarMass == 0. || planetOrbitPeriod == 0.
		window = 0.75;
	else
		window = pyWotan.t14(R_s = Measurements.value(stellarRadius), M_s = Measurements.value(stellarMass), P = Measurements.value(planetOrbitPeriod));
	end
	sortedTimes::Vector{Measurement{Float64}}, sortedFluxes::Vector{Measurement{Float64}} = libYukiBasicSort2VectorsByTheFirstOne(times, fluxes);
	detrendedFluxes, _ = pyWotan.flatten(Measurements.value.(sortedTimes), Measurements.value.(sortedFluxes), window_length = window, method = "biweight", return_trend = true);
	return sortedTimes, detrendedFluxes .± (Measurements.uncertainty.(sortedFluxes) ./ Measurements.value.(sortedFluxes));
end

# Load lightcurve from file. 
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiTransitLoadLightCurveFromJLD2(stellarName)::Tuple{Vector{Measurement{Float64}}, Vector{Measurement{Float64}}}
	@load "$stellarName.jld2" times fluxes
	return times, fluxes;
end

# Save lightcurve to file. 
# Dependency: Measurements, JLD2.
# TODO: Validate & Example.
function libYukiTransitSaveLightCurveToJLD2(stellarName::String, times::Vector{Measurement{Float64}}, fluxes::Vector{Measurement{Float64}})
	@save "$stellarName.jld2" times fluxes
end
