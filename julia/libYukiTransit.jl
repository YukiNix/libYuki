using JLD2, Measurements, Dates, Tables, DataFrames;
using Turing, Transits, Random;

# Load lightcurve from file. 
# Dependency: Measurements, JLD2.
function libYukiTransitLoadLightCurveFromJLD2(stellarName)::Tuple{Vector{Measurement{Float64}}, Vector{Measurement{Float64}}}
	@load "$stellarName.jld2" times fluxes
	return times, fluxes;
end

# Save lightcurve to file. 
# Dependency: Measurements, JLD2.
function libYukiTransitSaveLightCurveToJLD2(stellarName::String, times::Vector{Measurement{Float64}}, fluxes::Vector{Measurement{Float64}})
	@save "$stellarName.jld2" times fluxes
end
