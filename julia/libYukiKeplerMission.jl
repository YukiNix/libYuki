using HTTP, DataFrames, Dates, Tables, PyCall, Measurements, CSV, JLD2;

include("libYukiBasic.jl")
include("libYukiTransit.jl")

# Kepler Time BJD REFI: BJD Time = Kepler Time + 2454833.0;
const libYukiKeplerMissionBJDREFI::Float64 = 2454833.0;

# Convert Time of Kepler Lightcurve to BJD.
# Dependency: Measurements.
# TODO: Validate & Example.
function libYukiKeplerMissionConvertLightcurveTimeToBJD(time::Vector{Measurement{Float64}}, flux::Vector{Measurement{Float64}})
	return time .+ libYukiKeplerMissionBJDREFI, flux;
end

# Convert Kepler planets list to convenient format with uncertainty.
# Dependency: DataFrames, Measurements.
# TODO: Validate & Example.
function libYukiKeplerMissionConvertExoplanetInformation(keplerConfirmedPlanetsFrame)

	convertedExoplanetInformationFrame = DataFrame(
		planetName = Union{String, Missing}[],
		hostName = Union{String, Missing}[],
		systemDistance = Union{Measurement{Float64}, Missing}[],					# pc
		planetOrbitPeriod = Union{Measurement{Float64}, Missing}[],					# day
		planetOrbitSemiMajorAxis = Union{Measurement{Float64}, Missing}[],			# AU
		planetOrbitEccentricity = Union{Measurement{Float64}, Missing}[],
		planetRadius = Union{Measurement{Float64}, Missing}[],						# Earth Radius
		stellarRadius = Union{Measurement{Float64}, Missing}[],						# Sun Radius
		planetMass = Union{Measurement{Float64}, Missing}[],						# Earth Mass
		planetMassProvenance = Union{String, Missing}[],
		stellarMass = Union{Measurement{Float64}, Missing}[],						# Sun Mass
		stellarEffectiveTemperature = Union{Measurement{Float64}, Missing}[],		# K
		stellarMetallicity = Union{Measurement{Float64}, Missing}[],				# dex
		stellarMetallicityRatio = Union{String, Missing}[],
		planetOrbitInclination = Union{Measurement{Float64}, Missing}[],			# deg
		planetTransitTimeConjunction = Union{Measurement{Float64}, Missing}[],		# day
		planetTransitTimeReference = Union{String, Missing}[]
	);

	for index in 1 : length(keplerConfirmedPlanetsFrame.pl_name)
		planetName::Union{String, Missing} = keplerConfirmedPlanetsFrame.pl_name[index];
		hostName::Union{String, Missing} = keplerConfirmedPlanetsFrame.hostname[index];
		systemDistance = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.sy_dist[index], keplerConfirmedPlanetsFrame.sy_disterr1[index], keplerConfirmedPlanetsFrame.sy_disterr2[index]);
		planetOrbitPeriod = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_orbper[index], keplerConfirmedPlanetsFrame.pl_orbpererr1[index], keplerConfirmedPlanetsFrame.pl_orbpererr2[index]);
		planetOrbitSemiMajorAxis = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_orbsmax[index], keplerConfirmedPlanetsFrame.pl_orbsmaxerr1[index], keplerConfirmedPlanetsFrame.pl_orbsmaxerr2[index]);
		planetOrbitEccentricity = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_orbeccen[index], keplerConfirmedPlanetsFrame.pl_orbeccenerr1[index], keplerConfirmedPlanetsFrame.pl_orbeccenerr2[index]);
		planetRadius = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_rade[index], keplerConfirmedPlanetsFrame.pl_radeerr1[index], keplerConfirmedPlanetsFrame.pl_radeerr2[index]);
		stellarRadius = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.st_rad[index], keplerConfirmedPlanetsFrame.st_raderr1[index], keplerConfirmedPlanetsFrame.st_raderr2[index]);
		planetMass = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_bmasse[index], keplerConfirmedPlanetsFrame.pl_bmasseerr1[index], keplerConfirmedPlanetsFrame.pl_bmasseerr1[index] + keplerConfirmedPlanetsFrame.pl_bmasseerr2[index]);
		planetMassProvenance::Union{String, Missing} = keplerConfirmedPlanetsFrame.pl_bmassprov[index];
		stellarMass = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.st_mass[index], keplerConfirmedPlanetsFrame.st_masserr1[index], keplerConfirmedPlanetsFrame.st_masserr2[index]);
		stellarEffectiveTemperature = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.st_teff[index], keplerConfirmedPlanetsFrame.st_tefferr1[index], keplerConfirmedPlanetsFrame.st_tefferr2[index]);
		stellarMetallicity = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.st_met[index], keplerConfirmedPlanetsFrame.st_meterr1[index], keplerConfirmedPlanetsFrame.st_meterr2[index]);
		stellarMetallicityRatio::Union{String, Missing} = keplerConfirmedPlanetsFrame.st_metratio[index];
		planetOrbitInclination = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_orbincl[index], keplerConfirmedPlanetsFrame.pl_orbinclerr1[index], keplerConfirmedPlanetsFrame.pl_orbinclerr2[index]);

		planetTransitTimeConjunction = libYukiBasicMeasurementWithMissing(keplerConfirmedPlanetsFrame.pl_tranmid[index], keplerConfirmedPlanetsFrame.pl_tranmiderr1[index], keplerConfirmedPlanetsFrame.pl_tranmiderr2[index]);
		planetTransitTimeReference::Union{String, Missing} = keplerConfirmedPlanetsFrame.pl_tranmid_systemref[index];

		push!(convertedExoplanetInformationFrame, (; 
			planetName = planetName, 
			hostName = hostName,
			systemDistance = systemDistance,
			planetOrbitPeriod = planetOrbitPeriod,
			planetOrbitSemiMajorAxis = planetOrbitSemiMajorAxis,
			planetOrbitEccentricity = planetOrbitEccentricity,
			planetRadius = planetRadius,
			stellarRadius = stellarRadius,
			planetMass = planetMass,
			planetMassProvenance = planetMassProvenance,
			stellarMass = stellarMass,
			stellarEffectiveTemperature = stellarEffectiveTemperature,
			stellarMetallicity = stellarMetallicity,
			stellarMetallicityRatio = stellarMetallicityRatio,
			planetOrbitInclination = planetOrbitInclination,
			planetTransitTimeConjunction = planetTransitTimeConjunction,
			planetTransitTimeReference = planetTransitTimeReference
		));

	end

	return convertedExoplanetInformationFrame;
end

# Load Kepler confirmed planets list.
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiKeplerMissionLoadConfirmedExoplanetInformation()
	@load "KeplerComfirmedPlanets.jld2" keplerConfirmedPlanetsData;
	return keplerConfirmedPlanetsData;
end

# Save Kepler confirmed planets list.
# Dependency: JLD2.
# TODO: Validate & Example.
function libYukiKeplerMissionSaveConfirmedExoplanetInformation(keplerConfirmedPlanetsData)
	@save "KeplerComfirmedPlanets.jld2" keplerConfirmedPlanetsData;
end

# Get Kepler confirmed planets list from NASA Exoplanet Archive.
# Dependency: HTTP, DataFrames, Dates, CSV.
# TODO: Validate & Example.
function libYukiKeplerMissionGetConfirmedExoplanetInfomation()
	TAPBaseServiceURL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=";
    TAPDataTypes = "pl_name,hostname,sy_dist,sy_disterr1,sy_disterr2,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_rade,pl_radeerr1,pl_radeerr2,st_rad,st_raderr1,st_raderr2,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_bmassprov,st_mass,st_masserr1,st_masserr2,st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,st_metratio,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_tranmid_systemref";
    TAPDatabaseName = "pscomppars";

	println("#INFO:[" * string(now()) * "] Waiting server response.");
    TAPResponse = HTTP.get(TAPBaseServiceURL * "select+" * TAPDataTypes * "+from+" * TAPDatabaseName * "+where+pl_name+like+'Kepler%'&format=csv");
    keplerConfirmedPlanetsData = CSV.File(TAPResponse.body);

	println("#INFO:[" * string(now()) * "] Get " * string(length(keplerConfirmedPlanetsData)) * " confirmed exoplanets from Kepler.");

	keplerConfirmedPlanetsFrame = DataFrame(keplerConfirmedPlanetsData);

	return keplerConfirmedPlanetsFrame;
end

# Download Kepler lightcurve by LightKurve(from Python). 
# Dependency: PyCall, HTTP, DataFrames, Dates, Measurements, Tables.
# TODO: Validate & Example.
function libYukiKeplerMissionDownloadLightCurve(stellarOriginalName::String, lightKurve)::Tuple{Vector{Measurement{Float64}}, Vector{Measurement{Float64}}}
	times::Vector{Measurement{Float64}} = [];
	fluxes::Vector{Measurement{Float64}} = [];

	println("#INFO:[" * string(now()) * "] Finding planetary system of " * stellarOriginalName * "...");

	searchResult = lightKurve.search_lightcurve(stellarOriginalName, author="Kepler");

	println("#INFO:[" * string(now()) * "] Found " * string(length(searchResult)) * " sectors/quarters of " * stellarOriginalName * "...");

	println("\t#INFO: Update found, updating...");
	searchResultLength = length(searchResult);
	for index in 1 : searchResultLength
		result = searchResult[index];

		println("\t#INFO:[" * string(now()) * "] Downloading " * string(index) * "/" * string(length(searchResult)) * " curve...");

		try
			subTimes::Vector{Measurement{Float64}} = [];
			subFluxes::Vector{Measurement{Float64}} = [];

			kurve = result.download();
			for (time, flux, flux_err, timecorr) in zip(kurve.time, kurve.flux, kurve.flux_err, kurve.timecorr)
				try
					realFlux = convert(Measurement{Float64}, flux[1] ± flux_err[1]);
					realTime = convert(Measurement{Float64}, time.value + timecorr[1][1]);
					if isnan(realFlux) || isnan(realTime)
						continue;
					end

					append!(subTimes, realTime);
					append!(subFluxes, realFlux);
				catch err
					println(err);
				end
			end
			detrendedTime::Vector{Measurement{Float64}}, detrendedFluxes::Vector{Measurement{Float64}} = libYukiTransitDetrendByWotan(subTimes, subFluxes, 0. ± 0., 0. ± 0., 0. ± 0.);
			append!(times, detrendedTime);
			append!(fluxes, detrendedFluxes);
		catch err
			print("\t\t#ERRO: Curve download faild:");
			println(err);
			continue;
		end
	end

	return times, fluxes; 
end
