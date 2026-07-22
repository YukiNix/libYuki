using HTTP, DataFrames, Dates, Tables, PyCall, Measurements, CSV, JLD2, FITSIO;

include("libYukiBasic.jl")
include("libYukiConstant.jl")
include("libYukiMath.jl")
include("libYukiPhysics.jl")
include("libYukiAstronomyTransit.jl")

# Kepler Time BJD REFI: BJD Time = Kepler Time + 2454833.0;
const libYukiAstronomyKeplerMissionBJDREFI = 2454833.0;

"""
	libYukiAstronomyKeplerMissionLoadLightCurveFromFITS(
		KeplerID::Int, 
		mastKeplerDownloadPath::String
	)
Load the light curve data for a given Kepler ID from FITS files 
located in the specified directory. The function reads the time, 
flux, and flux error data from the FITS files and returns a 
`libYukiAstronomyTransitLightCurve` instance containing the 
combined data.
# Arguments
- `KeplerID`: The Kepler ID of the target star.
- `mastKeplerDownloadPath`: The path to the directory where the 
Kepler light curve FITS files are stored.
# Returns
- An instance of `libYukiAstronomyTransitLightCurve` containing 
the time, flux, and flux error data for the specified Kepler ID.
# Notes
- The function assumes that the FITS files are named in the format 
`kplr<KeplerID>-<quarter>_llc.fits` and are located in 
subdirectories of the specified `mastKeplerDownloadPath`.
- The time values are adjusted to Barycentric Julian Date (BJD) by 
adding the Kepler time reference value 
(`libYukiAstronomyKeplerMissionBJDREFI`).
- The FITS may cross multiple quarters, and the function will 
combine the data from all relevant FITS files, so that there could 
be some gaps in the time series and some stage between different 
quarters.
"""
function libYukiAstronomyKeplerMissionLoadLightCurveFromFITS(
	KeplerID::Int, 
	mastKeplerDownloadPath::String
)
	KeplerIDStr = lpad(KeplerID, 9, '0')
	targetDirectory = joinpath(
		mastKeplerDownloadPath,
		"kplr$(KeplerIDStr)_lc_Q111111111111111111"
	);

	fitsPaths::Vector{String} = String[];
	for (root, _, files) in walkdir(targetDirectory)
		for filename in files
			if startswith(filename, "kplr$(KeplerIDStr)-") && 
				endswith(filename, "_llc.fits")
				push!(fitsPaths, joinpath(root, filename));
			end
		end
	end
	sort!(fitsPaths);

	lightCurve = libYukiAstronomyTransitLightCurve();
	for fitsPath in fitsPaths
		time, timeCorr, flux, fluxErr = FITS(fitsPath, "r") do file
			table = file[2];
			return (
				read(table, "TIME"),
				read(table, "TIMECORR"),
				read(table, "PDCSAP_FLUX"),
				read(table, "PDCSAP_FLUX_ERR"),
			);
		end
		
		append!(lightCurve.time, time .+ timeCorr .+ 
			libYukiAstronomyKeplerMissionBJDREFI);
		append!(lightCurve.flux, flux);
		append!(lightCurve.fluxErr, fluxErr);
	end
	return lightCurve;
end

"""
	libYukiAstronomyKeplerMissionLoadConfirmedExoplanetInformation(
		saveFilePath::String
	)
Load the converted confirmed exoplanet information from the Kepler 
mission from a JLD2 file.
# Arguments
- `saveFilePath`: The path to the JLD2 file containing the converted 
confirmed exoplanet information.
# Returns
- A DataFrame containing the converted confirmed exoplanet 
information from the Kepler mission.
"""
function libYukiAstronomyKeplerMissionLoadConfirmedExoplanetInformation(
	saveFilePath::String = "KeplerConfirmedPlanets.jld2"
)
	@load saveFilePath convertedExoplanetInformationData;
	return convertedExoplanetInformationData;
end

"""
	libYukiAstronomyKeplerMissionSaveConfirmedExoplanetInformation(
		convertedExoplanetInformationData,
		saveFilePath::String
	)
Save the converted confirmed exoplanet information from the Kepler 
mission to a JLD2 file.
# Arguments
- `convertedExoplanetInformationData`: A DataFrame containing the 
converted confirmed exoplanet information from the Kepler mission.
- `saveFilePath`: The path to the JLD2 file where the data will be 
saved.
"""
function libYukiAstronomyKeplerMissionSaveConfirmedExoplanetInformation(
	convertedExoplanetInformationData, 
	saveFilePath::String = "KeplerConfirmedPlanets.jld2"
)
	@save saveFilePath convertedExoplanetInformationData;
end

"""
	libYukiAstronomyKeplerMissionConvertExoplanetInformation(
		keplerConfirmedPlanetsFrame, 
		keplerConfirmedPlanetsNameFrame
	)
Convert the confirmed exoplanet information from the Kepler mission 
into a structured DataFrame. The function takes two DataFrames as 
input: `keplerConfirmedPlanetsFrame`, which contains the confirmed 
exoplanet data, and `keplerConfirmedPlanetsNameFrame`, which contains 
the names of the confirmed exoplanets. The function processes the 
data and returns a new DataFrame with relevant parameters for each 
confirmed exoplanet, including its name, host star, distance,
orbital period, semi-major axis, eccentricity, radius, mass, and 
other properties.
# Arguments
- `keplerConfirmedPlanetsFrame`: A DataFrame containing the 
confirmed exoplanet data from the Kepler mission.
- `keplerConfirmedPlanetsNameFrame`: A DataFrame containing the 
names of the confirmed exoplanets from the Kepler mission.
# Returns
- A DataFrame containing the converted confirmed exoplanet 
information from the Kepler mission.
"""
function libYukiAstronomyKeplerMissionConvertExoplanetInformation(
	keplerConfirmedPlanetsFrame, 
	keplerConfirmedPlanetsNameFrame
)

	convertedExoplanetInformationData = DataFrame(
		planetName = Union{String, Missing}[],
		hostName = Union{String, Missing}[],
		KeplerID = Union{Int, Missing}[],
		KOIName = Union{String, Missing}[],
		systemDistance = Union{Measurement, Missing}[], # pc
		planetOrbitPeriod = Union{Measurement, Missing}[], # day
		planetOrbitSemiMajorAxis = Union{Measurement, Missing}[], # AU
		planetOrbitEccentricity = Union{Measurement, Missing}[],
		planetRadius = Union{Measurement, Missing}[], # Earth Radius
		stellarRadius = Union{Measurement, Missing}[], # Sun Radius
		planetMass = Union{Measurement, Missing}[], # Earth Mass
		planetMassProvenance = Union{String, Missing}[],
		stellarMass = Union{Measurement, Missing}[], # Sun Mass
		stellarEffectiveTemperature = Union{Measurement, Missing}[], # K
		stellarMetallicity = Union{Measurement, Missing}[], # dex
		stellarMetallicityRatio = Union{String, Missing}[],
		planetOrbitInclination = Union{Measurement, Missing}[], # deg
		planetTransitTimeConjunction = Union{Measurement, Missing}[], # day
		planetTransitTimeReference = Union{String, Missing}[]
	);

	for index in 1 : length(keplerConfirmedPlanetsFrame.pl_name)
		planetName::Union{String, Missing} = 
			keplerConfirmedPlanetsFrame.pl_name[index];
		hostName::Union{String, Missing} = 
			keplerConfirmedPlanetsFrame.hostname[index];
		planetNameIndex = findfirst(
			x -> x == planetName, 
			keplerConfirmedPlanetsNameFrame.pl_name);
		KeplerID::Union{Int, Missing} = 
			isempty(planetNameIndex) ? missing : 
			keplerConfirmedPlanetsNameFrame.kepid[planetNameIndex];
		KOIName::Union{String, Missing} = 
			isempty(planetNameIndex) ? missing : 
			keplerConfirmedPlanetsNameFrame.koi_name[planetNameIndex];
		systemDistance = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.sy_dist[index], 
				keplerConfirmedPlanetsFrame.sy_disterr1[index], 
				keplerConfirmedPlanetsFrame.sy_disterr2[index]
			);
		planetOrbitPeriod = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_orbper[index], 
				keplerConfirmedPlanetsFrame.pl_orbpererr1[index], 
				keplerConfirmedPlanetsFrame.pl_orbpererr2[index]
			);
		planetOrbitSemiMajorAxis = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_orbsmax[index], 
				keplerConfirmedPlanetsFrame.pl_orbsmaxerr1[index], 
				keplerConfirmedPlanetsFrame.pl_orbsmaxerr2[index]
			);
		planetOrbitEccentricity = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_orbeccen[index], 
				keplerConfirmedPlanetsFrame.pl_orbeccenerr1[index], 
				keplerConfirmedPlanetsFrame.pl_orbeccenerr2[index]
			);
		planetRadius = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_rade[index], 
				keplerConfirmedPlanetsFrame.pl_radeerr1[index], 
				keplerConfirmedPlanetsFrame.pl_radeerr2[index]
			);
		stellarRadius = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.st_rad[index], 
				keplerConfirmedPlanetsFrame.st_raderr1[index], 
				keplerConfirmedPlanetsFrame.st_raderr2[index]
			);
		planetMass = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_bmasse[index], 
				keplerConfirmedPlanetsFrame.pl_bmasseerr1[index], 
				keplerConfirmedPlanetsFrame.pl_bmasseerr1[index] + 
					keplerConfirmedPlanetsFrame.pl_bmasseerr2[index]
			);
		planetMassProvenance::Union{String, Missing} = 
			keplerConfirmedPlanetsFrame.pl_bmassprov[index];
		stellarMass = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.st_mass[index], 
				keplerConfirmedPlanetsFrame.st_masserr1[index], 
				keplerConfirmedPlanetsFrame.st_masserr2[index]
			);
		stellarEffectiveTemperature = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.st_teff[index], 
				keplerConfirmedPlanetsFrame.st_tefferr1[index], 
				keplerConfirmedPlanetsFrame.st_tefferr2[index]
			);
		stellarMetallicity = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.st_met[index], 
				keplerConfirmedPlanetsFrame.st_meterr1[index], 
				keplerConfirmedPlanetsFrame.st_meterr2[index]
			);
		stellarMetallicityRatio::Union{String, Missing} = 
			keplerConfirmedPlanetsFrame.st_metratio[index];
		planetOrbitInclination = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_orbincl[index], 
				keplerConfirmedPlanetsFrame.pl_orbinclerr1[index], 
				keplerConfirmedPlanetsFrame.pl_orbinclerr2[index]
			);
		planetTransitTimeConjunction = 
			libYukiBasicMeasurementWithMissing(
				keplerConfirmedPlanetsFrame.pl_tranmid[index], 
				keplerConfirmedPlanetsFrame.pl_tranmiderr1[index], 
				keplerConfirmedPlanetsFrame.pl_tranmiderr2[index]
			);
		planetTransitTimeReference::Union{String, Missing} = 
			keplerConfirmedPlanetsFrame.pl_tranmid_systemref[index];

		push!(convertedExoplanetInformationData, (;
			planetName = planetName, 
			hostName = hostName,
			KeplerID = KeplerID,
			KOIName = KOIName,
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
	return convertedExoplanetInformationData;
end

"""
	libYukiAstronomyKeplerMissionQueryConfirmedExoplanetInformation()
Query the confirmed exoplanet information from the Kepler mission 
using the NASA Exoplanet Archive TAP service. The function retrieves 
various parameters for each confirmed exoplanet, including its 
name, host star, distance, orbital period, semi-major axis, 
eccentricity, radius, mass, and other relevant properties.
# Returns
- A DataFrame containing the confirmed exoplanet information from 
the Kepler mission.
"""
function libYukiAstronomyKeplerMissionQueryConfirmedExoplanetInformation()
	TAPBaseServiceURL = "https://exoplanetarchive.ipac.caltech.edu" * 
		"/TAP/sync?query=";
	TAPDataTypes = "pl_name,hostname,sy_dist,sy_disterr1," * 
		"sy_disterr2,pl_orbper,pl_orbpererr1,pl_orbpererr2," * 
		"pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,pl_orbeccen," * 
		"pl_orbeccenerr1,pl_orbeccenerr2,pl_rade,pl_radeerr1," * 
		"pl_radeerr2,st_rad,st_raderr1,st_raderr2,pl_bmasse," * 
		"pl_bmasseerr1,pl_bmasseerr2,pl_bmassprov,st_mass," * 
		"st_masserr1,st_masserr2,st_teff,st_tefferr1,st_tefferr2," * 
		"st_met,st_meterr1,st_meterr2,st_metratio,pl_orbincl," * 
		"pl_orbinclerr1,pl_orbinclerr2,pl_tranmid,pl_tranmiderr1," * 
		"pl_tranmiderr2,pl_tranmid_systemref";
	TAPDatabaseName = "pscomppars";

	println("#INFO:[" * string(now()) * "] Waiting server response.");
	TAPResponse = HTTP.get(
		TAPBaseServiceURL * "select+" * TAPDataTypes * "+from+" * 
		TAPDatabaseName * "+where+pl_name+like+'Kepler%'&format=csv"
		);
	keplerConfirmedPlanetsData = CSV.File(TAPResponse.body);

	println("#INFO:[" * string(now()) * "] Query " * 
		string(length(keplerConfirmedPlanetsData)) * 
		" confirmed exoplanets from Kepler."
		);

	keplerConfirmedPlanetsFrame = DataFrame(keplerConfirmedPlanetsData);

	return keplerConfirmedPlanetsFrame;
end

"""
	libYukiAstronomyKeplerMissionQueryConfirmedExoplanetName()
Query the names of confirmed exoplanets from the Kepler mission using 
the NASA Exoplanet Archive TAP service. The function retrieves the 
Kepler ID, KOI name, and confirmed exoplanet name for each planet.
# Returns
- A DataFrame containing the Kepler ID, KOI name, and confirmed 
exoplanet name for each planet.
"""
function libYukiAstronomyKeplerMissionQueryConfirmedExoplanetName()
	TAPBaseServiceURL = "https://exoplanetarchive.ipac.caltech.edu" * 
		"/TAP/sync?query=";
	TAPDataTypes = "kepid,koi_name,kepler_name,pl_name";
	TAPDatabaseName = "keplernames";

	println("#INFO:[" * string(now()) * "] Waiting server response.");
	TAPResponse = HTTP.get(TAPBaseServiceURL * "select+" * 
		TAPDataTypes * "+from+" * TAPDatabaseName * 
		"+where+pl_name+like+'Kepler%'&format=csv");
	keplerConfirmedPlanetsNameData = CSV.File(TAPResponse.body);

	println("#INFO:[" * string(now()) * "] Query " * 
		string(length(keplerConfirmedPlanetsNameData)) * 
		" confirmed exoplanets from Kepler.");

	keplerConfirmedPlanetsNameFrame = DataFrame(
		keplerConfirmedPlanetsNameData
		);

	return keplerConfirmedPlanetsNameFrame;
end

"""
	libYukiAstronomyKeplerMissionDownloadLightCurveFITS(
		stellarOriginalName::String, 
		saveFilePath::String
	)
Download Kepler lightcurve by LightKurve(from Python) and save 
to FITS file.
# Arguments
- `stellarOriginalName`: The original name of the star for which 
to download the light curve.
- `saveFilePath`: The path where the downloaded FITS files will 
be saved.
"""
function libYukiAstronomyKeplerMissionDownloadLightCurveFITS(
	stellarOriginalName::String, 
	saveFilePath::String
)
	lightkurve = pyimport("lightkurve")
	searchResult = lightkurve.search_lightcurve(
		stellarOriginalName;
		mission = "Kepler",
		author = "Kepler",
		exptime = "long",
	);

	searchResult.download_all(
		quality_bitmask = "none",
		download_dir = saveFilePath,
	);
end
