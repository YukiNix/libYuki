"""
    libYukiAstronomyOrbit(;
        period::Real, 
        semiMajorAxis::Real, 
        inclination::Real, 
        eccentricity::Real, 
        argumentOfPeriapsis::Real, 
        longitudeOfAscendingNode::Real, 
        meanAnomalyAtEpoch::Real
    )
Create a new `libYukiAstronomyOrbit` instance with the specified 
parameters.
# Arguments
- `period`: Orbital period in days.
- `semiMajorAxis`: Semi-major axis in astronomical units (AU).
- `inclination`: Orbital inclination in degrees.
- `eccentricity`: Orbital eccentricity (0 for circular, 0 < e < 1 
for elliptical).
- `argumentOfPeriapsis`: Argument of periapsis in degrees.
- `longitudeOfAscendingNode`: Longitude of ascending node in degrees.
- `meanAnomalyAtEpoch`: Mean anomaly at epoch in degrees.
# Returns
- A new instance of `libYukiAstronomyOrbit`.
"""
mutable struct libYukiAstronomyOrbit
    period::Real
    semiMajorAxis::Real
    inclination::Real
    eccentricity::Real
    argumentOfPeriapsis::Real
    longitudeOfAscendingNode::Real
    meanAnomalyAtEpoch::Real
    libYukiAstronomyOrbit(;
        period::Real = NaN, 
        semiMajorAxis::Real = NaN, 
        inclination::Real = NaN, 
        eccentricity::Real = NaN, 
        argumentOfPeriapsis::Real = NaN, 
        longitudeOfAscendingNode::Real = NaN, 
        meanAnomalyAtEpoch::Real = NaN
    ) = new(
        period, 
        semiMajorAxis, 
        inclination, 
        eccentricity, 
        argumentOfPeriapsis, 
        longitudeOfAscendingNode, 
        meanAnomalyAtEpoch
    );
end

"""
    libYukiAstronomyBody(;
        name::String, 
        mass::Real, 
        radius::Real, 
        orbit::libYukiAstronomyOrbit
    )
Create a new `libYukiAstronomyBody` instance with the specified 
parameters.
# Arguments
- `name`: Name of the astronomical body.
- `mass`: Mass of the body in kilograms.
- `radius`: Radius of the body in meters.
- `orbit`: An instance of `libYukiAstronomyOrbit` representing 
the body's orbit.
# Returns
- A new instance of `libYukiAstronomyBody`.
"""
mutable struct libYukiAstronomyBody
    name::String
    mass::Real
    radius::Real
    orbit::libYukiAstronomyOrbit
    libYukiAstronomyBody(;
        name::String = "", 
        mass::Real = NaN, 
        radius::Real = NaN,
        orbit::libYukiAstronomyOrbit = libYukiAstronomyOrbit()
    ) = new(name, mass, radius, orbit);
end

include("libYukiAstronomyKeplerian.jl")
include("libYukiAstronomyTransit.jl")
include("libYukiAstronomyTransitLimbDarkening.jl")
include("libYukiAstronomyTransitVariation.jl")
include("libYukiAstronomyKeplerMission.jl")
