"""
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64, 
        inclination::Float64, 
        eccentricity::Float64, 
        argumentOfPeriapsis::Float64, 
        longitudeOfAscendingNode::Float64, 
        meanAnomalyAtEpoch::Float64
    )
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64, 
        inclination::Float64, 
        eccentricity::Float64
    )
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64, 
        inclination::Float64
    )
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64
    )
    libYukiAstronomyOrbit(
        period::Float64
    )
    libYukiAstronomyOrbit()
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
    period::Float64
    semiMajorAxis::Float64
    inclination::Float64
    eccentricity::Float64
    argumentOfPeriapsis::Float64
    longitudeOfAscendingNode::Float64
    meanAnomalyAtEpoch::Float64
    libYukiAstronomyOrbit() = new(
        NaN, NaN, NaN, NaN, NaN, NaN, NaN);
    libYukiAstronomyOrbit(
        period::Float64
    ) = new(period, 
        NaN, NaN, NaN, NaN, NaN, NaN);
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64
    ) = new(period, semiMajorAxis, 
        NaN, NaN, NaN, NaN, NaN);
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64,
        inclination::Float64
    ) = new(period, semiMajorAxis, inclination, 
        NaN, NaN, NaN, NaN);
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64, 
        inclination::Float64, 
        eccentricity::Float64
    ) = new(period, semiMajorAxis, inclination, eccentricity,
        NaN, NaN, NaN);
    libYukiAstronomyOrbit(
        period::Float64, 
        semiMajorAxis::Float64, 
        inclination::Float64, 
        eccentricity::Float64, 
        argumentOfPeriapsis::Float64, 
        longitudeOfAscendingNode::Float64, 
        meanAnomalyAtEpoch::Float64
    ) = new(
        period, semiMajorAxis, inclination, eccentricity, argumentOfPeriapsis, longitudeOfAscendingNode, meanAnomalyAtEpoch
    );
end

"""
    libYukiAstronomyBody(
        name::String, 
        mass::Float64, 
        radius::Float64, 
        orbit::libYukiAstronomyOrbit
    )
    libYukiAstronomyBody(
        name::String, 
        mass::Float64, 
        radius::Float64
    )
    libYukiAstronomyBody()
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
    mass::Float64
    radius::Float64
    orbit::libYukiAstronomyOrbit
    libYukiAstronomyBody() = new("", NaN, NaN, 
        libYukiAstronomyOrbit()
    );
    libYukiAstronomyBody(
        name::String, 
        mass::Float64, 
        radius::Float64
    ) = new(name, mass, radius, 
        libYukiAstronomyOrbit()
    );
    libYukiAstronomyBody(
        name::String, 
        mass::Float64, 
        radius::Float64,
        orbit::libYukiAstronomyOrbit
    ) = new(name, mass, radius, orbit);
end

include("libYukiAstronomyTransit.jl")
include("libYukiAstronomyTransitLimbDarkening.jl")
include("libYukiAstronomyTransitTTV.jl")
include("libYukiAstronomyKeplerMission.jl")
