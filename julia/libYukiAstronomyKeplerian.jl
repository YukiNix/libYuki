"""
    libYukiAstronomyKeplerianOrbitSemiMajorAxis(
        stellarMass::Real,
        planetOrbitPeriod::Real
    )
Compute the semi-major axis of a planet's orbit around a star using 
Kepler's third law.
# Arguments
- `stellarMass`: The mass of the star in solar masses.
- `planetOrbitPeriod`: The orbital period of the planet in days.
# Returns
- The semi-major axis of the planet's orbit in astronomical units.
"""
function libYukiAstronomyKeplerianOrbitSemiMajorAxis(
    stellarMass::Real,
    planetOrbitPeriod::Real
)
    stellarMassSIUnit = stellarMass * libYukiConstantValue(libYukiConstantSolarMass);
    planetOrbitPeriodSIUnit = planetOrbitPeriod * libYukiConstantDayTime;
    omega = 2 * π / planetOrbitPeriodSIUnit;
    semiMajorAxisSIUnit = (
        libYukiConstantValue(libYukiConstantGravitationalConstant) * 
            stellarMassSIUnit / 
        (omega ^ 2)
    ) ^ (1 / 3);
    return semiMajorAxisSIUnit / libYukiConstantValue(libYukiConstantAstronomicalUnit);
end
