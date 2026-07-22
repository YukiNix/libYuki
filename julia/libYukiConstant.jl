using Measurements

# Convenient constant

"""
    libYukiConstantZero
A constant representing the value zero.
"""
const libYukiConstantZero = 0.;

"""
    libYukiConstantOne
A constant representing the value one.
"""
const libYukiConstantOne = 1.;

"""
    libYukiConstant3DZeroVector
A constant representing a 3D vector with all components equal to zero.
"""
const libYukiConstant3DZeroVector = [0., 0., 0.];

"""
    libYukiConstantGravitationalConstant
The gravitational constant (SI Unit, 2022).
"""
const libYukiConstantGravitationalConstant = (6.67430 ± 0.00015) * 
    1.e-11;

"""
    libYukiConstantDayTime
The number of seconds in one day.
"""
const libYukiConstantDayTime = 86400.;

"""
    libYukiConstantYearTime
The number of seconds in one year (IUPAC 2025).
"""
const libYukiConstantYearTime = 31556925.9747; 

"""
    libYukiConstantAstronomicalUnit
A constant representing the astronomical unit (AU).
"""
const libYukiConstantAstronomicalUnit = 1.49597870700e11;

"""
    libYukiConstantParsec
A constant representing the parsec (pc).
"""
const libYukiConstantParsec = 3.0856775814913673e16;

"""
    libYukiConstantSolarMass
A constant representing the solar mass (SI base unit, 2016).
"""
const libYukiConstantSolarMass = (1.988475 ± 0.000092) * 1.e30;

"""
    libYukiConstantEarthMass
A constant representing the Earth mass (SI base unit, 2016).
"""
const libYukiConstantEarthMass = (5.9722 ± 0.0006) * 1.e24;

""" 
    libYukiConstantSolarRadius
A constant representing the solar radius (SI base unit, 2008).
"""
const libYukiConstantSolarRadius = (6.95660 ± 0.00140) * 1.e8;

"""
    libYukiConstantEarthRadiusEquatorial
A constant representing the Earth's equatorial radius 
(SI base unit, IAU).
"""
const libYukiConstantEarthRadiusEquatorial = 6.378100e6;

"""
    libYukiConstantEarthRadiusPolar
A constant representing the Earth's polar radius (SI base unit, IAU).
"""
const libYukiConstantEarthRadiusPloar = 6.356800e6;

"""
    libYukiConstantEarthRadius
A constant representing the Earth's mean radius (SI base unit, IAU).
"""
const libYukiConstantEarthRadius = (6.371 ± 0.010 ) * 1.e6;

"""
    libYukiConstantMarsRadius
A constant representing the Mars radius (SI base unit, 2007).
"""
const libYukiConstantMarsRadius = (3.3895 ± 0.0002) * 1.e6;

"""
    libYukiConstantJupiterRadius
A constant representing the Jupiter radius (SI base unit, 2026).
"""
const libYukiConstantJupiterRadius = (6.9886 ± 0.00004) * 1.e7;

"""
    libYukiConstantSaturnRadius
A constant representing the Saturn radius (SI base unit, 2020).
"""
const libYukiConstantSaturnRadius = (5.8232 ± 0.0006) * 1.e7;

"""
    libYukiConstantValue(
        c::Union{Float64, Measurement{Float64}}
    )
Return the value of a constant `c` as a `Float64`.
# Arguments
- `c`: A constant value, which can be a `Float64` or a 
`Measurement{Float64}`.
# Returns
- The value of the constant `c` as a `Float64`.
"""
libYukiConstantValue(
    c::Union{Float64, Measurement{Float64}}
) = Measurements.value(c);
