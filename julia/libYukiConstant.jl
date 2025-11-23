using Measurements

# Convenient constant
const libYukiConstantZero::Measurement{Float64} = 0. ± 0.;
const libYukiConstantOne::Measurement{Float64} = 1. ± 0.;
const libYukiConstant3DZeroVector::Vector{Measurement{Float64}} = [libYukiConstantZero, libYukiConstantZero, libYukiConstantZero];

# Fundamental physics
const libYukiConstantGravitationalConstant::Measurement{Float64} = (6.67430 ± 0.00015) * 1.e-11;    # gravitational constant (SI Unit, 2022)
const libYukiConstantDayTime::Measurement{Float64} = 24 * 60 * 60. ± 0; # one day time
const libYukiConstantYearTime::Measurement{Float64} = 31556925.9747 ± 0; # one year time (IUPAC 2025)

# Astrophysics
const libYukiConstantAstronomicalUnit::Measurement{Float64} = 1.49597870700 * 1.e11; # AU (defined)
const libYukiConstantSolarMass::Measurement{Float64} = (1.988475 ± 0.000092) * 1.e30;   # solar mass (SI base unit, 2016)
const libYukiConstantEarthMass::Measurement{Float64} = (5.9722 ± 0.0006) * 1.e24;   # earth mass (SI base unit, 2016)
const libYukiConstantSolarRadius::Measurement{Float64} = (6.95660 ± 0.00140) * 1.e8 # solar radius (SI base unit, 2008)
const libYukiConstantEarthRadiusEquatorial::Measurement{Float64} = (6.378100 ± 0.) * 1.e6 # earth equatorial radius (SI base unit, IAU)
const libYukiConstantEarthRadiusPloar::Measurement{Float64} = (6.356800 ± 0.) * 1.e6 # earth equatorial radius (SI base unit, IAU)
const libYukiConstantEarthRadius::Measurement{Float64} = (6.371 ± 0.010 ) * 1.e6 # earth radius (SI base unit, common used)
