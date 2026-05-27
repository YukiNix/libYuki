using Measurements

# Convenient constant
const libYukiConstantZero = 0.;
const libYukiConstantOne = 1.;
const libYukiConstant3DZeroVector = [0., 0., 0.];

# Fundamental physics
const libYukiConstantGravitationalConstant = (6.67430 ± 0.00015) * 1.e-11;    # gravitational constant (SI Unit, 2022)
const libYukiConstantDayTime = 86400.; # one day time
const libYukiConstantYearTime = 31556925.9747; # one year time (IUPAC 2025)

# Astrophysics
const libYukiConstantAstronomicalUnit = 1.49597870700e11; # AU (defined)
const libYukiConstantParsec = 3.0856775814913673e16; # parsec (defined)
const libYukiConstantSolarMass = (1.988475 ± 0.000092) * 1.e30;   # solar mass (SI base unit, 2016)
const libYukiConstantEarthMass = (5.9722 ± 0.0006) * 1.e24;   # earth mass (SI base unit, 2016)
const libYukiConstantSolarRadius = (6.95660 ± 0.00140) * 1.e8 # solar radius (SI base unit, 2008)
const libYukiConstantEarthRadiusEquatorial = 6.378100e6 # earth equatorial radius (SI base unit, IAU)
const libYukiConstantEarthRadiusPloar = 6.356800e6 # earth equatorial radius (SI base unit, IAU)
const libYukiConstantEarthRadius = (6.371 ± 0.010 ) * 1.e6 # earth radius (SI base unit, common used)
const libYukiConstantMarsRadius = (3.3895 ± 0.0002) * 1.e6 # mars radius (SI base unit, 2007)
const libYukiConstantJupiterRadius = (6.9886 ± 0.00004) * 1.e7 # jupiter radius (SI base unit, 2026)
const libYukiConstantSaturnRadius = (5.8232 ± 0.0006) * 1.e7 # saturn radius (SI base unit, 2020)

libYukiConstantValue(c) = Measurements.value(c);