using Turing, Orbits, Transits, Random, StableRNGs, JLD2;
using Statistics;

# Quadratic limb-darkening transit flux model.
# Dependency: Turing, Transits, Orbits.
# TODO: Validate & Example.
# @model function libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, f)
# 	u ~ Kipping13()
# 	σ ~ Truncated(Cauchy(0, 0.01), 0, Inf)
# 	fModel = libYukiAstronomyTransitQuadraticLimbDarkeningTransitFlux(t, orbitPeriod, transitDuration, planetStellarRadiusRatio, collect(u))
# 	for index in eachindex(t)
# 		f[index] ~ Normal(fModel[index], σ)
# 	end
# end

@model function libYukiAstronomyTransitLimbDarkeningDurationTransitFluxModel(
    time::AbstractArray{<:Real},
    orbitPeriod::Real,
    planetTransitDepth::Real,
    limbDarkeningFunc::AbstractLimbDark,
    flux::AbstractArray{<:Real}
)

end

"""
# Model
    libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
        lightCurve::libYukiAstronomyTransitLightCurve,
        transit::libYukiAstronomyTransit
    )
    libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
        time::AbstractArray{<:Real}, 
        orbitPeriod::Real, 
        planetTransitDuration::Real, 
        planetTransitDepth::Real, 
        flux::AbstractArray{<:Real}
    )
A Turing model for fitting transit flux data using a quadratic 
limb-darkening model. The model takes time, orbit period, 
transit duration, planet-to-stellar radius ratio, observed 
flux, and optional per-point flux uncertainty as inputs. It 
estimates limb-darkening coefficients and an additional jitter term.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing time and flux data.
- `transit`: An instance of `libYukiAstronomyTransit` containing 
transit parameters.
- `time`: A vector of time values at which the transit flux is observed.
- `orbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `planetTransitDepth`: The depth of the transit, representing the fractional decrease in flux during the transit.
- `flux`: A vector of observed flux values corresponding to the 
time values.
# Returns
- A Turing model for Bayesian inference of limb-darkening and 
noise parameters.
"""
libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
    lightCurve::libYukiAstronomyTransitLightCurve,
    transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
        lightCurve.time,
        transit.planet.orbit.period,
        transit.planetTransitDuration,
        transit.planetTransitDepth,
        lightCurve.flux;
        fluxErr = lightCurve.fluxErr
);
@model function libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
    time::AbstractArray{<:Real},
    orbitPeriod::Real,
    planetTransitDuration::Real,
    planetTransitDepth::Real,
    flux::AbstractArray{<:Real};
    fluxErr::AbstractArray{<:Real} = fill(0.01, length(flux))
)
    # Kipping (2013) reparameterization for quadratic limb darkening.
    q1 ~ Uniform(0, 1 - 1e-10)
    q2 ~ Uniform(0, 1 - 1e-10)
    logJitter ~ Normal(-6, 2)

    u1 = 2 * sqrt(q1) * q2;
    u2 = sqrt(q1) * (1 - 2 * q2);
    totalErr = sqrt.(fluxErr.^2 .+ exp(logJitter)^2);

    orbit, limbDarkeningFunc = 
        libYukiAstronomyTransitLimbDarkeningQuadratic(
            orbitPeriod, 
            planetTransitDuration, 
            [u1, u2]
        );
    modelLightCurve = libYukiAstronomyTransitLimbDarkeningTransitFlux(
        time, 
        planetTransitDepth,
        orbit, 
        limbDarkeningFunc
    );

    flux ~ MvNormal(modelLightCurve.flux, totalErr)
end

"""
    libYukiAstronomyTransitLimbDarkeningTransitFlux!(
        lightCurve::libYukiAstronomyTransitLightCurve,
        transit::libYukiAstronomyTransit
    )
    libYukiAstronomyTransitLimbDarkeningTransitFlux!(
        lightCurve::libYukiAstronomyTransitLightCurve, 
        planetTransitDepth::Real, 
        orbit::Orbits.AbstractOrbit, 
        limbDarkeningFunc::AbstractLimbDark
    )
    libYukiAstronomyTransitLimbDarkeningTransitFlux(
        time::AbstractVector{<:Real}, 
        planetTransitDepth::Real, 
        orbit::Orbits.AbstractOrbit, 
        limbDarkeningFunc::AbstractLimbDark
    )
Compute the transit flux for a given light curve and transit 
parameters, optionally using a specified limb darkening function. 
The results are stored in the provided light curve object.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing time and flux data.
- `transit`: An instance of `libYukiAstronomyTransit` containing 
transit parameters.
- `planetTransitDepth`: The depth of the transit, representing the fractional decrease in flux during the transit.
- `orbit`: An instance of `Orbits.AbstractOrbit` representing the 
orbital parameters of the planet.
- `limbDarkeningFunc`: An instance of `AbstractLimbDark` representing 
the limb darkening function to be used in the transit flux calculation.
# Returns
- The updated `lightCurve` object with computed flux values based 
on the transit parameters and limb darkening function.
"""
libYukiAstronomyTransitLimbDarkeningTransitFlux!(
    lightCurve::libYukiAstronomyTransitLightCurve,
    transit::libYukiAstronomyTransit
) = libYukiAstronomyTransitLimbDarkeningTransitFlux!(
    lightCurve, 
    transit.planetTransitDepth,
    SimpleOrbit(
        period = transit.planet.orbit.period,
        duration = transit.planetTransitDuration
    ),
    transit.limbDarkeningFunc
);
function libYukiAstronomyTransitLimbDarkeningTransitFlux!(
    lightCurve::libYukiAstronomyTransitLightCurve, 
    planetTransitDepth::Real, 
    orbit::Orbits.AbstractOrbit, 
    limbDarkeningFunc::AbstractLimbDark
)
    lightCurveGenerated = 
        libYukiAstronomyTransitLimbDarkeningTransitFlux(
            lightCurve.time, 
            planetTransitDepth, 
            orbit, 
            limbDarkeningFunc
        );
    lightCurve.flux = lightCurveGenerated.flux;
    return lightCurve;
end
function libYukiAstronomyTransitLimbDarkeningTransitFlux(
    time::AbstractVector{<:Real}, 
    planetTransitDepth::Real, 
    orbit::Orbits.AbstractOrbit, 
    limbDarkeningFunc::AbstractLimbDark
)
	return libYukiAstronomyTransitLightCurve(
        time = time,
        flux = limbDarkeningFunc.(
            Ref(orbit), 
            time, 
            sqrt(planetTransitDepth)
        )
    );
end

"""
    libYukiAstronomyTransitLimbDarkeningPolynomial(
        transit::libYukiAstronomyTransit,
        limbPolynomialDarkeningParameters::AbstractVector{<:Real}
    )
    libYukiAstronomyTransitLimbDarkeningPolynomial(
        orbitPeriod::Real, 
        planetTransitDuration::Real, 
        limbPolynomialDarkeningParameters::AbstractVector{<:Real}
    )
Create a simple orbit and polynomial limb darkening model for 
    transit light curve analysis.
# Arguments
- `transit`: An instance of `libYukiAstronomyTransit` containing 
transit parameters.
- `orbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `limbPolynomialDarkeningParameters`: A vector of coefficients for 
the polynomial limb darkening model.
# Returns
- A tuple containing a `SimpleOrbit` instance and a 
`PolynomialLimbDark` instance.
"""
function libYukiAstronomyTransitLimbDarkeningPolynomial(
    transit::libYukiAstronomyTransit,
    limbPolynomialDarkeningParameters::AbstractVector{<:Real}
)
    orbit, transit.limbDarkeningFunc = 
        libYukiAstronomyTransitLimbDarkeningPolynomial(
            transit.planet.orbit.period,
            transit.planetTransitDuration,
            limbPolynomialDarkeningParameters
        );
    return orbit, transit.limbDarkeningFunc;
end
function libYukiAstronomyTransitLimbDarkeningPolynomial(
    orbitPeriod::Real, 
    planetTransitDuration::Real, limbPolynomialDarkeningParameters::AbstractVector{<:Real}
) 
    return SimpleOrbit(
            period = orbitPeriod, 
            duration = planetTransitDuration
        ), 
        PolynomialLimbDark(limbPolynomialDarkeningParameters);
end

"""
    libYukiAstronomyTransitLimbDarkeningQuadratic(
        transit::libYukiAstronomyTransit,
        limbQuadraticDarkeningParameters::AbstractVector{<:Real}
    )
    libYukiAstronomyTransitLimbDarkeningQuadratic(
        transit::libYukiAstronomyTransit,
        chain::Chains
    )
    libYukiAstronomyTransitLimbDarkeningQuadratic(
        orbitPeriod::Real, 
        planetTransitDuration::Real, 
        limbQuadraticDarkeningParameters::AbstractVector{<:Real}
    )
    libYukiAstronomyTransitLimbDarkeningQuadratic(
        orbitPeriod::Real, 
        planetTransitDuration::Real, 
        chain::Chains
    )
Create a simple orbit and quadratic limb darkening model for 
transit light curve analysis.
# Arguments
- `transit`: An instance of `libYukiAstronomyTransit` containing 
transit parameters.
- `orbitPeriod`: The orbital period of the planet.
- `planetTransitDuration`: The duration of the transit.
- `limbQuadraticDarkeningParameters`: A vector of coefficients for 
the quadratic limb darkening model.
- `chain`: A `Chains` object containing MCMC samples from which
the limb darkening coefficients will be derived.
# Returns
- A tuple containing a `SimpleOrbit` instance and a `QuadLimbDark` 
instance.
"""
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    transit::libYukiAstronomyTransit,
    chain::Chains
)
    orbit, transit.limbDarkeningFunc = 
        libYukiAstronomyTransitLimbDarkeningQuadratic(
            transit.planet.orbit.period,
            transit.planetTransitDuration,
            chain
        );
    return orbit, transit.limbDarkeningFunc;
end
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    transit::libYukiAstronomyTransit,
    limbQuadraticDarkeningParameters::AbstractVector{<:Real}
)
    orbit, transit.limbDarkeningFunc =  
        libYukiAstronomyTransitLimbDarkeningQuadratic(
            transit.planet.orbit.period,
            transit.planetTransitDuration,
            limbQuadraticDarkeningParameters
        );
    return orbit, transit.limbDarkeningFunc;
end
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    orbitPeriod::Real, 
    planetTransitDuration::Real, 
    chain::Chains
)
    chainNames = Set(Symbol.(names(chain, :parameters)));
    if :u1 in chainNames && :u2 in chainNames
        u1 = median(chain[:u1]);
        u2 = median(chain[:u2]);
    elseif :q1 in chainNames && :q2 in chainNames
        q1 = median(chain[:q1]);
        q2 = median(chain[:q2]);
        u1 = 2 * sqrt(q1) * q2;
        u2 = sqrt(q1) * (1 - 2 * q2);
    else
        throw(ArgumentError("Chains must contain either" * 
            " (:u1, :u2) or (:q1, :q2)."))
    end

    return libYukiAstronomyTransitLimbDarkeningQuadratic(
        orbitPeriod, 
        planetTransitDuration,
        [u1, u2]
    );
end
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    orbitPeriod::Real, 
    planetTransitDuration::Real, 
    limbQuadraticDarkeningParameters::AbstractVector{<:Real}
) 
    return SimpleOrbit(
            period = orbitPeriod, 
            duration = planetTransitDuration
        ), 
        QuadLimbDark(limbQuadraticDarkeningParameters);
end
