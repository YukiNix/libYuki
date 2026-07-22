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

"""
# Model
    libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
        time, 
        orbitPeriod, 
        transitDuration, 
        planetStellarRadiusRatio, 
        flux,
        fluxErr = nothing
    )
A Turing model for fitting transit flux data using a quadratic 
limb-darkening model. The model takes time, orbit period, 
transit duration, planet-to-stellar radius ratio, observed 
flux, and optional per-point flux uncertainty as inputs. It 
estimates limb-darkening coefficients and an additional jitter term.
# Arguments
- `time`: A vector of time values at which the transit flux is observed.
- `orbitPeriod`: The orbital period of the planet.
- `transitDuration`: The duration of the transit.
- `planetStellarRadiusRatio`: The ratio of the planet's radius to 
the star's radius.
- `flux`: A vector of observed flux values corresponding to the 
time values.
- `fluxErr`: Optional per-point flux uncertainty. If provided, the 
model uses `sqrt(fluxErr[i]^2 + sigmaJitter^2)` as observational sigma.
# Returns
- A Turing model for Bayesian inference of limb-darkening and 
noise parameters.
"""
@model function libYukiAstronomyTransitLimbDarkeningQuadraticTransitFluxModel(
    time::AbstractArray{<:Real},
    orbitPeriod::Real,
    transitDuration::Real,
    planetStellarRadiusRatio::Real,
    flux::AbstractArray{<:Real}
)
    # Kipping (2013) reparameterization for quadratic limb darkening.
    q1 ~ Uniform(0, 1 - 1e-10)
    q2 ~ Uniform(0, 1 - 1e-10)
    u1 = 2 * sqrt(q1) * q2;
    u2 = sqrt(q1) * (1 - 2 * q2);

    orbit, limbDarkeningFunc = 
        libYukiAstronomyTransitLimbDarkeningQuadratic(
            orbitPeriod, 
            transitDuration, 
            [u1, u2]
        );
    fModel = libYukiAstronomyTransitLimbDarkeningTransitFlux(
        time, 
        planetStellarRadiusRatio, 
        orbit, 
        limbDarkeningFunc
    );
    for index in eachindex(time)
        flux[index] ~ Normal(fModel[index], 0.01)
    end
end

"""
    libYukiAstronomyTransitLimbDarkeningTransitFluxModelMCMCSample(
        transitFluxModel, 
        sampler = NUTS(), 
        modelTuringSamples::Int = 1000
    )
Perform MCMC sampling on the given transit flux model using the 
specified sampler and number of samples.
# Arguments
- `transitFluxModel`: A Turing model representing the transit flux.
- `sampler`: The MCMC sampler to use (default is NUTS).
- `modelTuringSamples`: The number of samples to draw from the posterior distribution (default is 1000).
# Returns
- A `Chains` object containing the MCMC samples from the posterior 
distribution of the model parameters.
"""
function libYukiAstronomyTransitLimbDarkeningTransitFluxModelMCMCSample(
    transitFluxModel, 
    sampler = NUTS(), 
    modelTuringSamples::Int = 1000
)
    modelTuringSamples > 0 || 
    throw(ArgumentError("modelTuringSamples must be a " * 
        "positive integer."))
    return sample(
        transitFluxModel, 
        sampler, 
        modelTuringSamples
    );
end

"""
    libYukiAstronomyTransitLimbDarkeningLoadMCMCChains(
        saveFilePath = "TransitLimbDarkeningMCMCChains.jld2"
    )
Load the MCMC chains for the transit flux model from a JLD2 file.
# Arguments
- `saveFilePath`: The path to the JLD2 file containing the MCMC 
chains (default is "TransitLimbDarkeningMCMCChains.jld2").
# Returns
- A `Chains` object containing the MCMC samples from the 
posterior distribution of the model parameters.
"""
function libYukiAstronomyTransitLimbDarkeningLoadMCMCChains(
    saveFilePath::String = "TransitLimbDarkeningMCMCChains.jld2"
)
    @load saveFilePath transitFluxModelChains;
    return transitFluxModelChains;
end

"""
    libYukiAstronomyTransitLimbDarkeningSaveMCMCChains(
        transitFluxModelChains, 
        saveFilePath = "TransitLimbDarkeningMCMCChains.jld2"
    )
Save the MCMC chains for the transit flux model to a JLD2 file.
# Arguments
- `transitFluxModelChains`: A `Chains` object containing the MCMC 
samples from the posterior distribution of the model parameters.
- `saveFilePath`: The path to the JLD2 file where the chains will 
be saved (default is "TransitLimbDarkeningMCMCChains.jld2").
"""
function libYukiAstronomyTransitLimbDarkeningSaveMCMCChains(
    transitFluxModelChains::Chains, 
    saveFilePath::String = "TransitLimbDarkeningMCMCChains.jld2"
)
    @save saveFilePath transitFluxModelChains;
end

"""
    libYukiAstronomyTransitLimbDarkeningTransitFlux(
        time, 
        planetStellarRadiusRatio, 
        orbit, 
        limbDarkeningFunc
    )
    libYukiAstronomyTransitLimbDarkeningTransitFlux!(
        lightCurve, 
        planetStellarRadiusRatio, 
        orbit, 
        limbDarkeningFunc
    )
Calculate the transit flux for a given set of time, planet-to-stellar 
radius ratio, orbit, and limb darkening function. The function can 
be called with either a vector of time or a 
`libYukiAstronomyTransitLightCurve` instance. In the latter case, 
the flux values in the light curve will be updated in place.
# Arguments
- `time`: A vector of time values at which to calculate the transit 
flux.
- `planetStellarRadiusRatio`: The ratio of the planet's radius to 
the star's radius.
- `orbit`: An instance of a subtype of `AbstractOrbit` representing 
the orbital parameters of the planet.
- `limbDarkeningFunc`: An instance of a subtype of `AbstractLimbDark` 
representing the limb darkening model to be used for the transit 
flux calculation.
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` 
containing time values and flux values to be updated.
# Returns
- A vector of calculated transit flux values corresponding to the 
input time.
- Or updates the flux values in the provided 
`libYukiAstronomyTransitLightCurve` instance in place.
"""
function libYukiAstronomyTransitLimbDarkeningTransitFlux(
    time::AbstractVector{<:Real}, 
    planetStellarRadiusRatio::Real, 
    orbit::Orbits.AbstractOrbit, 
    limbDarkeningFunc::AbstractLimbDark
)
	return limbDarkeningFunc.(
        Ref(orbit), 
        time, 
        planetStellarRadiusRatio
    );
end
function libYukiAstronomyTransitLimbDarkeningTransitFlux!(
    lightCurve::libYukiAstronomyTransitLightCurve, 
    planetStellarRadiusRatio::Real, 
    orbit::Orbits.AbstractOrbit, 
    limbDarkeningFunc::AbstractLimbDark
)
    lightCurve.flux .= libYukiAstronomyTransitLimbDarkeningTransitFlux(
        lightCurve.time, 
        planetStellarRadiusRatio, 
        orbit, 
        limbDarkeningFunc
    );
    return lightCurve;
end

"""
    libYukiAstronomyTransitLimbDarkeningPolynomial(
        orbitPeriod::Real, 
        transitDuration::Real, 
        limbPolynomialDarkeningParameters::AbstractVector{<:Real}
    )
Create a simple orbit and polynomial limb darkening model for 
    transit light curve analysis.
# Arguments
- `orbitPeriod`: The orbital period of the planet.
- `transitDuration`: The duration of the transit.
- `limbPolynomialDarkeningParameters`: A vector of coefficients for 
the polynomial limb darkening model.
# Returns
- A tuple containing a `SimpleOrbit` instance and a 
`PolynomialLimbDark` instance.
"""
function libYukiAstronomyTransitLimbDarkeningPolynomial(
    orbitPeriod::Real, 
    transitDuration::Real, limbPolynomialDarkeningParameters::AbstractVector{<:Real}
) 
    return SimpleOrbit(
            period = orbitPeriod, 
            duration = transitDuration
        ), 
        PolynomialLimbDark(limbPolynomialDarkeningParameters);
end

"""
    libYukiAstronomyTransitLimbDarkeningQuadratic(
        orbitPeriod::Real, 
        transitDuration::Real, 
        limbQuadraticDarkeningParameters::AbstractVector{<:Real}
    )
Create a simple orbit and quadratic limb darkening model for 
transit light curve analysis.
# Arguments
- `orbitPeriod`: The orbital period of the planet.
- `transitDuration`: The duration of the transit.
- `limbQuadraticDarkeningParameters`: A vector of coefficients for 
the quadratic limb darkening model.
# Returns
- A tuple containing a `SimpleOrbit` instance and a `QuadLimbDark` 
instance.
"""
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    orbitPeriod::Real, 
    transitDuration::Real, 
    limbQuadraticDarkeningParameters::AbstractVector{<:Real}
) 
    return SimpleOrbit(
            period = orbitPeriod, 
            duration = transitDuration
        ), 
        QuadLimbDark(limbQuadraticDarkeningParameters);
end
function libYukiAstronomyTransitLimbDarkeningQuadratic(
    orbitPeriod::Real, 
    transitDuration::Real, 
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
        transitDuration,
        [u1, u2]
    );
end
