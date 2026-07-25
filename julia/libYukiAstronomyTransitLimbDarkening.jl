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
    libYukiAstronomyTransitLimbDarkeningSingleTransitBrent(
        lightCurve::libYukiAstronomyTransitLightCurve,
        transit::libYukiAstronomyTransit;
        fitParam::Symbol = :T0,
        searchWindow::Real = 0.05,
        targetDeltaChi2::Real = 1.0,
        errorStepSize::Real = 1e-4
    )
Perform a single transit fit using Brent's method to optimize a specified transit parameter (T0, duration, or depth) based on the provided light curve and transit model. The function returns the best-fit value, associated errors, chi-squared value, and convergence status.
# Arguments
- `lightCurve`: An instance of `libYukiAstronomyTransitLightCurve` containing time and flux data.
- `transit`: An instance of `libYukiAstronomyTransit containing transit parameters.
- `fitParam`: A symbol indicating which parameter to fit (:T0, :duration, or :depth).
- `searchWindow`: A real number specifying the search window around the initial parameter value for optimization.
- `targetDeltaChi2`: A real number specifying the target delta chi-squared value for error estimation.
- `errorStepSize`: A real number specifying the step size for error estimation.
# Returns
- A named tuple containing:
  - `param`: The fitted parameter symbol.
  - `value`: The best-fit value of the parameter.
  - `lowerErr`: The lower error bound for the fitted parameter.
  - `upperErr`: The upper error bound for the fitted parameter.
  - `chi2`: The minimum chi-squared value achieved during optimization.
  - `converged`: A boolean indicating whether the optimization converged successfully.
"""
function libYukiAstronomyTransitLimbDarkeningSingleTransitBrent(
    lightCurve::libYukiAstronomyTransitLightCurve,
    transit::libYukiAstronomyTransit;
    fitParam::Symbol = :T0,
    searchWindow::Real = 0.05,
    targetDeltaChi2::Real = 1.0,
    errorStepSize::Real = 1e-4
)

    objFunction = function(x)
        t0 = (fitParam == :T0) ? 
            x : transit.planetTransitCentreTime;
        duration = (fitParam == :duration) ? 
            x : transit.planetTransitDuration.value;
        depth = (fitParam == :depth) ? 
            x : transit.planetTransitDepth.value;

        timeShifted = lightCurve.time .- (t0 - transit.planetTransitCentreTime);

        try
            lightCurveModel = 
                libYukiAstronomyTransitLimbDarkeningTransitFlux(
                    timeShifted, 
                    depth, 
                    SimpleOrbit(
                        period = transit.planet.orbit.period, 
                        duration = duration
                    ), 
                    transit.limbDarkeningFunc
                );
            
            chi2 = sum(
                (lightCurve.flux .- lightCurveModel.flux) .^ 2 ./ 
                lightCurve.fluxErr .^ 2
            );
            return chi2;
        catch
            return Inf;
        end
    end

    lowerBound, upperBound, initValue = if fitParam == :T0
        (
            transit.planetTransitCentreTime - searchWindow, 
            transit.planetTransitCentreTime + searchWindow, 
            transit.planetTransitCentreTime
        )
    elseif fitParam == :duration
        (
            max(
                1e-4, 
                transit.planetTransitDuration.value * 0.5
            ), 
            transit.planetTransitDuration.value * 1.5, 
            transit.planetTransitDuration.value
        )
    elseif fitParam == :depth
        (
            max(
                1e-5, 
                transit.planetTransitDepth.value * 0.2
            ), 
            min(
                1.0, 
                transit.planetTransitDepth.value * 2.0
            ), 
            transit.planetTransitDepth.value
        )
    else
        error(
            "Unsupported fitParam: $fitParam. " * 
            "Choose from :T0, :duration, :depth"
        )
    end

    result = optimize(objFunction, lowerBound, upperBound, Brent());
    bestValue = Optim.minimizer(result);
    minChi2 = Optim.minimum(result);
    converged = Optim.converged(result);

    val = bestValue
    while (objFunction(val) - minChi2) < targetDeltaChi2
        val += errorStepSize;
    end
    upperErr = val - bestValue;

    val = bestValue
    while (objFunction(val) - minChi2) < targetDeltaChi2
        val -= errorStepSize
    end
    lowerErr = bestValue - val;

    return (
        param = fitParam,
        value = bestValue,
        lowerErr = lowerErr,
        upperErr = upperErr,
        chi2 = minChi2,
        converged = converged
    );
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
        transit.planetTransitDuration.value,
        transit.planetTransitDepth.value,
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
    totalErr = sqrt.(fluxErr .^ 2 .+ exp(logJitter) ^ 2);

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
    transit.planetTransitDepth.value,
    SimpleOrbit(
        period = transit.planet.orbit.period,
        duration = transit.planetTransitDuration.value
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
            transit.planetTransitDuration.value,
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
            transit.planetTransitDuration.value,
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
            transit.planetTransitDuration.value,
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
