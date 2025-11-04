using Measurements

# Derive a Float64 value with uncertainty.
# Dependency: Measurements.
function libYukiBasicMeasurementWithMissing(value::Union{Float64, Missing}, valueErr1::Union{Float64, Missing}, valueErr2::Union{Float64, Missing})::Union{Measurement{Float64}, Missing}
    if ismissing(valueErr1) || ismissing(valueErr2)
        return value ± 0.;
    else
        return value ± ((valueErr1 + valueErr2) / 2);
    end
end
