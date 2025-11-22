using Measurements

# Sort 2 dimension vector by the 1st demension.
# Dependency: Measurements.
function libYukiBasicSort2VectorsByTheFirstOne(vec1::Vector{Measurement{Float64}}, vec2::Vector{Measurement{Float64}})::Tuple{Vector{Measurement{Float64}}, Vector{Measurement{Float64}}}
    if length(vec1) != length(vec2)
        error("ERR: libYuki Basic Sort 2 Vectors By The First One meets 2 different length vectors.")
    else
        sortedIndices::Vector{Int64} = sortperm(vec1);
        return vec1[sortedIndices], vec2[sortedIndices];
    end
end

# Derive a Float64 value with uncertainty.
# Dependency: Measurements.
function libYukiBasicMeasurementWithMissing(value::Union{Float64, Missing}, valueErr1::Union{Float64, Missing}, valueErr2::Union{Float64, Missing})::Union{Measurement{Float64}, Missing}
    if ismissing(valueErr1) || ismissing(valueErr2)
        return value ± 0.;
    else
        return value ± ((valueErr1 + valueErr2) / 2);
    end
end
