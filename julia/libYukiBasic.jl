using Measurements

# Sort elements by order vector.
function libYukiBasicSortElementsByOrderVector(orderVector, elementsVector)
	if length(orderVector) != length(elementsVector)
		error("ERR: libYuki Basic Sort Elements By Order Vector meets 2 different length vectors.")
	else
		sortedIndices = sortperm(orderVector);
		return orderVector[sortedIndices], elementsVector[sortedIndices];
	end
end

# Derive a Float64 value with uncertainty.
# Dependency: Measurements.
function libYukiBasicMeasurementWithMissing(value::Union{Real, Missing}, valueErr1::Union{Real, Missing}, valueErr2::Union{Real, Missing})::Union{Measurement, Missing}
	if ismissing(valueErr1) || ismissing(valueErr2)
		return value ± 0.;
	else
		return value ± ((valueErr1 + valueErr2) / 2);
	end
end
