using Measurements

"""
	libYukiBasicMeasurementWithMissing(
		value, 
		valueErr1, 
		valueErr2
	)
Calculate a measurement with missing error values.
# Arguments
- `value`: The measured value, which can be a real number or `missing`.
- `valueErr1`: The first error value, which can be a real number or 
`missing`.
- `valueErr2`: The second error value, which can be a real number or 
`missing`.
# Returns
- A `Measurement` object if both error values are provided, otherwise 
a measurement with zero error.
"""
function libYukiBasicMeasurementWithMissing(
	value::Union{Real, Missing}, 
	valueErr1::Union{Real, Missing}, 
	valueErr2::Union{Real, Missing}
)::Union{Measurement, Missing}
	if ismissing(valueErr1) || ismissing(valueErr2)
		return value ± 0.;
	else
		return value ± ((valueErr1 + valueErr2) / 2);
	end
end
