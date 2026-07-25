using Measurements

"""
	libYukiBasicMeasuredValue(;
		value::Union{Real, Missing} = missing,
		upperErr::Union{Real, Missing} = missing,
		lowerErr::Union{Real, Missing} = missing
	)
Define a measured value with optional upper and lower error values.
# Arguments
- `value`: The measured value, which can be a real number or `missing`.
- `upperErr`: The upper error value, which can be a real number or `missing`.
- `lowerErr`: The lower error value, which can be a real number or `missing	`.
# Returns
- A `libYukiBasicMeasuredValue` object containing the measured value and
its associated error values.
"""
mutable struct libYukiBasicMeasuredValue
	value::Union{Real, Missing}
	upperErr::Union{Real, Missing}
	lowerErr::Union{Real, Missing}
	libYukiBasicMeasuredValue(;
		value::Union{Real, Missing} = missing,
		upperErr::Union{Real, Missing} = missing,
		lowerErr::Union{Real, Missing} = missing
	) = new(value, upperErr, lowerErr)
end

"""
	libYukiBasicMeasurementWithMissing(
		value::Union{Real, Missing},
		valueErr1::Union{Real, Missing},
		valueErr2::Union{Real, Missing}
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
libYukiBasicMeasurementWithMissing(
	measured::libYukiBasicMeasuredValue
) = libYukiBasicMeasurementWithMissing(
	measured.value, 
	measured.upperErr, 
	measured.lowerErr
);
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
