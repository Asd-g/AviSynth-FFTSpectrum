#include "FFTSpectrum.h"
#include "complex_type.h"
#include "vcl_log_constants.h"
#include "vcl_utils.h"

void fill_fft_input_array_avx512(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept
{
    vcl_utils::fill_fft_input_array_templated<Vec16f, complex_float>(dstp, srcp, width, height, stride);
}

void calculate_absolute_values_avx512(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept
{
    vcl_utils::calculate_absolute_values_templated<Vec16f, true>(dstp, srcp, length);
}
