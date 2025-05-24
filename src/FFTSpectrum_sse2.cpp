#include "complex_type.h"
#include "FFTSpectrum.h"
#include "vcl_log_constants.h"
#include "vcl_utils.h"

void fill_fft_input_array_sse2(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept
{
    vcl_utils::fill_fft_input_array_templated<Vec4f, complex_float>(dstp, srcp, width, height, stride);
}

void calculate_absolute_values_sse2(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept
{
    vcl_utils::calculate_absolute_values_templated<Vec4f, false>(dstp, srcp, length);
}
