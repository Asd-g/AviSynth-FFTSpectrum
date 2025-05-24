#include <bit>
#include <cmath>

#include "FFTSpectrum.h"

constexpr float ONE = 1.0f;
constexpr float P0_5 = 0.5f;

constexpr uint32_t BITS_MIN_NORM_POS = 0x00800000u;
constexpr uint32_t BITS_INV_MANT_MASK = 0x807FFFFFu;
constexpr uint32_t BITS_0P5 = 0x3F000000u;
constexpr int32_t FLOAT_EXP_BIAS = 127;

constexpr float CEPHES_SQRTHF = 0.707106781186547524f;
constexpr float CEPHES_LOG_P0 = 7.0376836292E-2f;
constexpr float CEPHES_LOG_P1 = -1.1514610310E-1f;
constexpr float CEPHES_LOG_P2 = 1.1676998740E-1f;
constexpr float CEPHES_LOG_P3 = -1.2420140846E-1f;
constexpr float CEPHES_LOG_P4 = 1.4249322787E-1f;
constexpr float CEPHES_LOG_P5 = -1.6668057665E-1f;
constexpr float CEPHES_LOG_P6 = 2.0000714765E-1f;
constexpr float CEPHES_LOG_P7 = -2.4999993993E-1f;
constexpr float CEPHES_LOG_P8 = 3.3333331174E-1f;
constexpr float CEPHES_LOG_Q1 = -2.12194440e-4f;
constexpr float CEPHES_LOG_Q2 = 0.693359375f;

AVS_FORCEINLINE float log_ps_c(float x_input)
{
    if (x_input <= 0.0f || std::isnan(x_input))
        return std::numeric_limits<float>::quiet_NaN();

    if (std::isinf(x_input))
        return x_input;

    float x = std::max(x_input, std::bit_cast<float>(BITS_MIN_NORM_POS));

    uint32_t ix = std::bit_cast<uint32_t>(x);

    int32_t exp_val = static_cast<int32_t>((ix >> 23) & 0xFFu);

    ix = (ix & BITS_INV_MANT_MASK) | BITS_0P5;
    x = std::bit_cast<float>(ix);

    exp_val -= FLOAT_EXP_BIAS;
    float e_float = static_cast<float>(exp_val);
    e_float += ONE;

    if (x < CEPHES_SQRTHF)
    {
        e_float -= ONE;
        x = x + x - ONE;
    }
    else
        x = x - ONE;

    const float z = x * x;

    float poly_y = CEPHES_LOG_P0;
    poly_y = poly_y * x + CEPHES_LOG_P1;
    poly_y = poly_y * x + CEPHES_LOG_P2;
    poly_y = poly_y * x + CEPHES_LOG_P3;
    poly_y = poly_y * x + CEPHES_LOG_P4;
    poly_y = poly_y * x + CEPHES_LOG_P5;
    poly_y = poly_y * x + CEPHES_LOG_P6;
    poly_y = poly_y * x + CEPHES_LOG_P7;
    poly_y = poly_y * x + CEPHES_LOG_P8;
    poly_y = poly_y * x;

    poly_y = poly_y * z;

    float tmp = e_float * CEPHES_LOG_Q1;
    poly_y = poly_y + tmp;

    tmp = z * P0_5;
    poly_y = poly_y - tmp;

    x = x + poly_y;
    tmp = e_float * CEPHES_LOG_Q2;
    x = x + tmp;

    return x;
}

void fill_fft_input_array_c(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int src_stride) noexcept
{
    for (int y = 0; y < height; ++y)
    {
        const uint8_t* p_src = srcp + static_cast<ptrdiff_t>(y) * src_stride;
        complex_float* p_dst = dstp + static_cast<ptrdiff_t>(y) * width;
        for (int x = 0; x < width; ++x)
        {
            p_dst[x].re = static_cast<float>(p_src[x]);
            p_dst[x].im = 0.0f;
        }
    }
}

void calculate_absolute_values_c(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept
{
    for (int i = 0; i < length; ++i)
    {
        const float re = srcp[i].re;
        const float im = srcp[i].im;

        const float mag_sq = re * re + im * im;
        const float mag = sqrtf(mag_sq);

        const float val_to_log = mag + ONE;

        dstp[i] = log_ps_c(val_to_log);
        // Alternatively, for standard library log:
        // if (val_to_log <= 0.0f) { // logf domain check
        //     dstp_abs[i] = std::numeric_limits<float>::quiet_NaN();
        // } else {
        //     dstp_abs[i] = logf(val_to_log); // from <cmath>
        // }
    }
}
