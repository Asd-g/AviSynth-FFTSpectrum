#pragma once

#include <cstddef>
#include <memory>

#include <avisynth.h>
#include <fftw3.h>

#include "complex_type.h"

#ifndef STATIC_FFTW
#ifdef _WIN32
#ifndef WIN32_LEAND_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif // !WIN32_LEAND_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif // !NOMINMAX
#include <windows.h>
#else
#include <dlfcn.h>
#endif
typedef fftwf_plan (*fftwf_plan_dft_2d_type)(int n0, int n1, fftwf_complex* in, fftwf_complex* out, int sign, unsigned flags);
typedef void (*fftwf_destroy_plan_type)(fftwf_plan);
typedef void (*fftwf_execute_dft_type)(fftwf_plan, fftwf_complex*, fftwf_complex*);
#endif // STATIC_FFTW

AVS_FORCEINLINE void* aligned_malloc(size_t size, size_t align)
{
    void* result = [&]() {
#ifdef _WIN32
        return _aligned_malloc(size, align);
#else
        if (posix_memalign(&result, align, size))
            return result = nullptr;
        else
            return result;
#endif
    }();

    return result;
}

AVS_FORCEINLINE void aligned_free(void* ptr)
{
#ifdef _WIN32
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

template<typename T>
struct aligned_array_deleter
{
    void operator()(T* ptr) const noexcept
    {
        if (ptr)
            aligned_free(ptr);
    }
};

template<typename T>
using aligned_unique_ptr = std::unique_ptr<T[], aligned_array_deleter<T>>;

template<typename T>
inline aligned_unique_ptr<T> make_unique_aligned_array_fp(size_t num_elements, size_t alignment)
{
    if (num_elements == 0)
        return aligned_unique_ptr<T>(nullptr);

    T* ptr = static_cast<T*>(aligned_malloc(num_elements * sizeof(T), alignment));

    return aligned_unique_ptr<T>(ptr);
}

class FFTSpectrum : public GenericVideoFilter
{
public:
    FFTSpectrum(PClip _child, bool grid, int opt, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override;

    int __stdcall SetCacheHints(int cachehints, int frame_range) override
    {
        return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
    }

    ~FFTSpectrum();

private:
    bool m_grid;

    aligned_unique_ptr<complex_float> fft_in;
    aligned_unique_ptr<complex_float> fft_out;
    fftwf_plan p;
    aligned_unique_ptr<float> abs_array;

    bool has_at_least_v8;

#ifndef STATIC_FFTW
#ifdef _WIN32
    HINSTANCE fftw3_lib_handle;
#else
    void* fftw3_lib_handle;
#endif
    fftwf_plan_dft_2d_type fftwf_plan_dft_2d;
    fftwf_destroy_plan_type fftwf_destroy_plan;
    fftwf_execute_dft_type fftwf_execute_dft;
#endif // !STATIC_FFTW

    void (*fill_fft_input_array)(
        complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept;
    void (*calculate_absolute_values)(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept;
};

void fill_fft_input_array_c(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int src_stride) noexcept;
void calculate_absolute_values_c(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept;
void fill_fft_input_array_sse2(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept;
void calculate_absolute_values_sse2(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept;
void fill_fft_input_array_avx2(complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept;
void calculate_absolute_values_avx2(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept;
void fill_fft_input_array_avx512(
    complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept;
void calculate_absolute_values_avx512(float* __restrict dstp, const complex_float* __restrict srcp, int length) noexcept;
