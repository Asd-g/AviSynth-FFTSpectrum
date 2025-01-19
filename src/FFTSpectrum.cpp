#include <cstring>
#include <mutex>
#include <Windows.h>

#include "fftw3.h"
#include "avisynth.h"

#if defined(_MSC_VER)
#include <intrin.h>

#define USE_SSE_AUTO
#define __SSE4_2__
#define __x86_64__
#define SSE_MATHFUN_WITH_CODE
#include "sse_mathfun.h"
#undef SSE_MATHFUN_WITH_CODE
#undef __x86_64__
#undef __SSE4_2__
#undef USE_SSE_AUTO
#undef inline

#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>

#define USE_SSE4
#define SSE_MATHFUN_WITH_CODE
#include "sse_mathfun.h"

#endif

#ifdef LOADPLUGIN_DLL
typedef float fftwf_complex[2];
typedef fftwf_plan(*fftwf_plan_dft_2d_)(int n0, int n1, fftwf_complex* in, fftwf_complex* out, int sign, unsigned flags);
typedef void (*fftwf_destroy_plan_)(fftwf_plan);
typedef void (*fftwf_execute_dft_)(fftwf_plan, fftwf_complex*, fftwf_complex*);
#endif

static std::mutex fftwf_;

static void fill_fft_input_array(fftwf_complex* dstp, const uint8_t* srcp, int width, int height, int stride)
{

    const int mod16_width = width - (width % 16);
    __m128 sse_zero = _mm_setzero_ps();

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < mod16_width; x += 16)
        {
            __m128i in_buffer = _mm_load_si128((const __m128i*)srcp);

            for (int j = 0; j < 4; ++j)
            {
                __m128i epu32_buffer = _mm_cvtepu8_epi32(in_buffer);
                __m128 cvt_buffer = _mm_cvtepi32_ps(epu32_buffer);

                __m128 out_buffer = _mm_unpacklo_ps(cvt_buffer, sse_zero);
                __m128 out_buffer1 = _mm_unpackhi_ps(cvt_buffer, sse_zero);

                _mm_store_ps(reinterpret_cast<float*>(dstp), out_buffer);
                _mm_store_ps(reinterpret_cast<float*>(dstp + 2), out_buffer1);

                in_buffer = _mm_shuffle_epi32(in_buffer, _MM_SHUFFLE(0, 3, 2, 1));

                dstp += 4;
            }

            srcp += 16;
        }
        for (int x = mod16_width; x < width; ++x)
        {
            *dstp[0] = static_cast<float>(*srcp);
            *dstp[1] = 0.0;
            ++srcp;
            ++dstp;
        }

        srcp += static_cast<int64_t>(stride) - width;
    }
}

static void calculate_absolute_values(float* dstp, fftwf_complex* srcp, int length)
{
    const int mod4_length = length - (length % 4);
    __m128 sse_one = _mm_set_ps1(1.0f);

    for (int i = 0; i < mod4_length; i += 4)
    {
        __m128 in_buffer = _mm_load_ps(reinterpret_cast<float*>(srcp));
        __m128 in_buffer1 = _mm_load_ps(reinterpret_cast<float*>(srcp + 2));

        __m128 mul_buffer = _mm_mul_ps(in_buffer, in_buffer);
        __m128 mul_buffer1 = _mm_mul_ps(in_buffer1, in_buffer1);

        __m128 add_buffer = _mm_hadd_ps(mul_buffer, mul_buffer1);
        add_buffer = _mm_sqrt_ps(add_buffer);
        add_buffer = _mm_add_ps(add_buffer, sse_one);

        __m128 out_buffer = log_ps(add_buffer);

        _mm_store_ps(dstp, out_buffer);

        srcp += 4;
        dstp += 4;
    }

    for (int i = mod4_length; i < length; ++i)
        dstp[i] = logf(sqrtf(srcp[i][0] * srcp[i][0] + srcp[i][1] * srcp[i][1]) + 1.0);
}

static void draw_fft_spectrum(uint8_t* dstp, float* srcp, int width, int height, int stride)
{
    float max = 0.f;

    memset(dstp, 0, static_cast<int64_t>(stride) * height);

    for (int i = 1; i < height * width; ++i)
    {
        if (srcp[i] > max)
            max = srcp[i];
    }

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float buf = srcp[x + y * width] > max / 2 ? srcp[x + y * width] : 0;
            buf = 255 * buf / max;
            if (buf < 0) buf = 0;
            if (buf > 255) buf = 255;

            if (y < height / 2)
            {
                if (x < width / 2)
                    dstp[x + (width / 2) + stride * (y + height / 2)] = static_cast<uint8_t>(lrintf(buf));
                else
                    dstp[x - (width / 2) + stride * (y + height / 2)] = static_cast<uint8_t>(lrintf(buf));
            }
            else
            {
                if (x < width / 2)
                    dstp[x + (width / 2) + stride * (y - height / 2)] = static_cast<uint8_t>(lrintf(buf));
                else
                    dstp[x - (width / 2) + stride * (y - height / 2)] = static_cast<uint8_t>(lrintf(buf));
            }
        }
    }
}

static void draw_grid(uint8_t* buf, int width, int height, int stride)
{
    for (int x = (width / 2) % 100; x < width; x += 100)
    {
        for (int y = 0; y < height; ++y)
            buf[x + y * stride] = 255;
    }

    for (int y = (height / 2) % 100; y < height; y += 100)
    {
        for (int x = 0; x < width; ++x)
            buf[x + y * stride] = 255;
    }
}

class FFTSpectrum : public GenericVideoFilter
{
    bool _grid;

    fftwf_complex* fft_in;
    fftwf_complex* fft_out;
    fftwf_plan p;
    float* abs_array;

    bool has_at_least_v8;
    HINSTANCE fftw3_lib;
    fftwf_plan_dft_2d_ fftwf_plan_dft_2d;
    fftwf_destroy_plan_ fftwf_destroy_plan;
    fftwf_execute_dft_ fftwf_execute_dft;

public:
    FFTSpectrum(PClip _child, bool grid, IScriptEnvironment* env)
        : GenericVideoFilter(_child), _grid(grid)
    {
        if (vi.BitsPerComponent() != 8 || vi.IsRGB() || !vi.IsPlanar())
            env->ThrowError("FFTSpectrum: clip must be in YUV 8-bit planar format.");

#ifdef LOADPLUGIN_DLL
        fftw3_lib = LoadLibrary("libfftw3f-3.dll");
        if (!fftw3_lib)
            fftw3_lib = LoadLibrary("fftw3.dll");

        if (fftw3_lib)
        {
            fftwf_plan_dft_2d = reinterpret_cast<fftwf_plan_dft_2d_>(GetProcAddress(fftw3_lib, "fftwf_plan_dft_2d"));
            fftwf_destroy_plan = reinterpret_cast<fftwf_destroy_plan_>(GetProcAddress(fftw3_lib, "fftwf_destroy_plan"));
            fftwf_execute_dft = reinterpret_cast<fftwf_execute_dft_>(GetProcAddress(fftw3_lib, "fftwf_execute_dft"));
        }

        if (!fftw3_lib || !fftwf_plan_dft_2d || !fftwf_destroy_plan)
            env->ThrowError("FFTSpectrum: unable to load libfftw3f-3.dll or fftw3.dll.");
#endif

        const int64_t width = vi.width;
        const int height = vi.height;
        const int fft_size = sizeof(fftw_complex);
        const int float_size = sizeof(float);

        fft_in = static_cast<fftwf_complex*>(_aligned_malloc((width * height * fft_size), 32));
        if (fft_in)
            memset(fft_in, 0, (width * height * fft_size));
        else
            env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_in).");

        fft_out = static_cast<fftwf_complex*>(_aligned_malloc((width * height * fft_size), 32));
        if (fft_out)
            memset(fft_out, 0, (width * height * fft_size));
        else
            env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_out).");

        abs_array = static_cast<float*>(_aligned_malloc((width * height * float_size), 32));
        if (abs_array)
            memset(abs_array, 0, (width * height * float_size));
        else
            env->ThrowError("FFTSpectrum: _aligned_malloc failure (abs_array).");

        const std::lock_guard<std::mutex> lock(fftwf_);
        p = fftwf_plan_dft_2d(height, width, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);

        if (vi.NumComponents() > 1)
            vi.pixel_type = VideoInfo::CS_Y8;

        has_at_least_v8 = true;
        try { env->CheckVersion(8); }
        catch (const AvisynthError&)
        { has_at_least_v8 = false; }
    }

    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
    }

    ~FFTSpectrum()
    {
        _aligned_free(fft_in);
        _aligned_free(fft_out);
        _aligned_free(abs_array);
        fftwf_destroy_plan(p);
    }

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env)
    {

        PVideoFrame src = child->GetFrame(n, env);
        const int width = src->GetRowSize();
        const int height = src->GetHeight();

        fill_fft_input_array(fft_in, src->GetReadPtr(), width, height, src->GetPitch());

        fftwf_execute_dft(p, fft_in, fft_out);

        calculate_absolute_values(abs_array, fft_out, (width * height));

        PVideoFrame dst = has_at_least_v8 ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi);

        draw_fft_spectrum(dst->GetWritePtr(), abs_array, width, height, dst->GetPitch());

        if (_grid)
            draw_grid(dst->GetWritePtr(), width, height, dst->GetPitch());

        return dst;
    }
};

AVSValue __cdecl Create_FFTSpectrum(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new FFTSpectrum(
        args[0].AsClip(),
        args[1].AsBool(false),
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("FFTSpectrum", "c[grid]b", Create_FFTSpectrum, 0);
    return "FFTSpectrum";
}
