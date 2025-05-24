#include <cmath>
#include <cstring>
#include <mutex>

#include "FFTSpectrum.h"

static std::mutex fftwf_plan_mutex;

#ifndef STATIC_FFTW
template<typename T>
T load_symbol_portable(
#ifdef _WIN32
    HINSTANCE lib_handle,
#else
    void* lib_handle,
#endif
    const char* symbol_name)
{
#ifdef _WIN32
    return reinterpret_cast<T>(GetProcAddress(lib_handle, symbol_name));
#else
    dlerror();
    T func_ptr = reinterpret_cast<T>(dlsym(lib_handle, symbol_name));
    return func_ptr;
#endif
}
#endif

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
            if (buf < 0)
                buf = 0;
            if (buf > 255)
                buf = 255;

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

FFTSpectrum::FFTSpectrum(PClip _child, bool grid, int opt, IScriptEnvironment* env)
    : GenericVideoFilter(_child),
      m_grid(grid)
#ifndef STATIC_FFTW
      ,fftw3_lib_handle(nullptr),
      fftwf_plan_dft_2d(nullptr),
      fftwf_destroy_plan(nullptr),
      fftwf_execute_dft(nullptr)
#endif
{
    if (vi.BitsPerComponent() != 8 || vi.IsRGB() || !vi.IsPlanar())
        env->ThrowError("FFTSpectrum: clip must be in YUV 8-bit planar format.");

#ifndef STATIC_FFTW
#ifdef _WIN32
    const char* fftw_lib_names[] = {"libfftw3f-3.dll", "fftw3.dll"};
#elif defined(__APPLE__)
    const char* fftw_lib_names[] = {"libfftw3f.3.dylib", "libfftw3f.dylib", "libfftw3f.so.3", "libfftw3f.so"};
#else
    const char* fftw_lib_names[] = {"libfftw3f.so.3", "libfftw3f.so"};
#endif

    for (const char* name : fftw_lib_names)
    {
#ifdef _WIN32
        fftw3_lib_handle = LoadLibraryA(name);
#else
        fftw3_lib_handle = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
#endif
        if (fftw3_lib_handle)
            break;
    }

    if (!fftw3_lib_handle)
    {
        std::string attempted_names;
        for (size_t i = 0; i < (sizeof(fftw_lib_names) / sizeof(fftw_lib_names[0])); ++i)
        {
            attempted_names += fftw_lib_names[i];
            if (i < (sizeof(fftw_lib_names) / sizeof(fftw_lib_names[0])) - 1)
                attempted_names += ", ";
        }
#ifndef _WIN32
        const char* error_str = dlerror();
        env->ThrowError(
            "FFTSpectrum: unable to load FFTW3 library (tried: %s). Error: %s", attempted_names.c_str(), error_str ? error_str : "unknown");
#else
        env->ThrowError("FFTSpectrum: unable to load FFTW3 library (tried: %s). Error code: %lu", attempted_names.c_str(), GetLastError());
#endif
    }

    fftwf_plan_dft_2d = load_symbol_portable<fftwf_plan_dft_2d_type>(fftw3_lib_handle, "fftwf_plan_dft_2d");
    fftwf_destroy_plan = load_symbol_portable<fftwf_destroy_plan_type>(fftw3_lib_handle, "fftwf_destroy_plan");
    fftwf_execute_dft = load_symbol_portable<fftwf_execute_dft_type>(fftw3_lib_handle, "fftwf_execute_dft");

    if (!fftwf_plan_dft_2d || !fftwf_destroy_plan || !fftwf_execute_dft)
    {
        if (fftw3_lib_handle)
        {
#ifdef _WIN32
            FreeLibrary(fftw3_lib_handle);
#else
            dlclose(fftw3_lib_handle);
#endif
            fftw3_lib_handle = nullptr;
        }
        env->ThrowError("FFTSpectrum: unable to find required functions in FFTW3 library.");
    }
#endif

    const int64_t plane_size = vi.width * vi.height * sizeof(complex_float);

    fft_in = make_unique_aligned_array_fp<complex_float>(plane_size, 32);
    fft_out = make_unique_aligned_array_fp<complex_float>(plane_size, 32);
    abs_array = make_unique_aligned_array_fp<float>(vi.width * vi.height * sizeof(float), 32);

    if (!fft_in)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_in).");

    if (!fft_out)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_out).");

    if (!abs_array)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (abs_array).");

    {
        const std::lock_guard<std::mutex> lock(fftwf_plan_mutex);
        p = fftwf_plan_dft_2d(vi.height, vi.width, reinterpret_cast<fftwf_complex*>(fft_in.get()),
            reinterpret_cast<fftwf_complex*>(fft_out.get()), FFTW_FORWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    }

    if (vi.NumComponents() > 1)
        vi.pixel_type = VideoInfo::CS_Y8;

    if (opt < -1 || opt > 3)
        env->ThrowError("FFTSpectrum: opt must be between -1..3.");

    const bool avx512 = !!(env->GetCPUFlags() & CPUF_AVX512F);
    const bool avx2 = !!(env->GetCPUFlags() & CPUF_AVX2);
    const bool sse2 = !!(env->GetCPUFlags() & CPUF_SSE2);

    if (!avx512 && opt == 3)
        env->ThrowError("FFTSpectrum: opt=3 requires AVX512.");

    if (!avx2 && opt == 2)
        env->ThrowError("FFTSpectrum: opt=2 requires AVX2.");

    if (!sse2 && opt == 1)
        env->ThrowError("FFTSpectrum: opt=1 requires SSE2.");

    if ((avx512 && opt < 0) || opt == 3)
    {
        fill_fft_input_array = fill_fft_input_array_avx512;
        calculate_absolute_values = calculate_absolute_values_avx512;
    }
    else if ((avx2 && opt < 0) || opt == 2)
    {
        fill_fft_input_array = fill_fft_input_array_avx2;
        calculate_absolute_values = calculate_absolute_values_avx2;
    }
    else if ((sse2 && opt < 0) || opt == 1)
    {
        fill_fft_input_array = fill_fft_input_array_sse2;
        calculate_absolute_values = calculate_absolute_values_sse2;
    }
    else
    {
        fill_fft_input_array = fill_fft_input_array_c;
        calculate_absolute_values = calculate_absolute_values_c;
    }

    has_at_least_v8 = env->FunctionExists("propShow");
}

FFTSpectrum::~FFTSpectrum()
{
    {
        const std::lock_guard<std::mutex> lock(fftwf_plan_mutex);
        fftwf_destroy_plan(p);
    }

#ifndef STATIC_FFTW
    if (fftw3_lib_handle)
    {
#ifdef _WIN32
        FreeLibrary(fftw3_lib_handle);
#else
        dlclose(fftw3_lib_handle);
#endif
        fftw3_lib_handle = nullptr;
    }
#endif
}

PVideoFrame __stdcall FFTSpectrum::GetFrame(int n, IScriptEnvironment* env)
{

    PVideoFrame src = child->GetFrame(n, env);
    const int width = src->GetRowSize();
    const int height = src->GetHeight();

    fill_fft_input_array(fft_in.get(), src->GetReadPtr(), width, height, src->GetPitch());

    fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(fft_in.get()), reinterpret_cast<fftwf_complex*>(fft_out.get()));

    calculate_absolute_values(abs_array.get(), fft_out.get(), (width * height));

    PVideoFrame dst = has_at_least_v8 ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi);

    draw_fft_spectrum(dst->GetWritePtr(), abs_array.get(), width, height, dst->GetPitch());

    if (m_grid)
        draw_grid(dst->GetWritePtr(), width, height, dst->GetPitch());

    return dst;
}

AVSValue __cdecl Create_FFTSpectrum(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new FFTSpectrum(args[0].AsClip(), args[1].AsBool(false), args[2].AsInt(1), env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("FFTSpectrum", "c[grid]b[opt]i", Create_FFTSpectrum, 0);
    return "FFTSpectrum";
}
