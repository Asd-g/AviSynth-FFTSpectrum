#include <algorithm>
#include <cmath>
#include <cstring>
#include <mutex>
#include <span>

#include "FFTSpectrum.h"

inline std::mutex& GetFFTWMutex()
{
    static std::mutex m;
    return m;
}

class GlobalLockGuard
{
    IScriptEnvironment* env;
    const char* name;
    bool acquired;
    bool is_legacy;

public:
    GlobalLockGuard(IScriptEnvironment* _env, const char* _name, bool use_avs_lock)
        : env(_env), name(_name), acquired(false), is_legacy(false)
    {
        if (!name)
            return;
        if (use_avs_lock)
        {
            acquired = env->AcquireGlobalLock(name);
            if (acquired)
                return;
            else
                env->ThrowError("FFTSpectrum: env is corrupted or the lock name is invalid.");
        }
        if (strcmp(name, "fftw") == 0)
        {
            GetFFTWMutex().lock();
            acquired = true;
            is_legacy = true;
        }
    }

    ~GlobalLockGuard()
    {
        if (acquired)
        {
            if (is_legacy)
                GetFFTWMutex().unlock();
            else
                env->ReleaseGlobalLock(name);
        }
    }

    GlobalLockGuard(const GlobalLockGuard&) = delete;
    GlobalLockGuard& operator=(const GlobalLockGuard&) = delete;
};

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
    std::fill_n(dstp, static_cast<size_t>(stride) * height, 0);

    const std::span<const float> src_span(srcp, static_cast<size_t>(height) * width);
    const float max = *std::ranges::max_element(src_span);

    if (max <= 0.0f)
        return;

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
      ,
      fftw3_lib_handle(nullptr),
      fftwf_plan_dft_2d(nullptr),
      fftwf_destroy_plan(nullptr),
      fftwf_execute_dft(nullptr),
      has_at_least_v8(env->FunctionExists("propShow")),
      has_at_least_v12(true),
      m_env(env)
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

    if (opt < -1 || opt > 3)
        env->ThrowError("FFTSpectrum: opt must be between -1..3.");

    const bool avx512 = !!(env->GetCPUFlags() & CPUF_AVX512F) && (opt < 0 || opt == 3);
    const bool avx2 = !!(env->GetCPUFlags() & CPUF_AVX2) && (opt < 0 || opt == 2);
    const bool sse2 = !!(env->GetCPUFlags() & CPUF_SSE2) && (opt < 0 || opt == 1);

    if (!avx512 && opt == 3)
        env->ThrowError("FFTSpectrum: opt=3 requires AVX512.");

    if (!avx2 && opt == 2)
        env->ThrowError("FFTSpectrum: opt=2 requires AVX2.");

    if (!sse2 && opt == 1)
        env->ThrowError("FFTSpectrum: opt=1 requires SSE2.");

    if (avx512)
    {
        fill_fft_input_array = fill_fft_input_array_avx512;
        calculate_absolute_values = calculate_absolute_values_avx512;
    }
    else if (avx2)
    {
        fill_fft_input_array = fill_fft_input_array_avx2;
        calculate_absolute_values = calculate_absolute_values_avx2;
    }
    else if (sse2)
    {
        fill_fft_input_array = fill_fft_input_array_sse2;
        calculate_absolute_values = calculate_absolute_values_sse2;
    }
    else
    {
        fill_fft_input_array = fill_fft_input_array_c;
        calculate_absolute_values = calculate_absolute_values_c;
    }

    const int alignment = (avx512) ? 64 : 32;

    fft_in = make_unique_aligned_array_fp<complex_float>(plane_size, alignment);
    fft_out = make_unique_aligned_array_fp<complex_float>(plane_size, alignment);
    abs_array = make_unique_aligned_array_fp<float>(vi.width * vi.height * sizeof(float), alignment);

    if (!fft_in)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_in).");

    if (!fft_out)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (fft_out).");

    if (!abs_array)
        env->ThrowError("FFTSpectrum: _aligned_malloc failure (abs_array).");

    try
    {
        env->CheckVersion(12);
    }
    catch (const AvisynthError&)
    {
        has_at_least_v12 = false;
    }

    {
        GlobalLockGuard fftw_lock(env, "fftw", has_at_least_v12);
        p = fftwf_plan_dft_2d(vi.height, vi.width, reinterpret_cast<fftwf_complex*>(fft_in.get()),
            reinterpret_cast<fftwf_complex*>(fft_out.get()), FFTW_FORWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    }

    if (!p)
        env->ThrowError("FFTSpectrum: FFTW failed to create a plan.");

    if (vi.NumComponents() > 1)
        vi.pixel_type = VideoInfo::CS_Y8;
}

FFTSpectrum::~FFTSpectrum()
{
    {
        GlobalLockGuard fftw_lock(m_env, "fftw", has_at_least_v12);
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
