#pragma once

#include <cmath>
#include <type_traits>

#include "../VCL2/vectorclass.h"
#include "complex_type.h"
#include <avs/config.h>

namespace vcl_utils
{
    template<typename float_vector_type>
    AVS_FORCEINLINE static float_vector_type load_n_uint8_to_float(const uint8_t* p)
    {
        if constexpr (std::is_same_v<float_vector_type, Vec4f>)
            return to_float(Vec4i().load_4uc(p));
#if INSTRSET >= 8
        else if constexpr (std::is_same_v<float_vector_type, Vec8f>)
            return to_float(Vec8i().load_8uc(p));
#endif
#if INSTRSET >= 10
        else if constexpr (std::is_same_v<float_vector_type, Vec16f>)
            return to_float(Vec16i().load_16uc(p));
#endif
    }

    template<typename float_vector_type, typename complex_float>
    AVS_FORCEINLINE void fill_fft_input_array_templated(
        complex_float* __restrict dstp, const uint8_t* __restrict srcp, int width, int height, int stride) noexcept
    {
        constexpr int uint8_per_native_vector_load = float_vector_type::size();
        constexpr int ops_in_unrolled_loop = 4;
        constexpr int uint8_in_unrolled_loop = ops_in_unrolled_loop * uint8_per_native_vector_load;

        const int mod_width_unrolled = width - (width % uint8_in_unrolled_loop);

        for (int y = 0; y < height; ++y)
        {
            const uint8_t* p_src = srcp + static_cast<ptrdiff_t>(y) * stride;
            complex_float* p_dst_base = dstp + static_cast<ptrdiff_t>(y) * width;

            for (int x = 0; x < mod_width_unrolled; x += uint8_in_unrolled_loop)
            {
                float_vector_type src_parts[ops_in_unrolled_loop] = {vcl_utils::load_n_uint8_to_float<float_vector_type>(p_src + x),
                    vcl_utils::load_n_uint8_to_float<float_vector_type>(p_src + x + uint8_per_native_vector_load),
                    vcl_utils::load_n_uint8_to_float<float_vector_type>(p_src + x + 2 * uint8_per_native_vector_load),
                    vcl_utils::load_n_uint8_to_float<float_vector_type>(p_src + x + 3 * uint8_per_native_vector_load)};

                const int complex_offset_in_dst[ops_in_unrolled_loop] = {
                    x, x + uint8_per_native_vector_load, x + 2 * uint8_per_native_vector_load, x + 3 * uint8_per_native_vector_load};

                if constexpr (std::is_same_v<float_vector_type, Vec4f>)
                {
                    Vec4f out[8] = {permute4<0, -1, 1, -1>(src_parts[0]), permute4<2, -1, 3, -1>(src_parts[0]),
                        permute4<0, -1, 1, -1>(src_parts[1]), permute4<2, -1, 3, -1>(src_parts[1]), permute4<0, -1, 1, -1>(src_parts[2]),
                        permute4<2, -1, 3, -1>(src_parts[2]), permute4<0, -1, 1, -1>(src_parts[3]), permute4<2, -1, 3, -1>(src_parts[3])};

                    out[0].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0]));
                    out[1].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0] + 2));

                    out[2].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1]));
                    out[3].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1] + 2));

                    out[4].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2]));
                    out[5].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2] + 2));

                    out[6].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3]));
                    out[7].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3] + 2));
                }
#if INSTRSET >= 8
                else if constexpr (std::is_same_v<float_vector_type, Vec8f>)
                {
                    Vec8f out[8] = {permute8<0, -1, 1, -1, 2, -1, 3, -1>(src_parts[0]), permute8<4, -1, 5, -1, 6, -1, 7, -1>(src_parts[0]),
                        permute8<0, -1, 1, -1, 2, -1, 3, -1>(src_parts[1]), permute8<4, -1, 5, -1, 6, -1, 7, -1>(src_parts[1]),
                        permute8<0, -1, 1, -1, 2, -1, 3, -1>(src_parts[2]), permute8<4, -1, 5, -1, 6, -1, 7, -1>(src_parts[2]),
                        permute8<0, -1, 1, -1, 2, -1, 3, -1>(src_parts[3]), permute8<4, -1, 5, -1, 6, -1, 7, -1>(src_parts[3])};

                    out[0].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0]));
                    out[1].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0] + 4));

                    out[2].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1]));
                    out[3].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1] + 4));

                    out[4].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2]));
                    out[5].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2] + 4));

                    out[6].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3]));
                    out[7].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3] + 4));
                }
#endif
#if INSTRSET >= 10
                else if constexpr (std::is_same_v<float_vector_type, Vec16f>)
                {
                    Vec16f out[8] = {permute16<0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1>(src_parts[0]),
                        permute16<8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1>(src_parts[0]),
                        permute16<0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1>(src_parts[1]),
                        permute16<8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1>(src_parts[1]),
                        permute16<0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1>(src_parts[2]),
                        permute16<8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1>(src_parts[2]),
                        permute16<0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1>(src_parts[3]),
                        permute16<8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1>(src_parts[3])};

                    out[0].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0]));
                    out[1].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[0] + 8));

                    out[2].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1]));
                    out[3].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[1] + 8));

                    out[4].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2]));
                    out[5].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[2] + 8));

                    out[6].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3]));
                    out[7].store_a(reinterpret_cast<float*>(p_dst_base + complex_offset_in_dst[3] + 8));
                }
#endif
            }

            for (int x = mod_width_unrolled; x < width; ++x)
            {
                p_dst_base[x].re = static_cast<float>(p_src[x]);
                p_dst_base[x].im = 0.0f;
            }
        }
    }

    template<typename float_vector_type>
    AVS_FORCEINLINE static void load_deinterleaved(
        float_vector_type& real_parts, float_vector_type& imag_parts, const float* interleaved_data)
    {
        if constexpr (std::is_same_v<float_vector_type, Vec4f>)
        {
            // Load 8 floats (4 complex numbers)
            Vec4f low_chunk = Vec4f().load(interleaved_data + 0);  // (r0, i0, r1, i1)
            Vec4f high_chunk = Vec4f().load(interleaved_data + 4); // (r2, i2, r3, i3)
            // De-interleave
            // real_parts = (r0, r1, r2, r3)
            // imag_parts = (i0, i1, i2, i3)
            real_parts = blend4<0, 2, 0 + 4, 2 + 4>(low_chunk, high_chunk); // (low[0],low[2],high[0],high[2])
            imag_parts = blend4<1, 3, 1 + 4, 3 + 4>(low_chunk, high_chunk); // (low[1],low[3],high[1],high[3])
        }
#if INSTRSET >= 8
        else if constexpr (std::is_same_v<float_vector_type, Vec8f>)
        {
            // Load 16 floats (8 complex numbers)
            Vec8f chunk0 = Vec8f().load(interleaved_data + 0); // r0,i0,r1,i1,r2,i2,r3,i3
            Vec8f chunk1 = Vec8f().load(interleaved_data + 8); // r4,i4,r5,i5,r6,i6,r7,i7
            // De-interleave
            real_parts = blend8<0, 2, 4, 6, 0 + 8, 2 + 8, 4 + 8, 6 + 8>(chunk0, chunk1);
            imag_parts = blend8<1, 3, 5, 7, 1 + 8, 3 + 8, 5 + 8, 7 + 8>(chunk0, chunk1);
        }
#endif
#if INSTRSET >= 10
        else if constexpr (std::is_same_v<float_vector_type, Vec16f>)
        {
            // Load 32 floats (16 complex numbers)
            Vec16f chunk0 = Vec16f().load(interleaved_data + 0);  // r0..i7
            Vec16f chunk1 = Vec16f().load(interleaved_data + 16); // r8..i15
            // De-interleave
            real_parts =
                blend16<0, 2, 4, 6, 8, 10, 12, 14, 0 + 16, 2 + 16, 4 + 16, 6 + 16, 8 + 16, 10 + 16, 12 + 16, 14 + 16>(chunk0, chunk1);
            imag_parts =
                blend16<1, 3, 5, 7, 9, 11, 13, 15, 1 + 16, 3 + 16, 5 + 16, 7 + 16, 9 + 16, 11 + 16, 13 + 16, 15 + 16>(chunk0, chunk1);
        }
#endif
    }

    template<typename float_vector_type, bool intermediate_vectors>
    AVS_FORCEINLINE void calculate_absolute_values_templated(
        float* __restrict dstp, const complex_float* __restrict src, int length) noexcept
    {
        constexpr int complex_per_native_vector_op = float_vector_type::size();
        constexpr int ops_in_unrolled_loop = 4;
        constexpr int complex_in_unrolled_loop = ops_in_unrolled_loop * complex_per_native_vector_op;

        const int mod_length_unrolled = length - (length % complex_in_unrolled_loop);
        const float_vector_type vcl_one(1.0f);

        const float* srcp_float = reinterpret_cast<const float*>(src);

        for (int i = 0; i < mod_length_unrolled; i += complex_in_unrolled_loop)
        {
            if constexpr (intermediate_vectors)
            {

                float_vector_type real_parts_arr[ops_in_unrolled_loop];
                float_vector_type imag_parts_arr[ops_in_unrolled_loop];

                load_deinterleaved<float_vector_type>(real_parts_arr[0], imag_parts_arr[0], srcp_float + i * 2);
                load_deinterleaved<float_vector_type>(
                    real_parts_arr[1], imag_parts_arr[1], srcp_float + (i + complex_per_native_vector_op) * 2);
                load_deinterleaved<float_vector_type>(
                    real_parts_arr[2], imag_parts_arr[2], srcp_float + (i + 2 * complex_per_native_vector_op) * 2);
                load_deinterleaved<float_vector_type>(
                    real_parts_arr[3], imag_parts_arr[3], srcp_float + (i + 3 * complex_per_native_vector_op) * 2);

                float_vector_type abs_val_sq_arr[ops_in_unrolled_loop] = {
                    mul_add(real_parts_arr[0], real_parts_arr[0], imag_parts_arr[0] * imag_parts_arr[0]),
                    mul_add(real_parts_arr[1], real_parts_arr[1], imag_parts_arr[1] * imag_parts_arr[1]),
                    mul_add(real_parts_arr[2], real_parts_arr[2], imag_parts_arr[2] * imag_parts_arr[2]),
                    mul_add(real_parts_arr[3], real_parts_arr[3], imag_parts_arr[3] * imag_parts_arr[3])};

                float_vector_type log_abs_val_arr[ops_in_unrolled_loop] = {
                    log_ps_vcl_generic<float_vector_type>(sqrt(abs_val_sq_arr[0]) + vcl_one),
                    log_ps_vcl_generic<float_vector_type>(sqrt(abs_val_sq_arr[1]) + vcl_one),
                    log_ps_vcl_generic<float_vector_type>(sqrt(abs_val_sq_arr[2]) + vcl_one),
                    log_ps_vcl_generic<float_vector_type>(sqrt(abs_val_sq_arr[3]) + vcl_one)};

                log_abs_val_arr[0].store_a(dstp + i);
                log_abs_val_arr[1].store_a(dstp + i + complex_per_native_vector_op);
                log_abs_val_arr[2].store_a(dstp + i + 2 * complex_per_native_vector_op);
                log_abs_val_arr[3].store_a(dstp + i + 3 * complex_per_native_vector_op);
            }
            else
            {
                float_vector_type real_parts_arr[ops_in_unrolled_loop];
                float_vector_type imag_parts_arr[ops_in_unrolled_loop];
                float_vector_type abs_val_sq_arr[ops_in_unrolled_loop];
                float_vector_type abs_val_arr[ops_in_unrolled_loop];
                float_vector_type log_abs_val_arr[ops_in_unrolled_loop];

                const float* current_src_float0 = srcp_float + i * 2;
                load_deinterleaved<float_vector_type>(real_parts_arr[0], imag_parts_arr[0], current_src_float0);
                abs_val_sq_arr[0] = mul_add(real_parts_arr[0], real_parts_arr[0], imag_parts_arr[0] * imag_parts_arr[0]);
                abs_val_arr[0] = sqrt(abs_val_sq_arr[0]);
                abs_val_arr[0] = abs_val_arr[0] + vcl_one;
                log_abs_val_arr[0] = log_ps_vcl_generic<float_vector_type>(abs_val_arr[0]);
                log_abs_val_arr[0].store_a(dstp + i);

                const int complex_start_idx1 = i + complex_per_native_vector_op;
                const float* current_src_float1 = srcp_float + complex_start_idx1 * 2;
                load_deinterleaved<float_vector_type>(real_parts_arr[1], imag_parts_arr[1], current_src_float1);
                abs_val_sq_arr[1] = mul_add(real_parts_arr[1], real_parts_arr[1], imag_parts_arr[1] * imag_parts_arr[1]);
                abs_val_arr[1] = sqrt(abs_val_sq_arr[1]);
                abs_val_arr[1] = abs_val_arr[1] + vcl_one;
                log_abs_val_arr[1] = log_ps_vcl_generic<float_vector_type>(abs_val_arr[1]);
                log_abs_val_arr[1].store_a(dstp + complex_start_idx1);

                const int complex_start_idx2 = i + 2 * complex_per_native_vector_op;
                const float* current_src_float2 = srcp_float + complex_start_idx2 * 2;
                load_deinterleaved<float_vector_type>(real_parts_arr[2], imag_parts_arr[2], current_src_float2);
                abs_val_sq_arr[2] = mul_add(real_parts_arr[2], real_parts_arr[2], imag_parts_arr[2] * imag_parts_arr[2]);
                abs_val_arr[2] = sqrt(abs_val_sq_arr[2]);
                abs_val_arr[2] = abs_val_arr[2] + vcl_one;
                log_abs_val_arr[2] = log_ps_vcl_generic<float_vector_type>(abs_val_arr[2]);
                log_abs_val_arr[2].store_a(dstp + complex_start_idx2);

                const int complex_start_idx3 = i + 3 * complex_per_native_vector_op;
                const float* current_src_float3 = srcp_float + complex_start_idx3 * 2;
                load_deinterleaved<float_vector_type>(real_parts_arr[3], imag_parts_arr[3], current_src_float3);
                abs_val_sq_arr[3] = mul_add(real_parts_arr[3], real_parts_arr[3], imag_parts_arr[3] * imag_parts_arr[3]);
                abs_val_arr[3] = sqrt(abs_val_sq_arr[3]);
                abs_val_arr[3] = abs_val_arr[3] + vcl_one;
                log_abs_val_arr[3] = log_ps_vcl_generic<float_vector_type>(abs_val_arr[3]);
                log_abs_val_arr[3].store_a(dstp + complex_start_idx3);
            }
        }

        for (int i = mod_length_unrolled; i < length; ++i)
        {
            float re = src[i].re;
            float im = src[i].im;
            dstp[i] = logf(sqrtf(re * re + im * im) + 1.0f);
        }
    }
} // namespace vcl_utils
