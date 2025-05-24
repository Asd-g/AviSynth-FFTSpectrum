#pragma once

#include "../VCL2/vectorclass.h"
#include <avs/config.h>

namespace vcl_log_constants
{
    template<typename float_vector_type>
    struct corresponding_int_vector;

    template<>
    struct corresponding_int_vector<Vec4f>
    {
        using type = Vec4i;
    };

    template<>
    struct corresponding_int_vector<Vec8f>
    {
        using type = Vec8i;
    };

    template<>
    struct corresponding_int_vector<Vec16f>
    {
        using type = Vec16i;
    };

    template<typename float_vector_type>
    using int_vec_t = typename corresponding_int_vector<float_vector_type>::type;

    template<typename float_vector_type>
    struct log_ps
    {
        using int_vec_type = int_vec_t<float_vector_type>;

        inline static const float_vector_type one = float_vector_type(1.0f);
        inline static const float_vector_type p0_5 = float_vector_type(0.5f);

        inline static const float_vector_type min_norm_pos = reinterpret_f(int_vec_type(0x00800000));
        inline static const float_vector_type inv_mant_mask = reinterpret_f(int_vec_type(~0x7f800000));

        inline static const int_vec_type pi32_0x7f = int_vec_type(0x7f);

        inline static const float_vector_type cephes_SQRTHF = float_vector_type(0.707106781186547524f);
        inline static const float_vector_type cephes_log_p0 = float_vector_type(7.0376836292E-2f);
        inline static const float_vector_type cephes_log_p1 = float_vector_type(-1.1514610310E-1f);
        inline static const float_vector_type cephes_log_p2 = float_vector_type(1.1676998740E-1f);
        inline static const float_vector_type cephes_log_p3 = float_vector_type(-1.2420140846E-1f);
        inline static const float_vector_type cephes_log_p4 = float_vector_type(1.4249322787E-1f);
        inline static const float_vector_type cephes_log_p5 = float_vector_type(-1.6668057665E-1f);
        inline static const float_vector_type cephes_log_p6 = float_vector_type(2.0000714765E-1f);
        inline static const float_vector_type cephes_log_p7 = float_vector_type(-2.4999993993E-1f);
        inline static const float_vector_type cephes_log_p8 = float_vector_type(3.3333331174E-1f);
        inline static const float_vector_type cephes_log_q1 = float_vector_type(-2.12194440e-4f);
        inline static const float_vector_type cephes_log_q2 = float_vector_type(0.693359375f);
    };

} // namespace vcl_log_constants

template<typename float_vector_type>
struct corresponding_bool_vector;

template<>
struct corresponding_bool_vector<Vec4f>
{
    using type = Vec4fb;
};

template<>
struct corresponding_bool_vector<Vec8f>
{
    using type = Vec8fb;
};

template<>
struct corresponding_bool_vector<Vec16f>
{
    using type = Vec16fb;
};

template<typename float_vector_type>
using bool_vec_t = typename corresponding_bool_vector<float_vector_type>::type;

template<typename float_vector_type>
struct corresponding_uint_vector;

template<>
struct corresponding_uint_vector<Vec4f>
{
    using type = Vec4ui;
};

template<>
struct corresponding_uint_vector<Vec8f>
{
    using type = Vec8ui;
};

template<>
struct corresponding_uint_vector<Vec16f>
{
    using type = Vec16ui;
};

template<typename float_vector_type>
using uint_vec_t = typename corresponding_uint_vector<float_vector_type>::type;

template<typename float_vector_type>
AVS_FORCEINLINE float_vector_type log_ps_vcl_generic(float_vector_type x_input)
{
    using constants = vcl_log_constants::log_ps<float_vector_type>;
    using int_vec_type = vcl_log_constants::int_vec_t<float_vector_type>;
    using uint_vec_type = uint_vec_t<float_vector_type>;
    using bool_vec_type = bool_vec_t<float_vector_type>;

    const float_vector_type& one = constants::one;

    bool_vec_type invalid_mask = (x_input <= 0.0f);
    float_vector_type x = max(x_input, constants::min_norm_pos);

    int_vec_type x_int_representation = reinterpret_i(x);
    int_vec_type emm0 = int_vec_type(uint_vec_type(x_int_representation) >> 23);

    int_vec_type inv_mant_mask_int = reinterpret_i(constants::inv_mant_mask);
    x_int_representation = x_int_representation & inv_mant_mask_int;

    int_vec_type const_0p5_int = reinterpret_i(constants::p0_5);
    x_int_representation = x_int_representation | const_0p5_int;

    x = reinterpret_f(x_int_representation);

    emm0 = emm0 - constants::pi32_0x7f;
    float_vector_type e_float = to_float(emm0);
    e_float = e_float + one;

    float_vector_type f;
    {
        bool_vec_type mask_lt_sqrthf = (x < constants::cephes_SQRTHF);
        float_vector_type x_if_lt = x;

        f = x - one;
        e_float = if_sub(mask_lt_sqrthf, e_float, one);
        f = select(mask_lt_sqrthf, f + x_if_lt, f);
    }

    float_vector_type poly_z, poly_y;
    poly_z = f * f;

    poly_y = constants::cephes_log_p0;
    poly_y = mul_add(poly_y, f, constants::cephes_log_p1);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p2);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p3);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p4);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p5);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p6);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p7);
    poly_y = mul_add(poly_y, f, constants::cephes_log_p8);
    poly_y = poly_y * f;
    poly_y = poly_y * poly_z;

    float_vector_type tmp;
    tmp = e_float * constants::cephes_log_q1;
    poly_y = poly_y + tmp;
    tmp = poly_z * constants::p0_5;
    poly_y = poly_y - tmp;
    tmp = e_float * constants::cephes_log_q2;
    f = f + poly_y;
    f = f + tmp;

    float_vector_type result = select(invalid_mask, nan_vec<float_vector_type>(0x101), f);
    return result;
}
