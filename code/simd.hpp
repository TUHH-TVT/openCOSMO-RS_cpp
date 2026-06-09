/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/

#pragma once

#if !defined(OPENCOMSORS_SIMD_FMA) && !defined(OPENCOMSORS_SIMD_AVX) && !defined(OPENCOMSORS_SIMD_SSE3) && !defined(OPENCOMSORS_SIMD_NEON) && !defined(OPENCOMSORS_SIMD_SCALAR)
#if defined(__FMA__)
#define OPENCOMSORS_SIMD_FMA 1
#elif defined(__AVX__)
#define OPENCOMSORS_SIMD_AVX 1
#elif defined(__SSE3__)
#define OPENCOMSORS_SIMD_SSE3 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#define OPENCOMSORS_SIMD_NEON 1
#else
#define OPENCOMSORS_SIMD_SCALAR 1
#endif
#endif

#if defined(OPENCOMSORS_SIMD_FMA) || defined(OPENCOMSORS_SIMD_AVX) || defined(OPENCOMSORS_SIMD_SSE3)
#include <immintrin.h>
#elif defined(OPENCOMSORS_SIMD_NEON)
#if defined(_MSC_VER) && defined(_M_ARM64)
#include <arm64_neon.h>
#else
#include <arm_neon.h>
#endif
#endif

template<typename T>
struct Simd;

#if defined(OPENCOMSORS_SIMD_FMA) || defined(OPENCOMSORS_SIMD_AVX)

template<>
struct Simd<float> {
    using Vec = __m256;
    static constexpr int lanes = 8;

    static Vec zero() { return _mm256_setzero_ps(); }
    static Vec load(const float* ptr) { return _mm256_loadu_ps(ptr); }
    static void store(float* ptr, Vec value) { _mm256_storeu_ps(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return _mm256_mul_ps(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return _mm256_add_ps(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) {
#if defined(OPENCOMSORS_SIMD_FMA)
        return _mm256_fmadd_ps(lhs, rhs, acc);
#else
        return add(mul(lhs, rhs), acc);
#endif
    }
    static float horizontal_sum(Vec value) {
        __m128 low = _mm256_castps256_ps128(value);
        __m128 high = _mm256_extractf128_ps(value, 1);
        __m128 sum = _mm_add_ps(low, high);
        __m128 shuf = _mm_movehdup_ps(sum);
        sum = _mm_add_ps(sum, shuf);
        shuf = _mm_movehl_ps(shuf, sum);
        sum = _mm_add_ss(sum, shuf);
        return _mm_cvtss_f32(sum);
    }
};

template<>
struct Simd<double> {
    using Vec = __m256d;
    static constexpr int lanes = 4;

    static Vec zero() { return _mm256_setzero_pd(); }
    static Vec load(const double* ptr) { return _mm256_loadu_pd(ptr); }
    static void store(double* ptr, Vec value) { _mm256_storeu_pd(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return _mm256_mul_pd(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return _mm256_add_pd(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) {
#if defined(OPENCOMSORS_SIMD_FMA)
        return _mm256_fmadd_pd(lhs, rhs, acc);
#else
        return add(mul(lhs, rhs), acc);
#endif
    }
    static double horizontal_sum(Vec value) {
        __m128d low = _mm256_castpd256_pd128(value);
        __m128d high = _mm256_extractf128_pd(value, 1);
        __m128d sum = _mm_add_pd(low, high);
        __m128d shuf = _mm_shuffle_pd(sum, sum, 0x1);
        sum = _mm_add_pd(sum, shuf);
        return _mm_cvtsd_f64(sum);
    }
};

#elif defined(OPENCOMSORS_SIMD_SSE3)

template<>
struct Simd<float> {
    using Vec = __m128;
    static constexpr int lanes = 4;

    static Vec zero() { return _mm_setzero_ps(); }
    static Vec load(const float* ptr) { return _mm_loadu_ps(ptr); }
    static void store(float* ptr, Vec value) { _mm_storeu_ps(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return _mm_mul_ps(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return _mm_add_ps(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return add(mul(lhs, rhs), acc); }
    static float horizontal_sum(Vec value) {
        __m128 shuf = _mm_movehdup_ps(value);
        __m128 sum = _mm_add_ps(value, shuf);
        shuf = _mm_movehl_ps(shuf, sum);
        sum = _mm_add_ss(sum, shuf);
        return _mm_cvtss_f32(sum);
    }
};

template<>
struct Simd<double> {
    using Vec = __m128d;
    static constexpr int lanes = 2;

    static Vec zero() { return _mm_setzero_pd(); }
    static Vec load(const double* ptr) { return _mm_loadu_pd(ptr); }
    static void store(double* ptr, Vec value) { _mm_storeu_pd(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return _mm_mul_pd(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return _mm_add_pd(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return add(mul(lhs, rhs), acc); }
    static double horizontal_sum(Vec value) {
        __m128d shuf = _mm_shuffle_pd(value, value, 0x1);
        __m128d sum = _mm_add_pd(value, shuf);
        return _mm_cvtsd_f64(sum);
    }
};

#elif defined(OPENCOMSORS_SIMD_NEON)

template<>
struct Simd<float> {
    using Vec = float32x4_t;
    static constexpr int lanes = 4;

    static Vec zero() { return vdupq_n_f32(0.0f); }
    static Vec load(const float* ptr) { return vld1q_f32(ptr); }
    static void store(float* ptr, Vec value) { vst1q_f32(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return vmulq_f32(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return vaddq_f32(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return vfmaq_f32(acc, lhs, rhs); }
    static float horizontal_sum(Vec value) { return vaddvq_f32(value); }
};

template<>
struct Simd<double> {
    using Vec = float64x2_t;
    static constexpr int lanes = 2;

    static Vec zero() { return vdupq_n_f64(0.0); }
    static Vec load(const double* ptr) { return vld1q_f64(ptr); }
    static void store(double* ptr, Vec value) { vst1q_f64(ptr, value); }
    static Vec mul(Vec lhs, Vec rhs) { return vmulq_f64(lhs, rhs); }
    static Vec add(Vec lhs, Vec rhs) { return vaddq_f64(lhs, rhs); }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return vfmaq_f64(acc, lhs, rhs); }
    static double horizontal_sum(Vec value) { return vaddvq_f64(value); }
};

#else

template<>
struct Simd<float> {
    using Vec = float;
    static constexpr int lanes = 1;

    static Vec zero() { return 0.0f; }
    static Vec load(const float* ptr) { return *ptr; }
    static void store(float* ptr, Vec value) { *ptr = value; }
    static Vec mul(Vec lhs, Vec rhs) { return lhs * rhs; }
    static Vec add(Vec lhs, Vec rhs) { return lhs + rhs; }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return add(mul(lhs, rhs), acc); }
    static float horizontal_sum(Vec value) { return value; }
};

template<>
struct Simd<double> {
    using Vec = double;
    static constexpr int lanes = 1;

    static Vec zero() { return 0.0; }
    static Vec load(const double* ptr) { return *ptr; }
    static void store(double* ptr, Vec value) { *ptr = value; }
    static Vec mul(Vec lhs, Vec rhs) { return lhs * rhs; }
    static Vec add(Vec lhs, Vec rhs) { return lhs + rhs; }
    static Vec fmadd(Vec lhs, Vec rhs, Vec acc) { return add(mul(lhs, rhs), acc); }
    static double horizontal_sum(Vec value) { return value; }
};

#endif
