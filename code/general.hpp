/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once

// CMake sets one OPENCOMSORS_SIMD_* macro. If the project is built outside
// CMake, choose a conservative backend from compiler-provided architecture macros.
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

// to change the OpenMP settings go to the project properties to:
// Configuration Properties > C/C++ > Language > Open MP Support
#if defined(_OPENMP) && defined(PRINT_DEBUG_INFO)
#error If parallelization through OPENMP is enbled and PRINT_DEBUG_INFO should not be set 
#endif

#if defined(_DEBUG) || defined(DEBUG)
std::string compilation_mode = " compile mode: DEBUG";
#pragma message(" compile mode: DEBUG")
#else
std::string compilation_mode = " compile mode: RELEASE";
#pragma message(" compile mode: RELEASE")
#endif

#if defined(_OPENMP)
std::string OPENMP_parallelization = "       openmp: activated";
#pragma message("       openmp: activated")
#else
std::string OPENMP_parallelization = "      openmp: deactivated";
#pragma message("       openmp: deactivated")
#endif

#if defined(OPENCOMSORS_SIMD_FMA)
#pragma message("vectorization: FMA")
std::string vectorization_level = "vectorization: FMA";
#elif defined(OPENCOMSORS_SIMD_AVX)
#pragma message("vectorization: AVX")
std::string vectorization_level = "vectorization: AVX";
#elif defined(OPENCOMSORS_SIMD_SSE3)
#pragma message("vectorization: SSE3")
std::string vectorization_level = "vectorization: SSE3";
#elif defined(OPENCOMSORS_SIMD_NEON)
#pragma message("vectorization: NEON")
std::string vectorization_level = "vectorization: NEON";
#else
#pragma message("vectorization: SCALAR")
std::string vectorization_level = "vectorization: SCALAR";
#endif

#if defined(__GNUC__)
// disable gcc compiler warning for comparing different types of numbers
#pragma GCC diagnostic ignored "-Wsign-compare"
// disable gcc compiler warning for reordering initializers
#pragma GCC diagnostic ignored "-Wreorder"
// disable gcc compiler warning for formatting of scanning text
#pragma GCC diagnostic ignored "-Wformat="
#endif


#include <Eigen/Dense>

#if defined(USE_DOUBLE)
using calcType = double;
using MatrixCalcType = Eigen::MatrixXd;
using VectorCalcType = Eigen::VectorXd;
std::string precision = "    precision: double";
#pragma message("    precision: double")
#else
using calcType = float;
using MatrixCalcType = Eigen::MatrixXf;
using VectorCalcType = Eigen::VectorXf;
std::string precision = "    precision: single";
#pragma message("    precision: single")
#endif

#include "types.hpp"

std::vector<std::shared_ptr<molecule>> molecules;
std::vector<calculation> calculations;
std::vector<std::string> warnings;

parameters param;
int n_ex = -1;
