/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once

// to change the vectorization level go to the project properties:
// Configuration Properties > C/C++ > Code Generation > Enable Enhanced Instruction Set > Dropdown
// please do not define in this file as compilation in linux will fail
// 
// SSE3 is the minimum vectorization level supported, the code will not run on a machine not supporting it
// MSVC compiler does not provide the __FMA__ macro, but every Intel processor having __AVX2__ also has __FMA__:
#if !defined(__FMA__) && defined(__AVX2__)
#define __FMA__ 1
#endif

#if !defined(__SSE3__) && !defined(__AVX__) && !defined(__FMA__)  
#error Please specify exactly one of the following compiler switches to set the vectorization level: SSE3, AVX, FMA, AVX2
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

#if defined(__FMA__)
#pragma message("vectorization: FMA")
std::string vectorization_level = "vectorization: FMA";
#elif defined(__AVX__)
#pragma message("vectorization: AVX")
std::string vectorization_level = "vectorization: AVX";
#else
#define __SSE3__ 1
#pragma message("vectorization: SSE3")
std::string vectorization_level = "vectorization: SSE3";
#endif

#if defined(__GNUC__)
// disable gcc compiler warning for comparing different types of numbers
#pragma GCC diagnostic ignored "-Wsign-compare"
// disable gcc compiler warning for reaordering initializers
#pragma GCC diagnostic ignored "-Wreorder"
// disable gcc compiler warning for formatting of scanning text
#pragma GCC diagnostic ignored "-Wformat="
#endif

#include "types.hpp"
#include "helper_functions.hpp"

std::vector<std::shared_ptr<molecule>> molecules;
std::vector<calculation> calculations;

parameters param;
int n_ex = -1;