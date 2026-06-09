cmake_minimum_required(VERSION 3.10)
project(openCOSMORS)
set (CMAKE_CXX_STANDARD 14)
set(ignoreMe "${BINARY}")
set(OPENCOMSORS_SIMD "AUTO" CACHE STRING "SIMD backend: AUTO, SCALAR, SSE3, AVX, FMA, NEON")
string(TOUPPER "${OPENCOMSORS_SIMD}" OPENCOMSORS_SIMD)
if(NOT OPENCOMSORS_SIMD MATCHES "^(AUTO|SCALAR|SSE3|AVX|FMA|NEON)$")
  message(FATAL_ERROR "OPENCOMSORS_SIMD must be one of AUTO, SCALAR, SSE3, AVX, FMA, NEON")
endif()

# configure pybind11 module
set(MY_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/code" 
    "${CMAKE_CURRENT_BINARY_DIR}"
)
set(MY_SYSTEM_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/eigen"
)

if(NOT DEFINED BINARY)
  message("-----------------CREATE A PYTHON MODULE-----------------")

  add_subdirectory(pybind11)
  pybind11_add_module(openCOSMORS "${CMAKE_CURRENT_SOURCE_DIR}/code/bindings_forPython.cpp")
  target_compile_definitions(openCOSMORS PRIVATE PYBIND11)

else()
  message("-----------------CREATE AN STANDALONE BINARY-----------------")
  add_executable(openCOSMORS "${CMAKE_CURRENT_SOURCE_DIR}/code/bindings_forCLI.cpp")
endif()

target_include_directories(openCOSMORS PUBLIC ${MY_INCLUDES})
target_include_directories(openCOSMORS SYSTEM PUBLIC ${MY_SYSTEM_INCLUDES})

if(OPENCOMSORS_SIMD STREQUAL "AUTO")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    set(OPENCOMSORS_SELECTED_SIMD "NEON")
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$")
    set(OPENCOMSORS_SELECTED_SIMD "FMA")
  else()
    set(OPENCOMSORS_SELECTED_SIMD "SCALAR")
  endif()
else()
  set(OPENCOMSORS_SELECTED_SIMD "${OPENCOMSORS_SIMD}")
endif()

message("openCOSMORS SIMD backend: ${OPENCOMSORS_SELECTED_SIMD}")
target_compile_definitions(openCOSMORS PRIVATE "OPENCOMSORS_SIMD_${OPENCOMSORS_SELECTED_SIMD}")

if(OPENCOMSORS_SELECTED_SIMD STREQUAL "SSE3")
  if(MSVC)
    target_compile_definitions(openCOSMORS PRIVATE __SSE3__)
  else()
    target_compile_options(openCOSMORS PRIVATE -msse3)
  endif()
elseif(OPENCOMSORS_SELECTED_SIMD STREQUAL "AVX")
  if(MSVC)
    target_compile_options(openCOSMORS PRIVATE /arch:AVX)
  else()
    target_compile_options(openCOSMORS PRIVATE -mavx)
  endif()
elseif(OPENCOMSORS_SELECTED_SIMD STREQUAL "FMA")
  if(MSVC)
    target_compile_options(openCOSMORS PRIVATE /arch:AVX2)
  else()
    target_compile_options(openCOSMORS PRIVATE -mavx -mfma)
  endif()
elseif(OPENCOMSORS_SELECTED_SIMD STREQUAL "NEON")
  if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    message(FATAL_ERROR "OPENCOMSORS_SIMD=NEON currently supports AArch64/ARM64 targets only")
  endif()
endif()

# activate OpenMP for release builds if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(openCOSMORS PUBLIC $<$<CONFIG:RELEASE>:OpenMP::OpenMP_CXX>)
else()
    message("OpenMP was not found, to improve performance for release builds use a compiler supporting it.")
endif()
