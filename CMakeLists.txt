project(openCOSMORS)
cmake_minimum_required(VERSION 3.0)
set (CMAKE_CXX_STANDARD 14)

# configure pybind11 module
set(MY_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/code" 
    "${CMAKE_CURRENT_SOURCE_DIR}/eigen" 
    "${CMAKE_CURRENT_BINARY_DIR}"
)
add_subdirectory(pybind11)
pybind11_add_module(openCOSMORS "${CMAKE_CURRENT_SOURCE_DIR}/code/bindings_forPython.cpp")
set_target_properties(openCOSMORS PROPERTIES COMPILE_FLAGS -DPYBIND11)
target_include_directories(openCOSMORS PUBLIC ${MY_INCLUDES})

# please adjust the instruction set to your hardware by specifying one of the following
# compile flags: __SSE3__, __AVX__, __FMA__
# if you do not know what your hardware supports you can use the program cpu-z to check
set_target_properties(openCOSMORS PROPERTIES COMPILE_FLAGS -D__SSE3__)

# activate OpenMP for release builds if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(openCOSMORS PUBLIC $<$<CONFIG:RELEASE>:OpenMP::OpenMP_CXX>)
else()
    message("OpenMP was not found, to improve performance for release builds use a compiler supporting it.")
endif()