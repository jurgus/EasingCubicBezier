cmake_minimum_required(VERSION 3.10)
project(EasingCubicBezier VERSION 0.9.0 LANGUAGES CXX)

#----------------------------------------------------------------#
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------#
option(USE_AVX2 "Enable AVX2 instructions" ON)
option(USE_AVX "Enable AVX instructions" OFF)
option(USE_SSE "Enable SSE4.2 instructions" OFF)
option(USE_FAST_MATH "Enable fast floating-point math" OFF)
option(USE_PRECISE_MATH "Enable precise floating-point math" ON)

#----------------------------------------------------------------#
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(USE_AVX2)
        add_compile_options(-mavx2)
    elseif(USE_AVX)
        add_compile_options(-mavx)
    elseif(USE_SSE)
        add_compile_options(-msse4.2)
    endif()
    if(USE_FAST_MATH)
        add_compile_options(-ffast-math)
    elseif(USE_PRECISE_MATH)
        add_compile_options(-fno-fast-math)
        add_compile_definitions(FP_PRECISE)
    endif()
    # add_compile_options(-march=native)
elseif(MSVC)
    add_compile_options(/Ob2 /Oi /Ot)
    if(USE_AVX2)
        add_compile_options(/arch:AVX2)
    elseif(USE_AVX)
        add_compile_options(/arch:AVX)
    elseif(USE_SSE)
        add_compile_options(/arch:SSE2)
    endif()
    if(USE_FAST_MATH)
        add_compile_options(/fp:fast)
    elseif(USE_PRECISE_MATH)
        add_compile_options(/fp:precise)
        add_compile_definitions(FP_PRECISE)
    endif()
endif()
#----------------------------------------------------------------#
add_subdirectory(external/benchmark)

add_library(EasingCubicBezier INTERFACE)

target_include_directories(EasingCubicBezier
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

add_executable(benchmark_tests tests/benchmark.cpp tests/cubic_bezier.hpp )
target_link_libraries(benchmark_tests PRIVATE benchmark::benchmark EasingCubicBezier)

enable_testing()
add_executable(test_cubic_bezier tests/test_cubic_bezier.cpp  tests/cubic_bezier.hpp)
target_link_libraries(test_cubic_bezier PRIVATE EasingCubicBezier)
