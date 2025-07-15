cmake_minimum_required(VERSION 3.10)
project(EasingCubicBezier VERSION 0.9.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(EasingCubicBezier INTERFACE)

target_include_directories(EasingCubicBezier
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)


enable_testing()
add_executable(test_cubic_bezier tests/test_cubic_bezier.cpp)
target_link_libraries(test_cubic_bezier PRIVATE EasingCubicBezier)
