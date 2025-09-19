# Easing Cubic Bézier Function

Easy to integrate Cubic Bézier curve as Easing function library.

**License:** MIT. See [LICENSE.txt](LICENSE.txt)

# Problem

Cubic Bézier curves are widely used in animation systems to interpolate parameter values in keyframe animation techniques.
This library transforms Bézier curves into an explicit easing function form, eliminating the need to solve cubic polynomial equations at runtime—whether analytically or numerically.
Instead, it enables efficient parameter interpolation using transcendental functions, resulting in faster and more predictable execution.

# Library

The library is designed with C++20 and newer standards in mind, taking advantage of modern language features for clarity and performance. If needed, support for earlier versions of C++ can be added to ensure broader compatibility.

- Self-contained C++ library in single header file. No external dependencies other than standard C and C++ library.
- The entire concept is encapsulated in a single template class: `EasingCubicBezier<T>`. This class handles the interpolation of parameters used in the keyframe method. The interpolation of parameters follows the same principles as standard Bézier curve evaluation.
- The constructor takes the X and Y coordinates of the 4 control points of the cubic Bézier curve.
- Once instantiated, you can call the `evaluate` function with a parameter `t`, which should lie between x0 (the X coordinate of the first control point, representing the start time of the frame) and x3 (the X coordinate of the fourth control point, representing the end time). 
- The function returns the interpolated value at time t, based on the shape of the Bézier curve.
- The library can be used in two modes:
    - **PRECISE** – uses transcendental and irrational functions implemented in the <cmath> library, ensuring high accuracy.
    - **FAST** – uses fast approximations for selected transcendental and irrational functions, allowing for better performance at the expense of precision.
- Platform-independent, but developed and tested on Windows using Visual Studio and was tested in Linux using Clang.

## How to use

After calling either `find_package` or `add_subdirectory` simply link the library.
This automatically handles configuring the include directory. Example:

```cmake
find_package(EasingCubicBezier CONFIG REQUIRED)
target_link_libraries(YourGameEngine PRIVATE EasingCubicBezier)
```

For more info on using CMake visit the official [CMake documentation](https://cmake.org/cmake/help/latest/index.html).
