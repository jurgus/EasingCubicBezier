#include <benchmark/benchmark.h>
#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include "cubic_bezier.hpp"
#include <easingCubicBezier.hpp>

template<typename T>
T precise_cbrt(T phi)
{
    return std::cbrt(phi);
}
template<typename T>
T fast_cbrt_2(T phi)
{
    return fast_cbrt(phi);
}

template<typename T>
T precise_cosh_acosh(T phi)
{
    return phi >= T(1) ? std::cosh(std::acosh(phi) / T(3))
        : std::cos(std::acos(std::max(T(-1), phi)) / T(3));
}
template<typename T>
T fast_cosh_acosh(T phi)
{
    constexpr T one_third = T(1) / T(3);
    return phi >= T(1) ? std::cosh(fast_acosh(phi) * one_third)
        : std::cos(fast_acos(std::max(T(-1), phi)) * one_third);
}
template<typename T>
T fast_cosh_acosh_2(T phi)
{
    return phi >= T(1) ? std::cosh(fast_acosh(phi) / T(3))
        : std::cos(fast_acos(std::max(T(-1), phi)) / T(3));
}

template<typename T>
T precise_sinh_asinh(T phi)
{
    constexpr T one_third = T(1) / T(3);
    return std::sinh(std::asinh(phi) / T(3));
}
template<typename T>
T fast_sinh_asinh(T phi)
{
    constexpr T one_third = T(1) / T(3);
    return std::sinh(fast_asinh(phi) * one_third);
}
template<typename T>
T fast_sinh_asinh_2(T phi)
{
    return std::sinh(fast_asinh(phi) / T(3));
}

template<typename T>
T precise_cos_acos(T phi)
{
    constexpr T one_third = T(1) / T(3);
    constexpr T two_third_pi = -T(2) / T(3) * std::numbers::pi_v<T>;
    return std::cos(std::acos(std::clamp(phi, T(-1), T(1))) / T(3) + two_third_pi);
}
template<typename T>
T fast_cos_acos(T phi)
{
    constexpr T one_third = T(1) / T(3);
    constexpr T two_third_pi = -T(2) / T(3) * std::numbers::pi_v<T>;
    return std::cos(std::fma(fast_acos(std::clamp(phi, T(-1), T(1))), one_third, two_third_pi));
}
template<typename T>
T fast_cos_acos_2(T phi)
{
    constexpr T two_third_pi = -T(2) / T(3) * std::numbers::pi_v<T>;
    return std::cos(fast_acos(std::clamp(phi, T(-1), T(1))) / T(3) + two_third_pi);
}

constexpr int nSteps = 1001;

template<typename T>
static void BM_CorrectnessFunction(benchmark::State& state)
{
    T phi = state.range(0) / (T)10000;
    for (auto _ : state)
    {
        T max1 = T(0);
        T max2 = T(0);
        for (int i = 0; i < nSteps; ++i)
        {
            T t = i / T(nSteps - 1);
            auto result1 = complex_cosh_acosh_div3(t);
            auto result2 = real_cosh_acosh_div3(t);
            auto result3 = real_cosh_acosh_div3_2(t);
            max1 = std::max(max1, std::abs(result1 - result2));
            max2 = std::max(max2, std::abs(result1 - result3));
        }
        std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10)
            << "max1 " << max1 << '\n'
            << "max2 " << max2 << '\n';
    }
}

template<typename T>
static void Testing_Precise_Cos(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(precise_cos_acos(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Cos(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_cos_acos(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Cos_2(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_cos_acos_2(i / T(nSteps - 1)));
}

template<typename T>
static void Testing_Precise_Cosh(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(precise_cosh_acosh(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Cosh(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_cosh_acosh(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Cosh_2(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_cosh_acosh_2(i / T(nSteps - 1)));
}

template<typename T>
static void Testing_Precise_Sinh(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(precise_sinh_asinh(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Sinh(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_sinh_asinh(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Sinh_2(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_sinh_asinh_2(i / T(nSteps - 1)));
}

template<typename T>
static void Testing_Precise_Cbrt(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(precise_cbrt(i / T(nSteps - 1)));
}
template<typename T>
static void Testing_Fast_Cbrt(benchmark::State& state)
{
    for (auto _ : state)
        for (int i = 0; i < nSteps; ++i)
            benchmark::DoNotOptimize(fast_cbrt_2(i / T(nSteps - 1)));
}

#if 0
//BENCHMARK(BM_CorrectnessFunction<float>)->Name("Correctness/phi float")->Arg(0)->Iterations(1);
BENCHMARK(Testing_Precise_Cos<float>)->Name("Precise_Cos")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Cos<float>)->Name("Fast_Cos")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Cos_2<float>)->Name("Fast_Cos_2")->Iterations(1'000'000);
BENCHMARK(Testing_Precise_Cosh<float>)->Name("Precise_Cosh")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Cosh<float>)->Name("Fast_Cosh")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Cosh_2<float>)->Name("Fast_Cosh_2")->Iterations(1'000'000);
BENCHMARK(Testing_Precise_Sinh<float>)->Name("Precise_Sinh")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Sinh<float>)->Name("Fast_Sinh")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Sinh_2<float>)->Name("Fast_Sinh_2")->Iterations(1'000'000);
BENCHMARK(Testing_Precise_Cbrt<float>)->Name("Precise_Cbrt")->Iterations(1'000'000);
BENCHMARK(Testing_Fast_Cbrt<float>)->Name("Fast_Cbrt")->Iterations(1'000'000);
//BENCHMARK(BM_CorrectnessFunction<double>)->Name("Correctness/phi double")->Apply(CustomArguments)->Iterations(1);
//BENCHMARK(BM_Optimized_Pos<double>)->Name("Optimized/phi double")->Apply(CustomArguments)->Iterations(1'000'000);
//BENCHMARK(BM_Optimized_Pos2<double>)->Name("Optimized2/phi double")->Apply(CustomArguments)->Iterations(1'000'000);
//BENCHMARK(BM_Complex_Pos<double>)->Name("Complex  /phi double")->Apply(CustomArguments)->Iterations(1'000'000);
#endif

#if 0
static std::array<float, 4> P_Y = { 0.1f, 0.5f, 0.5f, 1.0f };
static std::array<float, 4> P_X0 = { 0.0f, 2.0f / 6.0f, 4.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X1 = { 0.0f, 1.0f / 6.0f, 3.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X2 = { 0.0f, 0.0f / 6.0f, 0.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X3 = { 0.0f, 1.0f / 6.0f, 5.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X4 = { 0.0f, 3.0f / 6.0f, 4.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X5 = { 0.0f, 5.0f / 6.0f, 6.0f / 6.0f, 1.0f };
static std::array<float, 4> P_X6 = { 0.0f, 0.0f / 6.0f, 1.0f / 6.0f, 1.0f };

static EasingCubicBezierf ecb_f_arr[] =
{
    EasingCubicBezierf(P_X0, P_Y),
    EasingCubicBezierf(P_X1, P_Y),
    EasingCubicBezierf(P_X2, P_Y),
    EasingCubicBezierf(P_X3, P_Y),
    EasingCubicBezierf(P_X4, P_Y),
    EasingCubicBezierf(P_X5, P_Y),
    EasingCubicBezierf(P_X6, P_Y),
};
static std::array<float, 4> P_X_arr[] =
{
    P_X0,
    P_X1,
    P_X2,
    P_X3,
    P_X4,
    P_X5,
    P_X6,
};

static std::array<double, 4> dP_Y = { 0.0, 0.5, 0.5, 1.0 };
static std::array<double, 4> dP_X0 = { 0.0, 2.0 / 6.0, 4.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X1 = { 0.0, 1.0 / 6.0, 3.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X2 = { 0.0, 0.0 / 6.0, 0.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X3 = { 0.0, 1.0 / 6.0, 5.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X4 = { 0.0, 3.0 / 6.0, 4.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X5 = { 0.0, 5.0 / 6.0, 6.0 / 6.0, 1.0 };
static std::array<double, 4> dP_X6 = { 0.0, 0.0 / 6.0, 1.0 / 6.0, 1.0 };

static EasingCubicBezierd decb_f_arr[] =
{
    EasingCubicBezierd(dP_X0, dP_Y),
    EasingCubicBezierd(dP_X1, dP_Y),
    EasingCubicBezierd(dP_X2, dP_Y),
    EasingCubicBezierd(dP_X3, dP_Y),
    EasingCubicBezierd(dP_X4, dP_Y),
    EasingCubicBezierd(dP_X5, dP_Y),
    EasingCubicBezierd(dP_X6, dP_Y),
};
static std::array<double, 4> dP_X_arr[] =
{
    dP_X0,
    dP_X1,
    dP_X2,
    dP_X3,
    dP_X4,
    dP_X5,
    dP_X6,
};

template<typename T>
static void BM_Easing_Bezier_Cubic_Correctness(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    for (auto _ : state)
    {
        T max1 = T(0);
        T max2 = T(0);
        T max3 = T(0);
        T max4 = T(0);
        std::ostringstream vals;
        vals << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
        for (int i = 0; i < steps; ++i)
        {
            T t = i / T(steps - 1);
            T x0, x1, x2, x3, x4;
            if constexpr (std::is_same_v<T, float>)
            {
                x0 = ecb_f_arr[n].evaluate(t);
                x1 = BlenderOrg::fcurve_eval_keyframes_interpolate(P_X_arr[n], P_Y, t);
                x2 = BlenderOptim::fcurve_eval_keyframes_interpolate(P_X_arr[n], P_Y, t);
                x3 = BlenderNumeric::fcurve_eval_keyframes_interpolate1(P_X_arr[n], P_Y, t);
                x4 = BlenderNumeric::fcurve_eval_keyframes_interpolate2(P_X_arr[n], P_Y, t);
            }
            if constexpr (std::is_same_v<T, double>)
            {
                x0 = decb_f_arr[n].evaluate(t);
                x1 = BlenderOrg::fcurve_eval_keyframes_interpolate(dP_X_arr[n], dP_Y, t);
                x2 = BlenderOptim::fcurve_eval_keyframes_interpolate(dP_X_arr[n], dP_Y, t);
                x3 = BlenderNumeric::fcurve_eval_keyframes_interpolate1(dP_X_arr[n], dP_Y, t);
                x4 = BlenderNumeric::fcurve_eval_keyframes_interpolate2(dP_X_arr[n], dP_Y, t);
            }
            vals << t << ';' << x0 << ';' << x1 << ';' << x2 << ';' << x3 << ';' << x4 << '\n';
            max1 = std::max(max1, std::abs(x0 - x1));
            max2 = std::max(max2, std::abs(x0 - x2));
            max3 = std::max(max3, std::abs(x0 - x3));
            max4 = std::max(max3, std::abs(x0 - x4));
        }
        std::cout << vals.str() << std::endl;
        std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10) 
            << "max1 " << max1 << '\n'
            << "max2 " << max2 << '\n'
            << "max3 " << max3 << '\n'
            << "max4 " << max4 << '\n';
    }
}

template<typename T>
static void BM_Easing_Bezier_Cubic(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    EasingCubicBezier<T> data;
    if constexpr (std::is_same_v<T, float>)
        data = ecb_f_arr[n];
    if constexpr (std::is_same_v<T, double>)
        data = decb_f_arr[n];
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(data.evaluate(t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    std::array<T, 4> data;
    std::array<T, 4> dataY;
    if constexpr (std::is_same_v<T, float>)
    {
        data = P_X_arr[n];
        dataY = P_Y;
    }
    if constexpr (std::is_same_v<T, double>)
    {
        data = dP_X_arr[n];
        dataY = dP_Y;
    }
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderOrg::fcurve_eval_keyframes_interpolate(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Optim(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    std::array<T, 4> data;
    std::array<T, 4> dataY;
    if constexpr (std::is_same_v<T, float>)
    {
        data = P_X_arr[n];
        dataY = P_Y;
    }
    if constexpr (std::is_same_v<T, double>)
    {
        data = dP_X_arr[n];
        dataY = dP_Y;
    }
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderOptim::fcurve_eval_keyframes_interpolate(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Numeric1(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    std::array<T, 4> data;
    std::array<T, 4> dataY;
    if constexpr (std::is_same_v<T, float>)
    {
        data = P_X_arr[n];
        dataY = P_Y;
    }
    if constexpr (std::is_same_v<T, double>)
    {
        data = dP_X_arr[n];
        dataY = dP_Y;
    }
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderNumeric::fcurve_eval_keyframes_interpolate1(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Numeric2(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t steps = state.range(1);
    std::array<T, 4> data;
    std::array<T, 4> dataY;
    if constexpr (std::is_same_v<T, float>)
    {
        data = P_X_arr[n];
        dataY = P_Y;
    }
    if constexpr (std::is_same_v<T, double>)
    {
        data = dP_X_arr[n];
        dataY = dP_Y;
    }
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderNumeric::fcurve_eval_keyframes_interpolate2(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}

constexpr int Iterations = 1'00'000;
constexpr int kSteps = 1001;

//#define REP ->Repetitions(5)->ReportAggregatesOnly(true)
#define REP
#define ARG ->Args({0, kSteps})->Args({1, kSteps})->Args({2, kSteps})->Args({3, kSteps})->Args({4, kSteps})->Args({5, kSteps})->Args({6, kSteps})
//#define ARG ->Args({4, kSteps})

BENCHMARK(BM_Easing_Bezier_Cubic_Correctness<float>)
->Name("Correctness/phi float")
->Iterations(1)
ARG
;
BENCHMARK(BM_Easing_Bezier_Cubic<float>)
->Name("EasingCubicBezierf ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender_Numeric1<float>)
->Name("BlenderNumeric ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender_Numeric2<float>)
->Name("BlenderNumeric ")
->Iterations(Iterations)
REP
ARG
;
//BENCHMARK(BM_Easing_Bezier_Blender<float>)
//->Name("Blender ")
//->Iterations(Iterations)
//REP
//ARG
//;
//BENCHMARK(BM_Easing_Bezier_Blender_Optim<float>)
//->Name("BlenderOptim ")
//->Iterations(Iterations)
//REP
//ARG
//;

//BENCHMARK(BM_Easing_Bezier_Cubic_Correctness<double>)
//->Name("Correctness/phi double")
//->Iterations(1)
//ARG
//
//BENCHMARK(BM_Easing_Bezier_Cubic<double>)
//->Name("EasingCubicBezierd ")
//->Iterations(Iterations)
//REP
//ARG
//;
//BENCHMARK(BM_Easing_Bezier_Blender<double>)
//->Name("Blender ")
//->Iterations(Iterations)
//REP
//ARG
//;
//BENCHMARK(BM_Easing_Bezier_Blender_Optim<double>)
//->Name("BlenderOptim ")
//->Iterations(Iterations)
//REP
//ARG
//;
//BENCHMARK(BM_Easing_Bezier_Blender_Numeric1<double>)
//->Name("BlenderNumeric ")
//->Iterations(Iterations)
//REP
//ARG
//;
#endif

#if 1
constexpr int Iterations = 1'0'000;
constexpr int kSteps = 1001;

#define K T(15.0)
#define W 15

template<typename T>
static void BM_Easing_Bezier_Cubic(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t m = state.range(1);
    int64_t steps = state.range(2);
    std::array<T, 4> P_Y = { T(0.0), T(0.5), T(0.5), T(1.0) };
    std::array<T, 4> P_X = { T(0.0), T(n) / K, T(m) / K, T(1.0) };
    EasingCubicBezier<T> data(P_X, P_Y);
//    std::cout << data.mL << "  " << data.mK << '\n';
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(data.evaluate(t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Numeric1(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t m = state.range(1);
    int64_t steps = state.range(2);
    std::array<T, 4> dataY = { T(0.0), T(0.5), T(0.5), T(1.0) };
    std::array<T, 4> data = { T(0.0), T(n) / K, T(m) / K, T(1.0) };
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderNumeric::fcurve_eval_keyframes_interpolate1(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Numeric2(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t m = state.range(1);
    int64_t steps = state.range(2);
    std::array<T, 4> dataY = { T(0.0), T(0.5), T(0.5), T(1.0) };
    std::array<T, 4> data = { T(0.0), T(n) / K, T(m) / K, T(1.0) };
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderNumeric::fcurve_eval_keyframes_interpolate2(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t m = state.range(1);
    int64_t steps = state.range(2);
    std::array<T, 4> dataY = { T(0.0), T(0.5), T(0.5), T(1.0) };
    std::array<T, 4> data = { T(0.0), T(n) / K, T(m) / K, T(1.0) };
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderOrg::fcurve_eval_keyframes_interpolate(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
template<typename T>
static void BM_Easing_Bezier_Blender_Optim(benchmark::State& state)
{
    int64_t n = state.range(0);
    int64_t m = state.range(1);
    int64_t steps = state.range(2);
    std::array<T, 4> dataY = { T(0.0), T(0.5), T(0.5), T(1.0) };
    std::array<T, 4> data = { T(0.0), T(n) / K, T(m) / K, T(1.0) };
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(steps);
        for (int i = 0; i < steps; ++i)
        {
            T t = static_cast<T>(i) / static_cast<T>(steps - 1);
            benchmark::DoNotOptimize(t);
            benchmark::DoNotOptimize(BlenderOptim::fcurve_eval_keyframes_interpolate(data, dataY, t));
        }
        benchmark::ClobberMemory();
    }
}
void CustomArguments(benchmark::internal::Benchmark* b)
{
    for (int n = 0; n <= W; ++n)
    {
        for (int m = 0; m <= W; ++m)
        {
            b->Args({ n, m, kSteps });
        }
    }
}

//#define REP ->Repetitions(5)->ReportAggregatesOnly(true)
#define REP
#define ARG ->Apply(CustomArguments)
//#define ARG ->Args({4, kSteps})

BENCHMARK(BM_Easing_Bezier_Cubic<float>)
->Name("EasingCubicBezierf ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender_Numeric1<float>)
->Name("Blender Numeric1 ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender_Numeric2<float>)
->Name("Blender Numeric2 ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender<float>)
->Name("Blender Orginal ")
->Iterations(Iterations)
REP
ARG
;
BENCHMARK(BM_Easing_Bezier_Blender_Optim<float>)
->Name("Blender Optim ")
->Iterations(Iterations)
REP
ARG
;
#endif

//UseManualTime()
BENCHMARK_MAIN();
