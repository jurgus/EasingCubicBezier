#include <benchmark/benchmark.h>
#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include "cubic_bezier.hpp"
#include <easingCubicBezier.hpp>

static std::vector<int64_t> args =
{
    -100000,
    -10000,
    -1000,
    -100,
    -10,
    -1,
    +0,
    +1,
    +10,
    +100,
    +1000,
    +10000,
    +100000
};

static void CustomArguments(benchmark::internal::Benchmark* b)
{
    for (auto i : args)
        b->Arg(i);
}

template<typename T>
T complex_cosh_acosh_div3(T phi)
{
    constexpr T one_third = T(1) / T(3);
    return std::cosh(std::acosh(phi) * one_third);
}

template<typename T>
T real_cosh_acosh_div3(T phi)
{
    constexpr T one_third = T(1) / T(3);
    return std::cos(fast_acosh(phi) * one_third);
}

template<typename T>
T real_cosh_acosh_div3_2(T phi)
{
    return std::cos(fast_acosh(phi) / T(3));
}

template<typename T>
static void BM_CorrectnessFunction(benchmark::State& state)
{
    T phi = state.range(0) / (T)10000;
    for (auto _ : state)
    {
        auto result1 = complex_cosh_acosh_div3(phi);
        auto result2 = real_cosh_acosh_div3(phi);
        if (std::abs(result1 - result2) > std::numeric_limits<double>::epsilon())
        {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(std::numeric_limits<float>::digits10) << result1 - result2;
            state.SkipWithError("Output mismatch! " + oss.str());
            break;
        }
    }
}

template<typename T>
static void BM_Optimized_Pos(benchmark::State& state)
{
    T phi = state.range(0) / (T)10000;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(real_cosh_acosh_div3(phi));
    }
}
template<typename T>
static void BM_Optimized_Pos2(benchmark::State& state)
{
    T phi = state.range(0) / (T)10000;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(real_cosh_acosh_div3_2(phi));
    }
}
template<typename T>
static void BM_Complex_Pos(benchmark::State& state)
{
    T phi = state.range(0) / (T)10000;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(complex_cosh_acosh_div3(phi));
    }
}

#if 0
//BENCHMARK(BM_CorrectnessFunction<float>)->Name("Correctness/phi float")->Apply(CustomArguments)->Iterations(1);
BENCHMARK(BM_Optimized_Pos<float>)->Name("Optimized/phi float")->Apply(CustomArguments)->Iterations(1'000'000);
BENCHMARK(BM_Optimized_Pos2<float>)->Name("Optimized2/phi float")->Apply(CustomArguments)->Iterations(1'000'000);
BENCHMARK(BM_Complex_Pos<float>)->Name("Complex  /phi float")->Apply(CustomArguments)->Iterations(1'000'000);
//BENCHMARK(BM_CorrectnessFunction<double>)->Name("Correctness/phi double")->Apply(CustomArguments)->Iterations(1);
BENCHMARK(BM_Optimized_Pos<double>)->Name("Optimized/phi double")->Apply(CustomArguments)->Iterations(1'000'000);
BENCHMARK(BM_Optimized_Pos2<double>)->Name("Optimized2/phi double")->Apply(CustomArguments)->Iterations(1'000'000);
BENCHMARK(BM_Complex_Pos<double>)->Name("Complex  /phi double")->Apply(CustomArguments)->Iterations(1'000'000);
#endif

static std::array<float, 4> P_Y = { 0.0f, 0.5f, 0.5f, 1.0f };
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
constexpr int kSteps = 1001;

template<typename T>
static void BM_Easing_Bezier_Cubic_Correctness(benchmark::State& state)
{
    int64_t n = state.range(0);
    for (auto _ : state)
    {
        T max1 = T(0);
        T max2 = T(0);
        T max3 = T(0);
        std::ostringstream vals;
        //vals << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
        for (int i = 0; i < kSteps; ++i)
        {
            T t = i / T(kSteps - 1);
            T x0 = ecb_f_arr[n].evaluate(t);
            T x1 = BlenderOrg::fcurve_eval_keyframes_interpolate(P_X_arr[n], P_Y, t);
            T x2 = BlenderOptim::fcurve_eval_keyframes_interpolate(P_X_arr[n], P_Y, t);
            T x3 = BlenderNumeric::fcurve_eval_keyframes_interpolate1(P_X_arr[n], P_Y, t);
            // vals << t << ';' << x0 << ';' << x1 << ';' << x2 << ';' << x3 << '\n';
            if (std::abs(x0 - x1) > std::numeric_limits<T>::epsilon())
            {
                max1 = std::max(max1, std::abs(x0 - x1));
            }
            if (std::abs(x0 - x2) > std::numeric_limits<T>::epsilon())
            {
                max2 = std::max(max1, std::abs(x0 - x2));
            }
            if (std::abs(x0 - x3) > std::numeric_limits<T>::epsilon())
            {
                max3 = std::max(max1, std::abs(x0 - x3));
            }
        }
        std::cout << vals.str() << std::endl;
        std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10) 
            << "max1 " << max1 << '\n'
            << "max2 " << max2 << '\n'
            << "max3 " << max3 << '\n';
        benchmark::ClobberMemory();
    }
}

template<typename T>
static void BM_Easing_Bezier_Cubic(benchmark::State& state)
{
    int64_t n = state.range(0);
    const auto& data = ecb_f_arr[n];
    for (auto _ : state)
    {
        for (int i = 0; i < kSteps; ++i)
        {
            float t = i / float(kSteps - 1);
            benchmark::DoNotOptimize(data.evaluate(t));
        }
        benchmark::ClobberMemory();
    }
    //state.SetItemsProcessed(state.iterations() * kSteps);
}
template<typename T>
static void BM_Easing_Bezier_Blender(benchmark::State& state)
{
    int64_t n = state.range(0);
    const auto& data = P_X_arr[n];
    for (auto _ : state)
    {
        for (int i = 0; i < kSteps; ++i)
        {
            float t = i / float(kSteps - 1);
            benchmark::DoNotOptimize(BlenderOrg::fcurve_eval_keyframes_interpolate(data, P_Y, t));
        }
        benchmark::ClobberMemory();
    }
    //state.SetItemsProcessed(state.iterations() * kSteps);
}
template<typename T>
static void BM_Easing_Bezier_Blender_Optim(benchmark::State& state)
{
    int64_t n = state.range(0);
    const auto& data = P_X_arr[n];
    for (auto _ : state)
    {
        for (int i = 0; i < kSteps; ++i)
        {
            float t = i / float(kSteps - 1);
            benchmark::DoNotOptimize(BlenderOptim::fcurve_eval_keyframes_interpolate(data, P_Y, t));
        }
        benchmark::ClobberMemory();
    }
    //state.SetItemsProcessed(state.iterations() * kSteps);
}
template<typename T>
static void BM_Easing_Bezier_Blender_Numeric1(benchmark::State& state)
{
    int64_t n = state.range(0);
    const auto& data = P_X_arr[n];
    for (auto _ : state)
    {
        for (int i = 0; i < kSteps; ++i)
        {
            float t = i / float(kSteps - 1);
            benchmark::DoNotOptimize(BlenderNumeric::fcurve_eval_keyframes_interpolate1(data, P_Y, t));
        }
        benchmark::ClobberMemory();
    }
    //state.SetItemsProcessed(state.iterations() * kSteps);
}

#if 1

// #define REP ->Repetitions(5)->ReportAggregatesOnly(true)
#define REP

constexpr int Iterations = 1'00'000;

BENCHMARK(BM_Easing_Bezier_Cubic_Correctness<float>)
->Name("Correctness/phi float")
->Iterations(1)
->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6);

BENCHMARK(BM_Easing_Bezier_Cubic<float>)
->Name("EasingCubicBezierf ")
->Iterations(Iterations)
REP
->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)
;
BENCHMARK(BM_Easing_Bezier_Blender<float>)
->Name("Blender ")
->Iterations(Iterations)
REP
->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)
;
BENCHMARK(BM_Easing_Bezier_Blender_Optim<float>)
->Name("BlenderOptim ")
->Iterations(Iterations)
REP
->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)
;
BENCHMARK(BM_Easing_Bezier_Blender_Numeric1<float>)
->Name("BlenderNumeric ")
->Iterations(Iterations)
REP
->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)
;
#endif
//UseManualTime()
BENCHMARK_MAIN();
