#pragma once

#include <cmath>
#include <array>
#include <algorithm>
#include <numbers>
#include <limits>

#if defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#else
#define FORCE_INLINE inline
#endif

#if defined(__cpp_lib_unreachable) && __cpp_lib_unreachable >= 202202L
#include <utility>
#define UNREACHABLE() std::unreachable()
#else
#if defined(_MSC_VER) && !defined(__clang__) // MSVC
#define UNREACHABLE() __assume(false)
#else // GCC, Clang, others
#define UNREACHABLE() __builtin_unreachable()
#endif
#endif

template<typename T>
inline bool isZero(T x, T epsilon = std::numeric_limits<T>::epsilon()) noexcept
{
    return std::fabs(x) < epsilon;
}

template<typename T>
inline T fast_cbrt(T x) noexcept
{
    constexpr T one_third = T(1) / T(3);
    return x == T(0) ? T(0) : std::copysign(T(1), x) * std::exp(std::log(std::fabs(x)) * one_third);
}

template<typename T>
inline T fast_asinh(T x) noexcept
{
#if defined(_MSC_VER) && !defined(__AVX2__)
    return std::log(x + std::sqrt(x * x + T(1)));
#else
    return std::log(x + std::sqrt(std::fma(x, x, T(1))));
#endif
}

template<typename T>
inline T fast_acosh(T x) noexcept
{
#if defined(_MSC_VER) && !defined(__AVX2__)
    return std::log(x + std::sqrt(x * x - T(1)));
#else
    return std::log(x + std::sqrt(std::fma(x, x, -T(1))));
#endif
}

template<typename T>
class EasingCubicBezier
{
public:
    enum TYPE
    {
        NONE   = 0x0,
        P3     = 0x1,
        X2     = 0x2,
        X3P0   = 0x3,
        X3COS  = 0x4,
        X3SINH = 0x5,
        X3COSH = 0x6
    };

private:
public:
    TYPE mType;
    T mA;
    T mB;
    T mC;
    T mD;
    T mL;
    T mK;

public:
    static T DEFAULT_EPSILON;

    EasingCubicBezier() noexcept = default;

    EasingCubicBezier(TYPE type, T a, T b, T c, T d, T l, T k) noexcept
        : mType(type)
        , mA(a)
        , mB(b)
        , mC(c)
        , mD(d)
        , mL(l)
        , mK(k)
    {
    }

    EasingCubicBezier(const std::array<T, 4>& P_X, const std::array<T, 4>& P_Y) noexcept
    {
        makeEasingFromBezier(P_X, P_Y);
    }

    EasingCubicBezier(T x0, T y0, T x1, T y1, T x2, T y2, T x3, T y3) noexcept
    {
        makeEasingFromBezier({x0, x1, x2, x3}, {y0, y1, y2, y3});
    }
    
    T innerFunc(T x) const noexcept
    {
        return std::fma(mL, x, mK);
    }

    T evaluate(T x) const noexcept
    {
        constexpr T one_third = T(1) / T(3);
        constexpr T two_third_pi = -T(2) / T(3) * std::numbers::pi_v<T>;
#if defined(_MSC_VER) && !defined(__AVX2__)
        T phi = mL * x + mK;
#else
        T phi = std::fma(mL, x, mK);
#endif
        switch (mType)
        {
        case TYPE::P3:
            break;
        case TYPE::X2:
            phi = std::sqrt(std::max(T(0), phi));
            break;
        case TYPE::X3P0:
#ifdef FP_PRECISE
            phi = std::cbrt(phi);
#else
            phi = fast_cbrt(phi);
#endif
            break;
        case TYPE::X3COS:
#ifdef FP_PRECISE
            phi = std::cos(std::acos(std::clamp(phi, T(-1), T(1))) / T(3) + two_third_pi);
#else
#if defined(_MSC_VER) && !defined(__AVX2__)
            phi = std::cos(std::acos(std::clamp(phi, T(-1), T(1))) * one_third + two_third_pi);
#else
            phi = std::cos(std::fma(std::acos(std::clamp(phi, T(-1), T(1))), one_third, two_third_pi));
#endif
#endif
            break;
        case TYPE::X3SINH:
#ifdef FP_PRECISE
            phi = std::sinh(std::asinh(phi) / T(3));
#else
            phi = std::sinh(fast_asinh(phi) * one_third);
#endif
            break;
        case TYPE::X3COSH:
#ifdef FP_PRECISE
            phi = phi >= T(1) ? std::cosh(std::acosh(phi) / T(3))
                              : std::cos(std::acos(std::max(T(-1), phi)) / T(3));
#else
            phi = phi >= T(1) ? std::cosh(fast_acosh(phi) * one_third)
                              : std::cos(std::acos(std::max(T(-1), phi)) * one_third);
#endif
            break;
        default:
            UNREACHABLE();
            break;
        }
#if defined(_MSC_VER) && !defined(__AVX2__)
        return ((mA * phi + mB) * phi + mC) * phi + mD;
#else
        return std::fma(std::fma(std::fma(mA, phi, mB), phi, mC), phi, mD);
#endif
    }

    void makeEasingFromBezier(const std::array<T, 4>& P_X, const std::array<T, 4>& P_Y, T epsilon = DEFAULT_EPSILON) noexcept
    {
        std::array<T, 4> polyX = calculatePolynomial(P_X, P_X[0], P_X[3]);
        std::array<T, 4> polyY = calculatePolynomial(P_Y, P_X[0], P_X[3]);
        makeEasingFromPolynomial(polyX, polyY, epsilon);
    }

    void makeEasingFromPolynomial(const std::array<T, 4>& polyX, const std::array<T, 4>& polyY, T epsilon = DEFAULT_EPSILON) noexcept
    {
        const T b_a = polyX[1] / polyX[0];
        const T c_a = polyX[2] / polyX[0];
        const T d_a = polyX[3] / polyX[0];
        const T p = std::fma(b_a / T(-3), b_a, c_a);
        const T b_a2 = polyX[1] / (T(3) * polyX[0]);
        const T q = std::fma(std::fma(T(2) * b_a2, b_a2, -c_a), b_a2, d_a);
        const T signP = std::copysign(T(1), p);
        const T signB = std::copysign(T(1), polyX[1]);
        mType = calculateType(polyX, epsilon);
        T A = T(1), B = T(1), extra = T(1);
        switch (mType)
        {
        case TYPE::P3:
            A = T(1) / polyX[2];
            B = -(polyX[3] / polyX[2]);
            mL = T(1);
            mK = T(0);
            break;
        case TYPE::X2:
            A = signB;
            B = -(polyX[2] / (T(2) * polyX[1]));
            mL = T(1) / polyX[1];
            mK = std::fma(B, B, -polyX[3] / polyX[1]);
            break;
        case TYPE::X3P0:
            B = -b_a2;
            mL = T(1) / polyX[0];
            mK = -q;
            break;
        case TYPE::X3COSH:
            extra = signB;
            [[fallthrough]];
        case TYPE::X3COS:
            [[fallthrough]];
        case TYPE::X3SINH:
            A = T(-2) * signP * extra * std::sqrt(std::fabs(p) / T(3));
            B = -b_a2;
            mL = (T(3) * signP) / (p * polyX[0] * A);
            mK = (T(3) * signP) * (-q / p / A);
            break;
        default:
            UNREACHABLE();
            break;
        }
        mA = polyY[0] * A * A * A;
        mB = std::fma(T(3) * polyY[0], B, polyY[1]) * A * A;
        mC = std::fma(std::fma(T(3) * polyY[0], B, T(2) * polyY[1]), B, polyY[2]) * A;
        mD = std::fma(std::fma(std::fma(polyY[0], B, polyY[1]), B, polyY[2]), B, polyY[3]);
    }

    static TYPE calculateType(const std::array<T, 4>& polyX, T epsilon = DEFAULT_EPSILON) noexcept
    {
        const T b_a = polyX[1] / polyX[0];
        const T c_a = polyX[2] / polyX[0];
        const T d_a = polyX[3] / polyX[0];
        const T p = c_a - b_a * b_a / T(3);
        const T b_a2 = polyX[1] / (T(3) * polyX[0]);
        const T q = (T(2.0) * b_a2 * b_a2 - c_a) * b_a2 + d_a;
        TYPE type = TYPE::NONE;
        if (isZero(polyX[0], epsilon) && isZero(polyX[1], epsilon))
            type = TYPE::P3;
        else if (isZero(polyX[0], epsilon))
            type = TYPE::X2;
        else if (isZero(p, epsilon))
            type = TYPE::X3P0;
        else if (polyX[0] < T(0))
            type = TYPE::X3COS;
        else if (p > T(0))
            type = TYPE::X3SINH;
        else if (p < T(0))
            type = TYPE::X3COSH;
        else
            type = TYPE::NONE;
        return type;
    }

    static std::array<T, 4> rebasePointsTo01(const std::array<T, 4>& P) noexcept
    {
        const T diff = P[3] - P[0];
        return { T(0), (P[1] - P[0]) / diff, (P[2] - P[0]) / diff, T(1) };
    }

    static std::array<T, 4> calculatePolynomial(const std::array<T, 4>& P, T x0, T x3) noexcept
    {
        const std::array<T, 4> m0 = { T(-1), T(3), T(-3), T(1) };
        const std::array<T, 4> m1 = { T(3) * x3, T(-3) * std::fma(T(2), x3, x0), T(3) * std::fma(T(2), x0, x3), T(-3) * x0 };
        const std::array<T, 4> m2 = { T(-3) * x3 * x3, T(3) * x3 * std::fma(T(2), x0, x3), T(-3) * x0 * std::fma(T(2), x3, x0), T(3) * x0 * x0 };
        const std::array<T, 4> m3 = { x3 * x3 * x3, T(-3) * x0 * x3 * x3, T(3) * x0 * x0 * x3, -x0 * x0 * x0 };
        const T xRange = T(1) / ((x3 - x0) * (x3 - x0) * (x3 - x0));
        const std::array<T, 4> poly =
        {
            (m0[0] * P[0] + m0[1] * P[1] + m0[2] * P[2] + m0[3] * P[3]) * xRange,
            (m1[0] * P[0] + m1[1] * P[1] + m1[2] * P[2] + m1[3] * P[3]) * xRange,
            (m2[0] * P[0] + m2[1] * P[1] + m2[2] * P[2] + m2[3] * P[3]) * xRange,
            (m3[0] * P[0] + m3[1] * P[1] + m3[2] * P[2] + m3[3] * P[3]) * xRange
        };
        return poly;
    }

    static std::array<T, 4> calculatePolynomial01(const std::array<T, 4>& P) noexcept
    {
        const std::array<T, 4> m0 = { T(-1), T(3), T(-3), T(1) };
        const std::array<T, 4> m1 = { T(3), T(-6), T(3), T(0) };
        const std::array<T, 4> m2 = { T(-3), T(3), T(0), T(0) };
        const std::array<T, 4> m3 = { T(1), T(0), T(0), T(0) };
        const std::array<T, 4> poly =
        {
            (m0[0] * P[0] + m0[1] * P[1] + m0[2] * P[2] + m0[3] * P[3]),
            (m1[0] * P[0] + m1[1] * P[1] + m1[2] * P[2] + m1[3] * P[3]),
            (m2[0] * P[0] + m2[1] * P[1] + m2[2] * P[2] + m2[3] * P[3]),
            (m3[0] * P[0] + m3[1] * P[1] + m3[2] * P[2] + m3[3] * P[3])
        };
        return poly;
    }
};

template<typename T>
T EasingCubicBezier<T>::DEFAULT_EPSILON = std::numeric_limits<T>::epsilon() * T(10);

using EasingCubicBezierf = EasingCubicBezier<float>;
using EasingCubicBezierd = EasingCubicBezier<double>;
