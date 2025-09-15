#pragma once

#include <optional>
#include <array>
#include <cmath>
#include <cfloat>

namespace BlenderOrg
{
    inline double sqrt3d(double d)
    {
        if (d == 0.0)
        {
            return 0.0;
        }
        else if (d < 0.0)
        {
            return -std::exp(std::log(-d) / 3.0);
        }
        else
        {
            return std::exp(std::log(d) / 3.0);
        }
    }

    inline int solve_cubic(double c0, double c1, double c2, double c3, float* o)
    {
        const float SMALL = -1.0e-10f;
        const float END = 1.000001f;
        double a, b, c, p, q, d, t, phi;
        int nr = 0;
        if (c3 != 0.0)
        {
            a = c2 / c3;
            b = c1 / c3;
            c = c0 / c3;
            a = a / 3;
            p = b / 3 - a * a;
            q = (2 * a * a * a - a * b + c) / 2;
            d = q * q + p * p * p;
            if (d > 0.0)
            {
                t = std::sqrt(d);
                o[0] = (float)(sqrt3d(-q + t) + sqrt3d(-q - t) - a);
                if ((o[0] >= SMALL) && (o[0] <= END))
                    return 1;
                return 0;
            }
            if (d == 0.0)
            {
                t = sqrt3d(-q);
                o[0] = (float)(2 * t - a);
                if ((o[0] >= SMALL) && (o[0] <= END))
                    nr++;
                o[nr] = (float)(-t - a);
                if ((o[nr] >= SMALL) && (o[nr] <= END))
                    return nr + 1;
                return nr;
            }
            phi = std::acos(-q / std::sqrt(-(p * p * p)));
            t = std::sqrt(-p);
            p = std::cos(phi / 3);
            q = std::sqrt(3 - 3 * p * p);
            o[0] = (float)(2 * t * p - a);
            if ((o[0] >= SMALL) && (o[0] <= END))
                nr++;
            o[nr] = (float)(-t * (p + q) - a);
            if ((o[nr] >= SMALL) && (o[nr] <= END))
                nr++;
            o[nr] = (float)(-t * (p - q) - a);
            if ((o[nr] >= SMALL) && (o[nr] <= END))
                return nr + 1;
            return nr;
        }
        a = c2;
        b = c1;
        c = c0;
        if (a != 0.0)
        {
            /* Discriminant */
            p = b * b - 4 * a * c;
            if (p > 0) {
                p = std::sqrt(p);
                o[0] = (float)((-b - p) / (2 * a));
                if ((o[0] >= SMALL) && (o[0] <= END))
                    nr++;
                o[nr] = (float)((-b + p) / (2 * a));
                if ((o[nr] >= SMALL) && (o[nr] <= END))
                    return nr + 1;
                return nr;
            }
            if (p == 0)
            {
                o[0] = (float)(-b / (2 * a));
                if ((o[0] >= SMALL) && (o[0] <= END))
                    return 1;
            }
            return 0;
        }
        if (b != 0.0)
        {
            o[0] = (float)(-c / b);
            if ((o[0] >= SMALL) && (o[0] <= END))
                return 1;
            return 0;
        }
        if (c == 0.0)
        {
            o[0] = 0.0f;
            return 1;
        }
        return 0;
    }
    inline int findzero(float x, float q0, float q1, float q2, float q3, float* o)
    {
        const double c0 = q0 - x;
        const double c1 = 3.0 * (q1 - q0);
        const double c2 = 3.0 * (q0 - 2.0 * q1 + q2);
        const double c3 = q3 - q0 + 3.0 * (q1 - q2);
        return solve_cubic(c0, c1, c2, c3, o);
    }
    inline void berekeny(float f1, float f2, float f3, float f4, float* o, int b)
    {
        float c0 = f1;
        float c1 = 3.0f * (f2 - f1);
        float c2 = 3.0f * (f1 - 2.0f * f2 + f3);
        float c3 = f4 - f1 + 3.0f * (f2 - f3);
        for (int a = 0; a < b; a++)
        {
            float t = o[a];
            o[a] = c0 + t * c1 + t * t * c2 + t * t * t * c3;
        }
    }
    inline void BKE_fcurve_correct_bezpart(const float v1[2], float v2[2], float v3[2], const float v4[2])
    {
        float h1[2], h2[2], len1, len2, len, fac;

        /* Calculate handle deltas. */
        h1[0] = v1[0] - v2[0];
        h1[1] = v1[1] - v2[1];

        h2[0] = v4[0] - v3[0];
        h2[1] = v4[1] - v3[1];

        /* Calculate distances:
         * - len  = Span of time between keyframes.
         * - len1 = Length of handle of start key.
         * - len2 = Length of handle of end key.
         */
        len = v4[0] - v1[0];
        len1 = fabsf(h1[0]);
        len2 = fabsf(h2[0]);

        /* If the handles have no length, no need to do any corrections. */
        if ((len1 + len2) == 0.0f) {
            return;
        }

        /* To prevent looping or rewinding, handles cannot
         * exceed the adjacent key-frames time position. */
        if (len1 > len) {
            fac = len / len1;
            v2[0] = (v1[0] - fac * h1[0]);
            v2[1] = (v1[1] - fac * h1[1]);
        }

        if (len2 > len) {
            fac = len / len2;
            v3[0] = (v4[0] - fac * h2[0]);
            v3[1] = (v4[1] - fac * h2[1]);
        }
    }

    inline float fcurve_eval_keyframes_interpolate(std::array<float, 4> P_X, std::array<float, 4> P_Y, float evaltime)
    {
        /* Bezier interpolation. */
        /* (v1, v2) are the first keyframe and its 2nd handle. */
        float v1[2], v2[2], v3[2], v4[2], opl[32];
        v1[0] = P_X[0];
        v1[1] = P_Y[0];
        v2[0] = P_X[1];
        v2[1] = P_Y[1];
        /* (v3, v4) are the last keyframe's 1st handle + the last keyframe. */
        v3[0] = P_X[2];
        v3[1] = P_Y[2];
        v4[0] = P_X[3];
        v4[1] = P_Y[3];

        if (fabsf(v1[1] - v4[1]) < FLT_EPSILON && fabsf(v2[1] - v3[1]) < FLT_EPSILON &&
            fabsf(v3[1] - v4[1]) < FLT_EPSILON) {
            /* Optimization: If all the handles are flat/at the same values,
             * the value is simply the shared value (see T40372 -> F91346).
             */
            return v1[1];
        }
        /* Adjust handles so that they don't overlap (forming a loop). */
        BKE_fcurve_correct_bezpart(v1, v2, v3, v4);
        /* Try to get a value for this position - if failure, try another set of points. */
        if (!findzero(evaltime, v1[0], v2[0], v3[0], v4[0], opl))
        {
            return 0.0f;
        }
        berekeny(v1[1], v2[1], v3[1], v4[1], opl, 1);
        return opl[0];
    }
}

namespace BlenderOptim
{
    template<typename T>
    inline std::optional<T> solve_cubic(T c0, T c1, T c2, T c3)
    {
        const T ZERO = -std::numeric_limits<T>::min();
        const T ONE = std::nextafterf(T(1), T(2));
        T roots[3] =
        { 
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max() 
        };
        if (std::fabs(c3) >= std::numeric_limits<T>::epsilon())
        {
            T a = c2 / c3 / T(3);
            T b = c1 / c3;
            T c = c0 / c3;
            T p = b / T(3) - a * a;
            T q = (T(2) * a * a * a - a * b + c) / T(2);
            T d = q * q + p * p * p;
            if (std::fabs(d) < std::numeric_limits<T>::epsilon())
            {
                T t = std::cbrt(-q);
                roots[0] = T(2) * t - a;
                roots[1] = -t - a;
            }
            else if (d > T(0))
            {
                T t = std::sqrt(d);
                roots[0] = std::cbrt(-q + t) + std::cbrt(-q - t) - a;
            }
            else
            {
                T phi = std::acos(-q / std::sqrt(-(p * p * p)));
                T t = std::sqrt(-p);
                p = std::cos(phi / T(3));
                q = std::sqrt(T(3) - T(3) * p * p);
                roots[0] = T(2) * t * p - a;
                roots[1] = -t * (p + q) - a;
                roots[2] = -t * (p - q) - a;
            }
        }
        else
        {
            T a = c2;
            T b = c1;
            T c = c0;
            if (std::fabs(a) >= std::numeric_limits<T>::epsilon())
            {
                T p = b * b - T(4) * a * c;
                if (std::fabs(p) < std::numeric_limits<T>::epsilon())
                {
                    roots[0] = -b / (T(2) * a);
                }
                else if (p > T(0))
                {
                    p = std::sqrt(p);
                    roots[0] = (-b - p) / (T(2) * a);
                    roots[1] = (-b + p) / (T(2) * a);
                }
            }
            else if (std::fabs(b) >= std::numeric_limits<T>::epsilon())
            {
                roots[0] = -c / b;
            }
            else if (std::fabs(c) < std::numeric_limits<T>::epsilon())
            {
                roots[0] = T(0);
            }
        }
        if (roots[0] >= ZERO && roots[0] <= ONE)
            return std::make_optional(roots[0]);
        if (roots[1] >= ZERO && roots[1] <= ONE)
            return std::make_optional(roots[1]);
        if (roots[2] >= ZERO && roots[2] <= ONE)
            return std::make_optional(roots[2]);
        return {};
    }
    template<typename T>
    inline std::optional<T> findzero(T x, T q0, T q1, T q2, T q3)
    {
        const T c0 = q0 - x;
        const T c1 = T(3) * (q1 - q0);
        const T c2 = T(3) * (q0 - T(2) * q1 + q2);
        const T c3 = q3 - q0 + T(3) * (q1 - q2);
        return solve_cubic(c0, c1, c2, c3);
    }
    template<typename T>
    inline T berekeny(T f1, T f2, T f3, T f4, T t)
    {
        T c0 = f1;
        T c1 = T(3) * (f2 - f1);
        T c2 = T(3) * (f1 - T(2) * f2 + f3);
        T c3 = f4 - f1 + T(3) * (f2 - f3);
        return std::fma(std::fma(std::fma(c3, t, c2), t, c1), t, c0);
    }
    template<typename T>
    inline T fcurve_eval_keyframes_interpolate(std::array<T, 4> P_X, std::array<T, 4> P_Y, T t)
    {
        auto found = findzero(t, P_X[0], P_X[1], P_X[2], P_X[3]);
        return found.has_value() ? berekeny(P_Y[0], P_Y[1], P_Y[2], P_Y[3], found.value()) : T(0);
    }
}

namespace BlenderNumeric
{
    template<typename T>
    inline T cubic_f(T x, T a, T b, T c, T d)
    {
#if defined(_MSC_VER) && defined(__AVX2__)
        return std::fma(std::fma(std::fma(a, x, b), x, c), x, d);
#else
        return ((a * x + b) * x + c) * x + d;
#endif
    }
    template<typename T>
    inline T cubic_df(T x, T a, T b, T c)
    {
#if defined(_MSC_VER) && defined(__AVX2__)
        return std::fma(std::fma(T(3) * a, x, T(2) * b), x, c);
#else
        return (T(3) * a * x + T(2) * b) * x + c;
#endif
    }
    template<typename T>
    inline T cubic_ddf(T x, T a, T b)
    {
#if defined(_MSC_VER) && defined(__AVX2__)
        return std::fma(T(6) * a, x, T(2) * b);
#else
        return T(6) * a * x + T(2) * b;
#endif
    }

    template<typename T>
    T findRootNewtonRaphson(T a, T b, T c, T d, T x, int maxIter = 100)
    {
        while (maxIter-- > 0)
        {
            const T fx = cubic_f(x, a, b, c, d);
            const T fpx = cubic_df(x, a, b, c);
            const T x_new = x - fx / fpx;
            if (std::abs(x_new - x) < std::numeric_limits<T>::epsilon())
                return x_new;
            x = x_new;
        }
        return x;
    }
    template<typename T>
    T findRootHalley(T a, T b, T c, T d, T x, int maxIter = 100)
    {
        for (int i = 0; i < maxIter; ++i)
        {
            const T fx = cubic_f(x, a, b, c, d);
            const T fpx = cubic_df(x, a, b, c);
            const T fppx = cubic_ddf(x, a, b);
            const T denominator = fpx * fpx - T(0.5) * fx * fppx;
            if (std::fabs(denominator) < std::numeric_limits<T>::epsilon())
                break;
            const T x_new = x - (fx * fpx) / denominator;
            if (std::abs(x_new - x) < std::numeric_limits<T>::epsilon())
                return x_new;
            x = x_new;
        }
        return x;
    }
    template<typename T>
    T findRootLaguerre(T a, T b, T c, T d, T x, int maxIter = 100)
    {
        while (maxIter-- > 0)
        {
            T fx = cubic_f(x, a, b, c, d);
            if (std::fabs(fx) < std::numeric_limits<T>::epsilon())
                break;
            T G = cubic_df(x, a, b, c) / fx;
            T H = G * G - cubic_ddf(x, a, b) / fx;
            T sqrt_term = std::sqrt( T(2) * (T(3) * H - G * G));
            T denom1 = G + sqrt_term;
            T denom2 = G - sqrt_term;
            T denom = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;
            T a_n = T(3) / denom;
            x = x - a_n;
            if (std::fabs(a_n) < std::numeric_limits<T>::epsilon())
                break;
        }
        return x;
    }
    template<typename T>
    T findRootBisection(T a, T b, T c, T d, T left, T right, int maxIter = 100)
    {
        if (cubic_f(left, a, b, c, d) * cubic_f(right, a, b, c, d) >= T(0))
            return std::numeric_limits<T>::quiet_NaN();
        while ((right - left) >= std::numeric_limits<T>::epsilon() && --maxIter > 0)
        {
            T mid = (left + right) * T(0.5);
            T f_mid = cubic_f(mid, a, b, c, d);
            if (std::fabs(f_mid) < std::numeric_limits<T>::epsilon())
                return mid;
            if (cubic_f(left, a, b, c, d) * f_mid < T(0))
                right = mid;
            else
                left = mid;
        }
        return (left + right) * T(0.5);
    }
    template<typename T>
    T findRootBrents(T a, T b, T c, T d, T l, T r, int maxIter = 100)
    {
        T f_l = cubic_f(l, a, b, c, d);
        T f_r = cubic_f(r, a, b, c, d);
        if (f_l * f_r >= T(0))
            return std::numeric_limits<double>::quiet_NaN();
        if (std::fabs(f_l) < std::fabs(f_r))
        {
            std::swap(l, r);
            std::swap(f_l, f_r);
        }
        double new_l = l;
        double f_c = f_l;
        double s = r;
        double new_r = r;
        bool mflag = true;
        for (int iter = 0; iter < maxIter; ++iter)
        {
            if (f_l != f_c && f_r != f_c) // Inverse quadratic interpolation
                s = (l * f_r * f_c) / ((f_l - f_r) * (f_l - f_c)) +
                    (r * f_l * f_c) / ((f_r - f_l) * (f_r - f_c)) +
                    (new_l * f_l * f_r) / ((f_c - f_l) * (f_c - f_r));
            else // Secant method                
                s = r - f_r * (r - l) / (f_r - f_l);
            bool condition1 = (s < (T(3) * l + r) / T(4) || s > r);
            bool condition2 = mflag && (std::fabs(s - r) >= std::fabs(r - new_l) / T(2));
            bool condition3 = !mflag && (std::fabs(s - r) >= std::fabs(new_l - new_r) / T(2));
            bool condition4 = mflag && (std::fabs(r - new_l) < std::numeric_limits<T>::epsilon());
            bool condition5 = !mflag && (std::fabs(new_l - new_r) < std::numeric_limits<T>::epsilon());
            mflag = condition1 || condition2 || condition3 || condition4 || condition5;
            if (mflag)
                s = (l + r) / T(2);
            T f_s = cubic_f(s, a, b, c, d);
            new_r = new_l;
            new_l = r;
            f_c = f_r;
            if (f_l * f_s < T(0))
            {
                r = s;
                f_r = f_s;
            }
            else
            {
                l = s;
                f_l = f_s;
            }
            if (std::fabs(f_l) < std::fabs(f_r))
            {
                std::swap(l, r);
                std::swap(f_l, f_r);
            }
            if (std::fabs(r - l) < std::numeric_limits<T>::epsilon() || std::fabs(f_s) < std::numeric_limits<T>::epsilon())
                return s;
        }
        return s;
    }

    template<typename T>
    inline std::optional<T> findzero(T x, T q0, T q1, T q2, T q3)
    {
        const T c0 = q0 - x;
        const T c1 = T(3) * (q1 - q0);
        const T c2 = T(3) * (q0 - T(2) * q1 + q2);
        const T c3 = q3 - q0 + T(3) * (q1 - q2);
        //return findRootNewtonRaphson(c3, c2, c1, c0, x);
        return findRootHalley(c3, c2, c1, c0, x);
        //return findRootLaguerre(c3, c2, c1, c0, T(0.5));
        //return findRootBisection(c3, c2, c1, c0, T(0.5));
        //return findRootBrents(c3, c2, c1, c0, T(0.5));
    }
    template<typename T>
    inline T berekeny(T f1, T f2, T f3, T f4, T t)
    {
        T c0 = f1;
        T c1 = T(3) * (f2 - f1);
        T c2 = T(3) * (f1 - T(2) * f2 + f3);
        T c3 = f4 - f1 + T(3) * (f2 - f3);
        return std::fma(std::fma(std::fma(c3, t, c2), t, c1), t, c0);
    }
    template<typename T>
    inline T fcurve_eval_keyframes_interpolate1(std::array<T, 4> P_X, std::array<T, 4> P_Y, T t)
    {
        const T q0 = P_X[0];
        const T q1 = P_X[1];
        const T q2 = P_X[2];
        const T q3 = P_X[3];
        const T c0 = q0 - t;
        const T c1 = T(3) * (q1 - q0);
        const T c2 = T(3) * (q0 - T(2) * q1 + q2);
        const T c3 = q3 - q0 + T(3) * (q1 - q2);
        auto found = findRootNewtonRaphson(c3, c2, c1, c0, t);
        //auto found = findRootHalley(c3, c2, c1, c0, t);
        //auto found = findRootLaguerre(c3, c2, c1, c0, t);
        //auto found = findRootBisection(c3, c2, c1, c0, t);
        //auto found = findRootBrents(c3, c2, c1, c0, t);
        return berekeny(P_Y[0], P_Y[1], P_Y[2], P_Y[3], found);
    }
    template<typename T>
    inline T fcurve_eval_keyframes_interpolate2(std::array<T, 4> P_X, std::array<T, 4> P_Y, T t)
    {
        const T q0 = P_X[0];
        const T q1 = P_X[1];
        const T q2 = P_X[2];
        const T q3 = P_X[3];
        const T c0 = q0 - t;
        const T c1 = T(3) * (q1 - q0);
        const T c2 = T(3) * (q0 - T(2) * q1 + q2);
        const T c3 = q3 - q0 + T(3) * (q1 - q2);
        auto found = findRootHalley(c3, c2, c1, c0, t);
        return berekeny(P_Y[0], P_Y[1], P_Y[2], P_Y[3], found);
    }
    template<typename T>
    inline T fcurve_eval_keyframes_interpolate3(std::array<T, 4> P_X, std::array<T, 4> P_Y, T t)
    {
        const T q0 = P_X[0];
        const T q1 = P_X[1];
        const T q2 = P_X[2];
        const T q3 = P_X[3];
        const T c0 = q0 - t;
        const T c1 = T(3) * (q1 - q0);
        const T c2 = T(3) * (q0 - T(2) * q1 + q2);
        const T c3 = q3 - q0 + T(3) * (q1 - q2);
        auto found = findRootBisection(c3, c2, c1, c0, T(0), T(1));
        return berekeny(P_Y[0], P_Y[1], P_Y[2], P_Y[3], found);
    }
}
