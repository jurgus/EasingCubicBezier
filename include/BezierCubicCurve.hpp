/*
 * MIT License
 *
 * Copyright (c) 2025 £ukasz Izdebski
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <array>
#include "Polynomial.hpp"

template<typename T>
class BezierCubicCurve
{
    std::array<T, 4> x;
    std::array<T, 4> y;

public:
    BezierCubicCurve() noexcept = default;

    BezierCubicCurve(T x0, T y0, T x1, T y1, T x2, T y2, T x3, T y3) noexcept
    {
        x = { x0, x1, x2, x3 };
        y = { y0, y1, y2, y3 };
    }

    BezierCubicCurve(const std::array<T, 4>& _x, const std::array<T, 4>& _y) noexcept
    {
        x = _x;
        y = _y;
    }

    T getValX(T t) const noexcept
    {
        const T one = T(1) - t;
        return x[0] * one * one * one + T(3) * x[1] * one * one * t + T(3) * x[2] * one * t * t + x[3] * t * t * t;
    }

    T getValY(T t) const noexcept
    {
        const T one = T(1) - t;
        return y[0] * one * one * one + T(3) * y[1] * one * one * t + T(3) * y[2] * one * t * t + y[3] * t * t * t;
    }

    void setXPolynomial(const std::array<T, 4>& poly, T x0, T x3) noexcept
    {
        const T m0[4] = { x0 * x0 * x0, x0 * x0, x0, T(1) };
        const T m1[4] = { x0 * x0 * x3, x0 * (x0 + T(2) * x3) / T(3), (T(2) * x0 + x3) / T(3), T(1) };
        const T m2[4] = { x0 * x3 * x3, x3 * (T(2) * x0 + x3) / T(3), (x0 + T(2) * x3) / T(3), T(1) };
        const T m3[4] = { x3 * x3 * x3, x3 * x3, x3, T(1) };
        x[0] = (m0[0] * poly[3] + m0[1] * poly[2] + m0[2] * poly[1] + m0[3] * poly[0]);
        x[1] = (m1[0] * poly[3] + m1[1] * poly[2] + m1[2] * poly[1] + m1[3] * poly[0]);
        x[2] = (m2[0] * poly[3] + m2[1] * poly[2] + m2[2] * poly[1] + m2[3] * poly[0]);
        x[3] = (m3[0] * poly[3] + m3[1] * poly[2] + m3[2] * poly[1] + m3[3] * poly[0]);
    }

    void setYPolynomial(const std::array<T, 4>& poly, T x0, T x3) noexcept
    {
        const T m0[4] = { x0 * x0 * x0, x0 * x0, x0, T(1) };
        const T m1[4] = { x0 * x0 * x3, x0 * (x0 + T(2) * x3) / T(3), (T(2) * x0 + x3) / T(3), T(1) };
        const T m2[4] = { x0 * x3 * x3, x3 * (T(2) * x0 + x3) / T(3), (x0 + T(2) * x3) / T(3), T(1) };
        const T m3[4] = { x3 * x3 * x3, x3 * x3, x3, T(1) };
        y[0] = (m0[0] * poly[3] + m0[1] * poly[2] + m0[2] * poly[1] + m0[3] * poly[0]);
        y[1] = (m1[0] * poly[3] + m1[1] * poly[2] + m1[2] * poly[1] + m1[3] * poly[0]);
        y[2] = (m2[0] * poly[3] + m2[1] * poly[2] + m2[2] * poly[1] + m2[3] * poly[0]);
        y[3] = (m3[0] * poly[3] + m3[1] * poly[2] + m3[2] * poly[1] + m3[3] * poly[0]);
    }

    Polynomial<T> getXPolynomial() const noexcept
    {
        const T x0 = x[0];
        const T x3 = x[3];
        const T m0[4] = { T(-1), T(3), T(-3), T(1) };
        const T m1[4] = { T(3) * x3, T(-3) * (x0 + T(2) * x3), T(3) * (T(2) * x0 + x3), T(-3) * x0 };
        const T m2[4] = { T(-3) * x3 * x3, T(3) * x3 * (T(2) * x0 + x3), T(-3) * x0 * (x0 + T(2) * x3), T(3) * x0 * x0 };
        const T m3[4] = { x3 * x3 * x3, T(-3) * x0 * x3 * x3, T(3) * x0 * x0 * x3, -x0 * x0 * x0 };
        const T d = (x3 - x0) * (x3 - x0) * (x3 - x0);
        return
        {
            (m0[0] * x[0] + m0[1] * x[1] + m0[2] * x[2] + m0[3] * x[3]) / d,
            (m1[0] * x[0] + m1[1] * x[1] + m1[2] * x[2] + m1[3] * x[3]) / d,
            (m2[0] * x[0] + m2[1] * x[1] + m2[2] * x[2] + m2[3] * x[3]) / d,
            (m3[0] * x[0] + m3[1] * x[1] + m3[2] * x[2] + m3[3] * x[3]) / d
        };
    }

    Polynomial<T> getYPolynomial() const noexcept
    {
        const T x0 = x[0];
        const T x3 = x[3];
        const T m0[4] = { T(-1), T(3), T(-3), T(1) };
        const T m1[4] = { T(3) * x3, T(-3) * (x0 + T(2) * x3), T(3) * (T(2) * x0 + x3), T(-3) * x0 };
        const T m2[4] = { T(-3) * x3 * x3, T(3) * x3 * (T(2) * x0 + x3), T(-3) * x0 * (x0 + T(2) * x3), T(3) * x0 * x0 };
        const T m3[4] = { x3 * x3 * x3, T(-3) * x0 * x3 * x3, T(3) * x0 * x0 * x3, -x0 * x0 * x0 };
        const T d = (x3 - x0) * (x3 - x0) * (x3 - x0);
        return
        {
            (m0[0] * y[0] + m0[1] * y[1] + m0[2] * y[2] + m0[3] * y[3]) / d,
            (m1[0] * y[0] + m1[1] * y[1] + m1[2] * y[2] + m1[3] * y[3]) / d,
            (m2[0] * y[0] + m2[1] * y[1] + m2[2] * y[2] + m2[3] * y[3]) / d,
            (m3[0] * y[0] + m3[1] * y[1] + m3[2] * y[2] + m3[3] * y[3]) / d
        };
    }

    void split(const T& delta, BezierCubicCurve& lower, BezierCubicCurve& upper) noexcept
    {
    }

    BezierCubicCurve splitLower(const T& delta) noexcept
    {
        return BezierCubicCurve();
    }

    BezierCubicCurve splitUpper(const T& delta) noexcept
    {
        return BezierCubicCurve();
    }

    T valueX(const T& delta) noexcept
    {
        T R0, R1, R2, Q0, Q1, S0;
        calcX(R0, R1, R2, Q0, Q1, S0, delta);
        return S0;
    }

    T valueY(const T& delta) noexcept
    {
        T R0, R1, R2, Q0, Q1, S0;
        calcY(R0, R1, R2, Q0, Q1, S0, delta);
        return S0;
    }

    void calcX(T& R0, T& R1, T& R2, T& Q0, T& Q1, T& S0, const T& delta) noexcept
    {
        R0 = (x[1] - x[0]) * delta + x[0];
        R1 = (x[2] - x[1]) * delta + x[1];
        R2 = (x[3] - x[2]) * delta + x[2];
        Q0 = (R1 - R0) * delta + R0;
        Q1 = (R2 - R1) * delta + R1;
        S0 = (Q1 - Q0) * delta + Q0;
    }

    void calcY(T& R0, T& R1, T& R2, T& Q0, T& Q1, T& S0, const T& delta) noexcept
    {
        R0 = (y[1] - y[0]) * delta + y[0];
        R1 = (y[2] - y[1]) * delta + y[1];
        R2 = (y[3] - y[2]) * delta + y[2];
        Q0 = (R1 - R0) * delta + R0;
        Q1 = (R2 - R1) * delta + R1;
        S0 = (Q1 - Q0) * delta + Q0;
    }

    bool isValid() const noexcept
    {
        if (x[0] != std::numeric_limits<T>::max() || x[1] != std::numeric_limits<T>::max() ||
            x[2] != std::numeric_limits<T>::max() || x[3] != std::numeric_limits<T>::max() ||
            y[0] != std::numeric_limits<T>::max() || y[1] != std::numeric_limits<T>::max() ||
            y[2] != std::numeric_limits<T>::max() || y[3] != std::numeric_limits<T>::max())
        {
            return true;
        }
        return false;
    }
};

using BezierCubicCurvef = BezierCubicCurve<float>;
using BezierCubicCurved = BezierCubicCurve<double>;
