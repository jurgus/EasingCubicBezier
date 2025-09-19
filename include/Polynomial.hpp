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
#include <span>
#include <cstdint>
#include <cassert>

template<typename T>
class Polynomial
{
    uint32_t m_degree;
    std::array<T, 6> m_coefficients;
public:
    explicit Polynomial(uint32_t d = 0u) noexcept
        : m_coefficients{}
        , m_degree(d)
    {
    }
    Polynomial(const Polynomial&) noexcept = default;
    Polynomial& operator=(const Polynomial&) noexcept = default;

    Polynomial(std::initializer_list<T> coefficients) noexcept
        : m_degree(static_cast<uint32_t>(coefficients.size() - 1))
        , m_coefficients{}
    {
        assert(coefficients.size() >= 2 && coefficients.size() <= m_coefficients.size());
        auto it = coefficients.begin();
        for (uint32_t i = m_degree + 1; i-- > 0; )
            m_coefficients[i] = *it++;
    }
#if __cplusplus >= 202002L
    template<typename... Coefs>
        requires (sizeof...(Coefs) >= 2 && sizeof...(Coefs) <= 6)
#else
    template<typename... Coefs,
             typename = typename std::enable_if<(sizeof...(Coefs) >= 2) && (sizeof...(Coefs) <= 6)>::type>
#endif
        Polynomial(Coefs&&... coefficients) noexcept
        : m_degree(static_cast<uint32_t>(sizeof...(Coefs) - 1))
        , m_coefficients{}

    {
        T arr[] = { static_cast<T>(std::forward<Coefs>(coefficients))... };
        for (size_t i = 0; i < sizeof...(Coefs); ++i)
            m_coefficients[m_degree - i] = arr[i];
    }
    template<size_t N>
    Polynomial(const std::array<T, N>& arr) noexcept
        : m_degree(static_cast<uint32_t>(arr.size()))
        , m_coefficients{}
    {
        static_assert(N <= 6);
        for (size_t i = 0; i < arr.size(); ++i)
            m_coefficients[m_degree - i] = arr[i];
    }
    Polynomial(const T* ptr, size_t size, uint32_t degree) noexcept
        : m_degree(degree)
        , m_coefficients{}
    {
        assert(degree <= 6 && size <= degree);
        for (size_t i = 0; i < size; ++i)
            m_coefficients[m_degree - i] = ptr[i];
    }
    Polynomial(std::span<T> span, uint32_t degree) noexcept
        : m_degree(degree)
        , m_coefficients{}
    {
        assert(degree <= 6 && span.size() <= degree);
        for (size_t i = 0; i < span.size(); ++i)
            m_coefficients[m_degree - i] = span[i];
    }
#if __cplusplus >= 202002L
    template<typename... Coefs>
    requires (sizeof...(Coefs) >= 2 && sizeof...(Coefs) <= 6)
#else
    template<typename... Coefs,
             typename = typename std::enable_if<(sizeof...(Coefs) >= 2) && (sizeof...(Coefs) <= 6)>::type>
#endif
    void set(Coefs&&... coefficients) noexcept
    {
        m_degree = static_cast<uint32_t>(sizeof...(Coefs) - 1);
        m_coefficients = {};
        T tmp[] = { static_cast<T>(std::forward<Coefs>(coefficients))... };
        for (std::size_t i = 0; i < sizeof...(Coefs); ++i)
            m_coefficients[m_degree - i] = tmp[i];
    }
    template<size_t N>
    void set(const std::array<T, N>& arr) noexcept
    {
        assert(N <= 6);
        m_degree = static_cast<uint32_t>(arr.size());
        m_coefficients = {};
        for (size_t i = 0; i < arr.size(); ++i)
            m_coefficients[m_degree - i] = arr[i];
    }
    void set(const T* ptr, size_t size, uint32_t degree) noexcept
    {
        assert(degree <= 6 && size <= degree);
        m_degree = degree;
        m_coefficients = {};
        for (size_t i = 0; i <= size; ++i)
            m_coefficients[m_degree - i] = ptr[i];
    }
    void set(std::span<T> span, uint32_t degree) noexcept
    {
        assert(degree <= 6 && span.size() <= degree);
        for (size_t i = 0; i < span.size(); ++i)
            m_coefficients[m_degree - i] = span[i];
    }

    constexpr T* ptr() noexcept
    {
        return m_coefficients;
    }
    constexpr const T* ptr() const noexcept
    {
        return m_coefficients;
    }
    constexpr const T& operator[](size_t i) const noexcept
    {
        assert(i <= m_degree);
        return m_coefficients[i];
    }
    constexpr T& operator[](size_t i) noexcept
    {
        assert(i <= m_degree);
        return m_coefficients[i];
    }
    constexpr const T& a() const noexcept
    {
        return m_coefficients[m_degree];
    }
    constexpr T& a() noexcept
    {
        return m_coefficients[m_degree];
    }
    constexpr const T& b() const noexcept
    {
        return m_coefficients[m_degree - 1];
    }
    constexpr T& b() noexcept
    {
        return m_coefficients[m_degree - 1];
    }
    constexpr const T& c() const noexcept
    {
        return m_coefficients[m_degree - 2];
    }
    constexpr T& c() noexcept
    {
        return m_coefficients[m_degree - 2];
    }
    constexpr const T& d() const noexcept
    {
        return m_coefficients[m_degree - 3];
    }
    constexpr T& d() noexcept
    {
        return m_coefficients[m_degree - 3];
    }
    constexpr const T& e() const noexcept
    {
        return m_coefficients[m_degree - 4];
    }
    constexpr T& e() noexcept
    {
        return m_coefficients[m_degree - 4];
    }
    constexpr const T& f() const noexcept
    {
        return m_coefficients[m_degree - 5];
    }
    constexpr T& f() noexcept
    {
        return m_coefficients[m_degree - 5];
    }

    constexpr T operator()(T t) const noexcept
    {
        return evaluate(t);
    }
    constexpr T evaluate(T t) const noexcept
    {
        T ret = m_coefficients[m_degree];
        for (uint32_t i = m_degree; i > 0u; --i)
            ret = std::fma(ret, t, m_coefficients[i - 1u]);
        return ret;
    }
    constexpr T evaluateFirstDerivative(T t) const noexcept
    {
        T ret = T(0);
        for (uint32_t i = m_degree; i > 1; --i)
            ret += m_coefficients[i] * std::pow(t, i - 1) * i;
        ret += m_coefficients[1];
        return ret;
    }
    constexpr T evaluateSecondDerivative(T t) const noexcept
    {
        T ret = T(0);
        for (uint32_t i = m_degree; i > 2; --i)
            ret += m_coefficients[i] * std::pow(t, i - 2) * i * T(i - 1);
        ret += m_coefficients[2] * T(2);
        return ret;
    }

    constexpr T evaluateIntegral(T t) const noexcept
    {
        T ret = T(0);
        for (int32_t i = m_degree; i > 0; --i)
            ret += m_coefficients[i] * std::pow(t, i + 1) / T(i + 1);
        ret += m_coefficients[0] * t;
        return ret;
    }
};

using Polynomialf = Polynomial<float>;
using Polynomiald = Polynomial<double>;
