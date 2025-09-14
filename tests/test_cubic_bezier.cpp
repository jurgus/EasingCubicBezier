#include <easingCubicBezier.hpp>
#include <BezierCubicCurve.hpp>
#include <Polynomial.hpp>
#include <iostream>
#include <iomanip>
#include <utility>
#include "cubic_bezier.hpp"

#define _PI_ std::numbers::pi_v<float>

template<typename T>
static EasingCubicBezier<T> make_curve(const Polynomial<T>& pX, const Polynomial<T>& pY)
{
    EasingCubicBezier<T> ecb{};
    T b_a = pX.b() / pX.a();
    T c_a = pX.c() / pX.a();
    T d_a = pX.d() / pX.a();
    T p = c_a - b_a * b_a / T(3.0);
    T q = T(2.0 / 27.0) * b_a * b_a * b_a - T(1.0 / 3.0) * b_a * c_a + d_a;
    T s = T(4.0) * p * p * p + T(27.0) * q * q;
    typename EasingCubicBezier<T>::TYPE mType = EasingCubicBezier<T>::NONE;
    T alpha = T(0), beta = T(0), gamma = T(0), delta = T(0), l = T(0), k = T(0);

    if (isZero(pX.a()) && isZero(pX.b()))
    {
        mType = EasingCubicBezier<T>::P3;
        ecb = EasingCubicBezier<T>( mType, alpha, beta, gamma, delta, l, k );
    }
    else if (isZero(pX.a()))
    {
        mType = EasingCubicBezier<T>::X2;
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else if (isZero(p))
    {
        mType = EasingCubicBezier<T>::X3P0;
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else if (p > T(0.0))
    {
        mType = EasingCubicBezier<T>::X3SINH;
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else if (pX.a() < T(0.0))
    {
        const Polynomial<T>& xpoly = pX;
        const Polynomial<T>& ypoly = pY;
        T _u = T(1);
        mType = EasingCubicBezier<T>::X3COS;
        T b_a = xpoly.b() / xpoly.a();
        T c_a = xpoly.c() / xpoly.a();
        T d_a = xpoly.d() / xpoly.a();
        T p = c_a - b_a * b_a / T(3.0);
        T q = T(2.0 / 27.0) * b_a * b_a * b_a - T(1.0 / 3.0) * b_a * c_a + d_a;
        T A = T(2.0) * std::sqrt(-p / T(3.0));
        T B = -xpoly.b() / (T(3.0) * xpoly.a());
        T n = T(3.0) / (A * p);
        k = q;
        l = T(-1.0) / xpoly.a();
        k *= n;
        l *= n;
        T u = _u * T(2.0) * _PI_ / T(3.0);
        alpha = ypoly.a() * A * A * A;
        beta = T(3.0) * ypoly.a() * B * A * A + ypoly.b() * A * A;
        gamma = T(3.0) * ypoly.a() * B * B * A + T(2.0) * ypoly.b() * B * A + ypoly.c() * A;
        delta = ypoly.a() * B * B * B + ypoly.b() * B * B + ypoly.c() * B + ypoly.d();
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else if (p < T(0.0) && pX.b() > T(0.0))
    {
        mType = EasingCubicBezier<T>::X3COSH;
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else if (p < T(0.0) && pX.b() < T(0.0))
    {
        mType = EasingCubicBezier<T>::X3COSH;
        ecb = { mType, alpha, beta, gamma, delta, l, k };
    }
    else
    {
        assert(false);
    }
    return ecb;
}

static std::pair<double,double> test_easing_creation(std::array<float, 4> P_Y,
                                              std::array<float, 4> P_X,
                                              std::array<double, 4> P_Xd)
{
    float k = std::numeric_limits<float>::epsilon();
    std::array<double, 4> P_Yd = { (double)P_Y[0], (double)P_Y[1], (double)P_Y[2], (double)P_Y[3] };
    std::array<float, 4> P_Xa = { P_X[0], P_X[1] - k, P_X[2], P_X[3] };
    std::array<double, 4> P_Xad = { P_Xd[0], P_Xd[1] - k, P_Xd[2], P_Xd[3] };

    std::array<float, 4> y0 = EasingCubicBezierf::calculatePolynomial(P_Y, 0.0f, 1.0f);
    std::array<float, 4> a0 = EasingCubicBezierf::calculatePolynomial(P_X, 0.0f, 1.0f);
    std::array<float, 4> a0a = EasingCubicBezierf::calculatePolynomial(P_Xa, 0.0f, 1.0f);

    EasingCubicBezierf X;
    EasingCubicBezierf Xa;
    EasingCubicBezierf Xb;
    X.makeEasingFromPolynomial(a0, y0, std::numeric_limits<float>::epsilon());
    Xa.makeEasingFromPolynomial(a0a, y0, std::numeric_limits<float>::epsilon());
    Xb.makeEasingFromPolynomial(a0a, y0, std::numeric_limits<float>::epsilon() * 10.0f);

    std::array<double, 4> y0_d = EasingCubicBezierd::calculatePolynomial(P_Yd, 0.0, 1.0);
    std::array<double, 4> a0a_d = EasingCubicBezierd::calculatePolynomial(P_Xad, 0.0, 1.0);
    EasingCubicBezierd Xa_d;
    Xa_d.makeEasingFromPolynomial(a0a_d, y0_d, (double)std::numeric_limits<float>::epsilon()*10.0);
    float t1 = X.evaluate(0.0f);
    double mmaxa = 0.0;
    double mmaxb = 0.0;
    for (int i = 0; i <= 1024; ++i)
    {
        float t = i / 1024.0f;
        double ta = Xa.evaluate(t);
        double tb = Xb.evaluate(t);
        double d = Xa_d.evaluate(t);
        mmaxa = std::max(std::abs(ta - d), mmaxa);
        mmaxb = std::max(std::abs(tb - d), mmaxb);
    }
    return { mmaxa, mmaxb };
}

void test_easing_creation2()
{
    float k = std::numeric_limits<float>::epsilon()*1.0f;
    std::array<float, 4> P_Y = { 0.0f, 0.5f, 0.5f, 1.0f };
    std::array<double, 4> P_Y_d = { 0.0, 0.5, 0.5, 1.0 };

    std::array<float, 4> P_X0 = { 0.0f, 2.0f / 6.0f, 4.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X1 = { 0.0f, 1.0f / 6.0f, 3.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X2 = { 0.0f, 0.0f / 6.0f, 0.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X3 = { 0.0f, 1.0f / 6.0f, 5.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X4 = { 0.0f, 3.0f / 6.0f, 4.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X5 = { 0.0f, 0.0f / 6.0f, 1.0f / 6.0f, 1.0f };

    std::array<float, 4> P_X0a = { 0.0f, 2.0f / 6.0f - k, 4.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X1a = { 0.0f, 1.0f / 6.0f - k, 3.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X2a = { 0.0f, 0.0f / 6.0f - k, 0.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X3a = { 0.0f, 1.0f / 6.0f - k, 5.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X4a = { 0.0f, 3.0f / 6.0f - k, 4.0f / 6.0f, 1.0f };
    std::array<float, 4> P_X5a = { 0.0f, 0.0f / 6.0f - k, 1.0f / 6.0f, 1.0f };

    std::array<double, 4> P_X0_d = { 0.0, 2.0 / 6.0, 4.0 / 6.0, 1.0 };
    std::array<double, 4> P_X1_d = { 0.0, 1.0 / 6.0, 3.0 / 6.0, 1.0 };
    std::array<double, 4> P_X2_d = { 0.0, 0.0 / 6.0, 0.0 / 6.0, 1.0 };
    std::array<double, 4> P_X3_d = { 0.0, 1.0 / 6.0, 5.0 / 6.0, 1.0 };
    std::array<double, 4> P_X4_d = { 0.0, 3.0 / 6.0, 4.0 / 6.0, 1.0 };
    std::array<double, 4> P_X5_d = { 0.0, 0.0 / 6.0, 1.0 / 6.0, 1.0 };

    //std::array<float, 4> y0 = EasingCubicBezierf::calculatePolynomial(P_Y, 0.0f, 1.0f);

    std::array<float, 4> a0 = EasingCubicBezierf::calculatePolynomial(P_X0, 0.0f, 1.0f);
    std::array<float, 4> a1 = EasingCubicBezierf::calculatePolynomial(P_X1, 0.0f, 1.0f);
    std::array<float, 4> a2 = EasingCubicBezierf::calculatePolynomial(P_X2, 0.0f, 1.0f);
    std::array<float, 4> a3 = EasingCubicBezierf::calculatePolynomial(P_X3, 0.0f, 1.0f);
    std::array<float, 4> a4 = EasingCubicBezierf::calculatePolynomial(P_X4, 0.0f, 1.0f);
    std::array<float, 4> a5 = EasingCubicBezierf::calculatePolynomial(P_X5, 0.0f, 1.0f);

    std::array<float, 4> a0a = EasingCubicBezierf::calculatePolynomial(P_X0a, 0.0f, 1.0f);
    std::array<float, 4> a1a = EasingCubicBezierf::calculatePolynomial(P_X1a, 0.0f, 1.0f);
    std::array<float, 4> a2a = EasingCubicBezierf::calculatePolynomial(P_X2a, 0.0f, 1.0f);
    std::array<float, 4> a3a = EasingCubicBezierf::calculatePolynomial(P_X3a, 0.0f, 1.0f);
    std::array<float, 4> a4a = EasingCubicBezierf::calculatePolynomial(P_X4a, 0.0f, 1.0f);
    std::array<float, 4> a5a = EasingCubicBezierf::calculatePolynomial(P_X5a, 0.0f, 1.0f);

    constexpr float epsilon = std::numeric_limits<float>::epsilon() * 10.0f;

    using typePair = std::pair<EasingCubicBezierf::TYPE, EasingCubicBezierf::TYPE>;

    typePair tp0 = { EasingCubicBezierf::calculateType(a0), EasingCubicBezierf::calculateType(a0a, epsilon) };
    typePair tp1 = { EasingCubicBezierf::calculateType(a1), EasingCubicBezierf::calculateType(a1a, epsilon) };
    typePair tp2 = { EasingCubicBezierf::calculateType(a2), EasingCubicBezierf::calculateType(a2a, epsilon) };
    typePair tp3 = { EasingCubicBezierf::calculateType(a3), EasingCubicBezierf::calculateType(a3a, epsilon) };
    typePair tp4 = { EasingCubicBezierf::calculateType(a4), EasingCubicBezierf::calculateType(a4a, epsilon) };
    typePair tp5 = { EasingCubicBezierf::calculateType(a5), EasingCubicBezierf::calculateType(a5a, epsilon) };
    

    auto maxX0 = test_easing_creation(P_Y, P_X0, P_X0_d);
    auto maxX1 = test_easing_creation(P_Y, P_X1, P_X1_d);
    auto maxX2 = test_easing_creation(P_Y, P_X2, P_X2_d);
    auto maxX3 = test_easing_creation(P_Y, P_X3, P_X3_d);
    auto maxX4 = test_easing_creation(P_Y, P_X4, P_X4_d);
    auto maxX5 = test_easing_creation(P_Y, P_X5, P_X5_d);
}

int main(int argc, char** argv)
{
    test_easing_creation2();
}
