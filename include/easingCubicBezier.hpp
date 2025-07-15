#include <cmath>
#include <array>
#include <complex>
#include <algorithm>
#include <numbers>

template<typename T>
class EasingCubicBezier
{
	enum class TYPE : uint8_t
	{
		NONE   = 0x0,
		P3     = 0x1,
		X2     = 0x2,
		X3P0   = 0x3,
		X3COS  = 0x4,
		X3SINH = 0x5,
		X3COSH = 0x6
	};

	TYPE mType;
	T a;
	T b;
	T c;
	T d;
	T l;
	T k;

public:
	EasingCubicBezier() = default;

	EasingCubicBezier(const std::array<T, 4>& P_X, const std::array<T, 4>& P_Y) noexcept
	{
		makeEasing(P_X, P_Y);
	}

	T evaluate(T x) const noexcept;

	void makeEasing(const std::array<T, 4>& P_X, const std::array<T, 4>& P_Y) noexcept
	{
		const T x0 = P_X[0];
		const T x3 = P_X[3];
		const std::array<T, 4> m0 = { T(-1), T(3), T(-3), T(1) };
		const std::array<T, 4> m1 = { T(3) * x3, T(-3) * (x0 + T(2) * x3), T(3) * (T(2) * x0 + x3), T(-3) * x0 };
		const std::array<T, 4> m2 = { T(-3) * x3 * x3, T(3) * x3 * (T(2) * x0 + x3), T(-3) * x0 * (x0 + T(2) * x3), T(3) * x0 * x0 };
		const std::array<T, 4> m3 = { x3 * x3 * x3, T(-3) * x0 * x3 * x3, T(3) * x0 * x0 * x3, -x0 * x0 * x0 };
		const T xRange = T(1) / ((x3 - x0) * (x3 - x0) * (x3 - x0));
		const std::array<T, 4> pX =
		{
			(m0[0] * P_X[0] + m0[1] * P_X[1] + m0[2] * P_X[2] + m0[3] * P_X[3]) * xRange,
			(m1[0] * P_X[0] + m1[1] * P_X[1] + m1[2] * P_X[2] + m1[3] * P_X[3]) * xRange,
			(m2[0] * P_X[0] + m2[1] * P_X[1] + m2[2] * P_X[2] + m2[3] * P_X[3]) * xRange,
			(m3[0] * P_X[0] + m3[1] * P_X[1] + m3[2] * P_X[2] + m3[3] * P_X[3]) * xRange
		};
		const std::array<T, 4> pY =
		{
			(m0[0] * P_Y[0] + m0[1] * P_Y[1] + m0[2] * P_Y[2] + m0[3] * P_Y[3]) * xRange,
			(m1[0] * P_Y[0] + m1[1] * P_Y[1] + m1[2] * P_Y[2] + m1[3] * P_Y[3]) * xRange,
			(m2[0] * P_Y[0] + m2[1] * P_Y[1] + m2[2] * P_Y[2] + m2[3] * P_Y[3]) * xRange,
			(m3[0] * P_Y[0] + m3[1] * P_Y[1] + m3[2] * P_Y[2] + m3[3] * P_Y[3]) * xRange
		};
		const T b_a = pX[1] / pX[0];
		const T c_a = pX[2] / pX[0];
		const T d_a = pX[3] / pX[0];
		const T p = c_a - b_a * b_a / T(3);;
		const T q = T(2) / T(27) * b_a * b_a * b_a - b_a * c_a + d_a / T(3);
		if (pX[0] == T(0) && pX[1] == T(0))
			mType = TYPE::P3;
		else if (pX[0] == T(0))
			mType = TYPE::X2;
		else if (p == T(0))
			mType = TYPE::X3P0;
		else if (pX[0] < T(0))
			mType = TYPE::X3COS;
		else if (p > T(0))
			mType = TYPE::X3SINH;
		else if (p < T(0))
			mType = TYPE::X3COSH;
		else
			mType = TYPE::NONE;
		T A = 0, B = 0;
		switch (mType)
		{
		case TYPE::P3:
			A = T(1) / pX[2];
			B = -pX[3] / pX[2];
			k = T(0);
			l = T(1);
			break;
		case TYPE::X2:
			A = std::copysign(T(1), pX[1]);
			B = -pX[2] / (T(2) * pX[1]);
			k = pX[2] * pX[2] / (T(4) * pX[1] * pX[1]) - pX[3] / pX[1];
			l = T(1) / pX[1];
			break;
		case TYPE::X3P0:
			A = T(1);
			B = -pX[1] / (T(3) * pX[0]);
			k = -q;
			l = T(1) / pX[0];
			break;
		case TYPE::X3COS:
			A = T(2) * std::sqrt(-p / T(3));
			B = -pX[1] / (T(3) * pX[0]);
			k = T(3) * q / (A * p);
			l = T(-3) / (A * p * pX[0]);
			break;
		case TYPE::X3SINH:
			A = T(-2) * std::sqrt(p / T(3));
			B = -pX[1] / (T(3) * pX[0]);
			k = T(-3) * q / (A * p);
			l = T(3) / (A * p * pX[0]);
			break;
		case TYPE::X3COSH:
			A = T(-2) * std::sqrt(-p / T(3)) * std::copysign(T(1), -pX[1]);
			B = -pX[1] / (T(3) * pX[0]);
			k = T(3) * q / (A * p);
			l = T(-3) / (A * p * pX[0]);
			break;
		default:
			break;
		}
		a = pY[0] * A * A * A;
		b = T(3) * pY[0] * B * A * A + pY[1] * A * A;
		c = T(3) * pY[0] * B * B * A + T(2) * pY[1] * B * A + pY[2] * A;
		d = pY[0] * B * B * B + pY[1] * B * B + pY[2] * B + pY[3];
	}
};

template<typename T>
T EasingCubicBezier<T>::evaluate(T x) const noexcept
{
	T phi = l * x + k;
	switch (mType)
	{
	case TYPE::X2:
		phi = std::sqrt(phi < T(0) ? T(0) : phi);
		break;
	case TYPE::X3P0:
		phi = std::cbrt(phi);
		break;
	case TYPE::X3COS:
		phi = std::cos(std::acos(std::clamp(phi, T(-1), T(1))) / T(3) - T(2) / T(3) * std::numbers::pi_v<float>);
		break;
	case TYPE::X3COSH:
		phi = std::cosh(std::acosh(std::complex<T>(phi)) / T(3)).real();
		break;
	case TYPE::X3SINH:
		phi = std::sinh(std::asinh(phi) / T(3));
		break;
	default:
		break;
	}
	return ((a * phi + b) * phi + c) * phi + d;
}

using EasingCubicBezierf = EasingCubicBezier<float>;
using EasingCubicBezierd = EasingCubicBezier<double>;