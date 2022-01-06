#pragma once
#include <array>
#include "../Utils/Utils.h"

using namespace auxutils;

template<class R, int Dim>
class dual
{
private:
	/// <summary>
	/// "Real" part of the dual number
	/// </summary>
	R x{};

	/// <summary>
	/// multi dimensional dual part of the "dual" number
	/// </summary>
	std::array<R, Dim> d{};

	/// <summary>
	/// Scales the dual part by the given factor
	/// </summary>
	void scale_dual_part(const R& scale_factor)
	{
		for (int i = 0; i < Dim; i++)
			d[i] *= scale_factor;
	}

public:
	/// <summary>
	/// Getter the real component
	/// </summary>
	/// <returns></returns>
	R Real() const { return x; }

	/// <summary>
	/// Setter for the real component
	/// </summary>
	/// <returns></returns>
	R& Real() { return x; }

	/// <summary>
	/// Getter for the dual component 
	/// </summary>
	const std::array<R, Dim>& Dual() const
	{
		return d;
	}

	/// <summary>
	/// Setter for the dual component 
	/// </summary>
	/// <returns></returns>
	std::array<R, Dim>& Dual()
	{
		return d;
	}

	dual() = default;

	/// <summary>
	/// Constructor from a scalar type
	/// </summary>
	/// <param name="scalar">The scalar value to initialize the dual number</param>
	dual(const R& scalar) : x{ scalar }, d{}
	{}

	/// <summary>
	/// Constructor for "full" initialization
	/// </summary>
	dual(const R& scalar, const std::array<R, Dim>& dual_array) : x{ scalar }, d{ dual_array }
	{}

	/// <summary>
	/// Composite assignment operator +=
	/// </summary>
	dual<R, Dim>& operator +=(const dual<R, Dim>& rhs)
	{
		x += rhs.x;
		for (int i = 0; i < Dim; i++)
			d[i] += rhs.d[i];

		return *this;
	}

	/// <summary>
	/// Composite assignment operator -=
	/// </summary>
	dual<R, Dim>& operator -=(const dual<R, Dim>& rhs)
	{
		x -= rhs.x;
		for (int i = 0; i < Dim; i++)
			d[i] -= rhs.d[i];

		return *this;
	}

	/// <summary>
	/// Composite assignment operator *=
	/// </summary>
	dual<R, Dim>& operator *=(const dual<R, Dim>& rhs)
	{
		for (int i = 0; i < Dim; i++)
			d[i] = rhs.d[i]*x + d[i]*rhs.x;

		x *= rhs.x;

		return *this;
	}

	/// <summary>
	/// Composite assignment operator /=
	/// </summary>
	dual<R, Dim>& operator /=(const dual<R, Dim>& rhs)
	{
		const auto denom = R(1) / (rhs.x);
		const auto  denom_sqr = denom * denom;
		for (int i = 0; i < Dim; i++)
			d[i] = d[i] * denom - x * rhs.d[i] * denom_sqr;

		x *= denom;

		return *this;
	}

	/// <summary>
	/// Unary minus operator
	/// </summary>
	friend dual<R, Dim> operator -(dual<R, Dim> arg)
	{
		arg.x = -arg.x;
		arg.scale_dual_part(R(-1));

		return arg;
	}

	/// <summary>
	/// Binary "+" operator
	/// </summary>
	friend dual<R, Dim> operator +(dual<R, Dim> lhs, const dual<R, Dim>& rhs)
	{
		return lhs += rhs;
	}

	/// <summary>
	/// Binary "-" operator
	/// </summary>
	friend dual<R, Dim> operator -(dual<R, Dim> lhs, const dual<R, Dim>& rhs)
	{
		return lhs -= rhs;
	}

	/// <summary>
	/// Binary "*" operator
	/// </summary>
	friend dual<R, Dim> operator *(dual<R, Dim> lhs, const dual<R, Dim>& rhs)
	{
		return lhs *= rhs;
	}

	/// <summary>
	/// Binary "/" operator
	/// </summary>
	friend dual<R, Dim> operator /(dual<R, Dim> lhs, const dual<R, Dim>& rhs)
	{
		return lhs /= rhs;
	}

	/// <summary>
	/// Sin function
	/// </summary>
	friend dual<R, Dim> Sin(dual<R, Dim> arg)
	{
		arg.scale_dual_part(Cos(arg.x));
		arg.x = Sin(arg.x);

		return arg;
	}

	/// <summary>
	/// Sin function
	/// </summary>
	friend dual<R, Dim> Cos(dual<R, Dim> arg)
	{
		arg.scale_dual_part(-Sin(arg.x));
		arg.x = Cos(arg.x);

		return arg;
	}

	/// <sumary>
	/// Exponent function
	/// </summary>
	friend dual<R, Dim> Exp(dual<R, Dim> arg)
	{
		arg.x = Exp(arg.x);
		arg.scale_dual_part(arg.x);

		return arg;
	}

	///<summary>
	/// Hyperbolic sine function
	/// </summary>
	friend dual<R, Dim> Sinh(dual<R, Dim> arg)
	{
		arg.scale_dual_part(Cosh(arg.x));
		arg.x = Sinh(arg.x);

		return arg;
	}

	///<summary>
	/// Hyperbolic cosine function
	/// </summary>
	friend dual<R, Dim> Cosh(dual<R, Dim> arg)
	{
		arg.scale_dual_part(Sinh(arg.x));
		arg.x = Cosh(arg.x);

		return arg;
	}

	/// <summary>
	/// Square root function
	/// </summary>
	friend dual<R, Dim> Sqrt(dual<R, Dim> arg)
	{
		arg.x = Sqrt(arg.x);
		arg.scale_dual_part(R(1)/(R(2)*arg.x));

		return arg;
	}
};