#pragma once
#include "ProblemGeneralAbstract.h"

template <class T>
class Bvp_t_x : public ProblemGeneralAbstract<T>
{
	T _l;
	T _one_over_l;
	T _one_over_sqrt_l;

	T E(const T& x) const
	{
		return exp(-T(2)*x*(x+T(2))*_one_over_sqrt_l);
	}
public:
	Bvp_t_x(const T& lambda)
	{
		_l = lambda;
		_one_over_l = T(1)/lambda;
		_one_over_sqrt_l = T(1)/auxutils::Sqrt(lambda);
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new Bvp_t_x<T>(_l));
	}

protected:
	///Nonlineariti
	T Nonlin(const T& u_deriv, const T& u, const T& x) const override
	{
		return (4*(x+T(1))*(x+T(1))*_one_over_l + 2*_one_over_sqrt_l)*E(x)*u*u;
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u, const T& x) const override
	{
		return 2*(4*(x+T(1))*(x+T(1))*_one_over_l + 2*_one_over_sqrt_l)*E(x)*u;
	}

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Derivative of nonlinearity with respect to x
	T dNonlin_dx(const T& u_deriv, const T& u, const T& x) const override
	{
		const auto x_plus_1 = x+T(1);
		return -16*x_plus_1*x_plus_1*x_plus_1*_one_over_l*_one_over_sqrt_l*E(x)*u*u;
	}

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u, const T& x) const override
	{
		return 2*(4*(x+T(1))*(x+T(1))*_one_over_l + 2*_one_over_sqrt_l)*E(x);
	}

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to x
	T ddNonlin_ddx(const T& u_deriv, const T& u, const T& x) const override
	{
		const auto x_plus_1 = x+T(1);
		return 16*x_plus_1*x_plus_1*(4*x_plus_1*x_plus_1*_one_over_sqrt_l - T(3))*_one_over_l*_one_over_sqrt_l*E(x)*u*u;
	}

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to u and x
	T ddNonlin_dudx(const T& u_deriv, const T& u, const T& x) const override
	{
		const auto x_plus_1 = x+T(1);
		return -32*x_plus_1*x_plus_1*x_plus_1*_one_over_l*_one_over_sqrt_l*E(x)*u;
	}

	///Second derivative of nonlinearity with respect to u_deriv and x
	T ddNonlin_dxdu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}
};