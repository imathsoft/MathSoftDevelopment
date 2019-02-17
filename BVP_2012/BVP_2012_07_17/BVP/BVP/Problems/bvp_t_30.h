#pragma once
#include "ProblemGeneralAbstract.h"

template <class T>
class Bvp_t_30 : public ProblemGeneralAbstract<T>
{
	T lambda;
public:
	Bvp_t_30(const T lambda_)
	{
		lambda = lambda_;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new Bvp_t_30<T>(lambda));
	}

protected:
	///Nonlineariti
	T Nonlin(const T& u_deriv, const T& u, const T& x) const override
	{
		return (1-u_deriv)/lambda;
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return -1/lambda;
	}

	///Derivative of nonlinearity with respect to x
	T dNonlin_dx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to x
	T ddNonlin_ddx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u and x
	T ddNonlin_dudx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u_deriv and x
	T ddNonlin_dxdu_deriv(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

};