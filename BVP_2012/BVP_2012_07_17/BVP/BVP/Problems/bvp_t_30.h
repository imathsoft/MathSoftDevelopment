#pragma once
#include "ProblemAutnomousGeneralAbstract.h"

template <class T>
class Bvp_t_30 : public ProblemAutonomousGeneralAbstract<T>
{
	T lambda;
public:
	Bvp_t_30(const T lambda_)
	{
		lambda = lambda_;
	}

protected:
	///Nonlineariti
	T Nonlin(const T& u_deriv, const T& u) const override
	{
		return (1-u_deriv)/lambda;
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u) const override { return T(0); }

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u) const override
	{
		return -1/lambda;
	}

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u) const override { return T(0); };

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u) const override { return T(0); }
};