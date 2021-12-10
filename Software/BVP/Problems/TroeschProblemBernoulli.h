#pragma once

#include "..\Utils\AuxUtils.h"
#include "ProblemGeneralAbstract.h"
#include "..\SpecialFunctions/Sinhc.h"

using namespace SpecianFunctions;

template <class T>
class TroeschProblemBernoulli : public ProblemGeneralAbstract<T>
{
private:
	T l;
public:
	///Constructor
	TroeschProblemBernoulli(T lambda)
	{
		l = lambda;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new TroeschProblem<T>(l));
	}

protected:

	///Nonlineariti
	T Nonlin(const T& u_deriv, const T& u, const T& x) const override
	{
		return Sinhc<T>::Func(u, l);
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u, const T& x) const override
	{
		return Sinhc<T>::Deriv(u, l);;
	}

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Derivative of nonlinearity with respect to x
	T dNonlin_dx(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u, const T& x) const override
	{
		return Sinhc<T>::DDeriv(u, l);
	}

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to x
	T ddNonlin_ddx(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to u and x
	T ddNonlin_dudx(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}

	///Second derivative of nonlinearity with respect to u_deriv and x
	T ddNonlin_dxdu_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return T(0);
	}
};