#pragma once

#include "ProblemAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"
#include <math.h>

template <class T>
class AutonomousNonUniformProblem : public ProblemAutonomousAbstract<T>
{
private:
	T _alpha;
protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) const override
	{
		return exp(2*u) - _alpha*_alpha;
	}

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) const override
	{
		return 2*exp(2*u);
	}

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) const override
	{
		return 4*exp(2*u);
	}

	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) const override
	{
		return -exp(2*sin(_alpha*x))*sin(_alpha*x);
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x) const override
	{
		return -_alpha*(cos(_alpha*x) + sin(2*_alpha*x))*exp(2*sin(_alpha*x));
	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x) const override
	{
		return -_alpha*_alpha*(2*cos(2*_alpha*x) - sin(_alpha*x) + 2*cos(_alpha*x)*(cos(_alpha*x) + sin(2*_alpha*x)))*exp(2*sin(_alpha*x));
	}

public:

	AutonomousNonUniformProblem(const T& alpha)
	{
		_alpha = alpha;
	}
};