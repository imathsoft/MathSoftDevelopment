#pragma once

#include "ProblemNonAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"
#include <math.h>

template <class T>
class NonAutonomousNonUniformProblem : public ProblemNonAutonomousAbstract<T>
{
private:
	T _alpha;
protected:

	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u, const T& x) const override
	{
		return exp(2*u + cos(_alpha*x)) - _alpha*_alpha;
	}

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlinDu(const T& u, const T& x) const override
	{
		return 2*exp(2*u + cos(_alpha*x));
	}
	
	///Derivative of nonlinearity (with respect to x)
	virtual T dNonlinDx(const T& u, const T& x) const override
	{
		return -_alpha*sin(_alpha*x)*exp(2*u + cos(_alpha*x));
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDu(const T& u, const T& x) const override
	{
		return 4*exp(2*u + cos(_alpha*x));;
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDx(const T& u, const T& x) const override
	{
		return -2*_alpha*sin(_alpha*x)*exp(2*u + cos(_alpha*x));
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDxDx(const T& u, const T& x) const override
	{
		return _alpha*_alpha*(auxutils::sqr(sin(_alpha*x)) - cos(_alpha*x))*exp(2*u + cos(_alpha*x));
	}

	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) const override
	{
		return -exp(2*sin(_alpha*x) + cos(_alpha*x))*sin(_alpha*x);
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x) const override
	{
		T t1 = _alpha * x;
		T t2 = cos(t1);
		T t5 = sin(t1);
		T t10 = exp(2 * t5 + t2);
		return -(2 * _alpha * t2 - _alpha * t5) * t10 * t5 - t10 * _alpha * t2;

	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x) const override
	{
		T t1 = _alpha * _alpha;
		T t2 = _alpha * x;
		T t3 = sin(t2);
		T t6 = cos(t2);
		T t11 = exp(2 * t3 + t6);
		T t14 = _alpha * t6;
		T t17 = -_alpha * t3 + 2 * t14;
		T t18 = t17 * t17;
		return -(-2 * t1 * t3 - t1 * t6) * t11 * t3 - t18 * t11 * t3 - 2 * t17 * t11 * t14 + t11 * t1 * t3;

	}

public:

	NonAutonomousNonUniformProblem(const T& alpha)
	{
		_alpha = alpha;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new NonAutonomousNonUniformProblem<T>(_alpha));
	}
};