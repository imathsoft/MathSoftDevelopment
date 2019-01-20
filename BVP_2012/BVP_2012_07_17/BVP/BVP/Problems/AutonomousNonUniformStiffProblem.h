#pragma once

#include "ProblemAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"

template <class T>
class AutonomousNonUniformStiffProblem : public ProblemAutonomousAbstract<T>
{
private:
	T _lambda;
	T _oneOverSqrtLambda;
	T _oneOverLambda;
protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) override
	{
		return (1+u)*_oneOverLambda;
	}

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) override
	{
		return _oneOverLambda;
	}

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) override
	{
		return T(0);
	}

	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) override
	{
		return -exp(-2*x*_oneOverSqrtLambda)*_oneOverLambda;
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x)  override
	{
		return 2*exp(-2*x*_oneOverSqrtLambda)*_oneOverLambda*_oneOverSqrtLambda;
	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x)  override
	{
		return -4*exp(-2*x*_oneOverSqrtLambda)*_oneOverLambda*_oneOverLambda;
	}

public:

	AutonomousNonUniformStiffProblem(const T& lambda)
	{
		_lambda = lambda;
		_oneOverSqrtLambda = T(1)/auxutils::Sqrt(_lambda);
		_oneOverLambda = _oneOverSqrtLambda*_oneOverSqrtLambda;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new AutonomousNonUniformStiffProblem<T>(lambda));
	}

};