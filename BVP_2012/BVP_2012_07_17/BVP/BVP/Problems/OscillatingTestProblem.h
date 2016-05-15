#ifndef GUARD_OSCILLATING_TEST_PROBLEM
#define GUARD_OSCILLATING_TEST_PROBLEM

#include "ProblemAutonomousAbstract.h"

/// u''(x) = cos(u(x))*u(x)

template<class T>
class OscillatingTestProblem : public ProblemAutonomousAbstract<T>
{
	protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) override
	{
		return cos(u);
	}

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) override
	{
		return - sin(u);
	}

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) override
	{
		return - cos(u);
	}
};

#endif
