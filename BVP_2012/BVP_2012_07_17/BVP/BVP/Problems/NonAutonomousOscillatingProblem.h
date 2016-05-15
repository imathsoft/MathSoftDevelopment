#ifndef GUARD_NONAUTONOMOUS_OSCILLATING_PROBLEM
#define GUARD_NONAUTONOMOUS_OSCILLATING_PROBLEM

#include "ProblemNonAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"

template <class T>
class NonAutonomousOscillatingProblem : public ProblemNonAutonomousAbstract<T>
{
	protected:
	virtual T Nonlin(const T& u, const T& x) override
	{
		return  1 - (sin(x) + 1)*log(u); 
	}

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlinDu(const T& u, const T& x) override
	{
		return - (sin(x) + 1)/u; 
	}
	
	///Derivative of nonlinearity (with respect to x)
	virtual T dNonlinDx(const T& u, const T& x) override
	{
		return - cos(x)*log(u);
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDu(const T& u, const T& x) override
	{
		return (sin(x) + 1)/auxutils::sqr(u); 
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDx(const T& u, const T& x) override
	{
		return -cos(x)/u;
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDxDx(const T& u, const T& x) override
	{
		return sin(x)*log(u);
	}
};

#endif