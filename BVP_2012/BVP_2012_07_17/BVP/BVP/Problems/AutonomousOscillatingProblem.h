#ifndef GUARD_AUTONOMOUS_OSCILLATING_PROBLEM
#define GUARD_AUTONOMOUS_OSCILLATING_PROBLEM

#include "ProblemAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"

template <class T>
class AutonomousOscillatingProblem : public ProblemAutonomousAbstract<T>
{
	protected:
	virtual T Nonlin(const T& u) override
	{
		T ln = log(u);
		return  1 - ln*(ln + 1); 
	}

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlin(const T& u) override
	{
		return - (1 + 2*log(u))/u; 
	}
	
	///Second derivative of nonlinearity
	virtual T ddNonlin(const T& u) override
	{
		return (2*log(u) - 1)/auxutils::sqr(u); 
	}

	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) override
	{
		return T(0);
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x)  override
	{
		return T(0);
	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x)  override
	{
		return T(0);
	}
};

#endif