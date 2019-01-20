#ifndef GUARD_NONAUTONOMOUS_OSCILLATING_PROBLEM
#define GUARD_NONAUTONOMOUS_OSCILLATING_PROBLEM

#include "ProblemNonAutonomousAbstract.h"
#include "../Utils/AuxUtils.h"

template <class T>
class NonAutonomousOscillatingProblem : public ProblemNonAutonomousAbstract<T>
{
	protected:
	virtual T Nonlin(const T& u, const T& x) const override
	{
		return  1 - (sin(x) + 1)*log(u); 
	}

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlinDu(const T& u, const T& x) const override
	{
		return - (sin(x) + 1)/u; 
	}
	
	///Derivative of nonlinearity (with respect to x)
	virtual T dNonlinDx(const T& u, const T& x) const override
	{
		return - cos(x)*log(u);
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDu(const T& u, const T& x) const override
	{
		return (sin(x) + 1)/auxutils::sqr(u); 
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDx(const T& u, const T& x) const override
	{
		return -cos(x)/u;
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDxDx(const T& u, const T& x) const override
	{
		return sin(x)*log(u);
	}

	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) const override
	{
		return T(0);
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x) const override
	{
		return T(0);
	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x) const override
	{
		return T(0);
	}
public:
	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new NonAutonomousOscillatingProblem<T>());
	}
};

#endif