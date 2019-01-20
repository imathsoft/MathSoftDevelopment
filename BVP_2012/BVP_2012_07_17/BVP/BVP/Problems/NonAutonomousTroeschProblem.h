#ifndef GUARD_NONAUTONOMOUS_TROESCH_PROBLEM
#define GUARD_NONAUTONOMOUS_TROESCH_PROBLEM

#include "ProblemNonAutonomousAbstract.h"
#include "../SpecialFunctions/Sinhc.h"

using namespace SpecianFunctions;

template <class T>
class NonAutonomousTroeschProblem : public ProblemNonAutonomousAbstract<T>
{
	private:
		T l;
	public:
	///Constructor
	NonAutonomousTroeschProblem(T lambda)
	{
		l = lambda;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new NonAutonomousTroeschProblem<T>(l));
	}

	protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u, const T& x) override
	{
		return Sinhc<T>::Func(u, l)*sinh(x*l);
	}

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlinDu(const T& u, const T& x) override
	{
		return Sinhc<T>::Deriv(u, l)*sinh(x*l);
	}
	
	///Derivative of nonlinearity (with respect to x)
	virtual T dNonlinDx(const T& u, const T& x) override
	{
		return l*Sinhc<T>::Func(u, l)*cosh(x*l);
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDu(const T& u, const T& x) override
	{
		return Sinhc<T>::DDeriv(u, l)*sinh(x*l);
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDx(const T& u, const T& x) override
	{
		return l*Sinhc<T>::Deriv(u, l)*cosh(x*l);
	}

	///Second derivative of nonlinearity
	virtual T ddNonlinDxDx(const T& u, const T& x) override
	{
		return l*l*Sinhc<T>::Func(u, l)*sinh(x*l);
	}
};

#endif