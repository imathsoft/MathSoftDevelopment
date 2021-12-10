#ifndef GUARD_TROESCH_PROBLEM
#define GUARD_TROESCH_PROBLEM

#include "..\Utils\AuxUtils.h"
#include "ProblemAutonomousAbstract.h"
#include "..\SpecialFunctions/Sinhc.h"

using namespace SpecianFunctions;

template <class T>
class TroeschProblem : public ProblemAutonomousAbstract<T>
{
private:
	T l;
public:
	///Constructor
	TroeschProblem(T lambda)
	{
		l = lambda;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new TroeschProblem<T>(l));
	}

protected:
	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) const override
	{
		return T(0);
	}

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x) const  override
	{
		return T(0);
	}

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x) const  override
	{
		return T(0);
	}

	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) const override
	{		
		return Sinhc<T>::Func(u, l);
	}

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) const override
	{
		return Sinhc<T>::Deriv(u, l);
	}

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) const override
	{
		return Sinhc<T>::DDeriv(u, l);
	}
};
#endif 