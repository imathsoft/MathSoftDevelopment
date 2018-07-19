#pragma once

#include "ProblemAbstract.h"

///An abstract class (interface) that represents a non-uniform term of a differential equation 
template <class T>
class ProblemNonUniformAbstract : public ProblemAbstract<T>
{
protected:
	///Right hand side function (non-uniformity term)
	virtual T Phi(const T& x) = 0;

	///The first order derivative of the non-unoformity term
	virtual T dPhi(const T& x) = 0;

	///The second order derivative of the non-unoformity term
	virtual T ddPhi(const T& x) = 0;

public:
	///Returns value of E coefficient at the given point x
	virtual T GetECoeff(const T& x) override
	{
		return dPhi(x);
	}

	///Returns value of F coefficient at the given point x
	virtual T GetFCoeff(const T& x) override
	{
		return Phi(x);
	}
	
	///Returns derivative of E coefficient with respect to x
	virtual T GetdEdX(const T& x) override
	{
		return ddPhi(x);
	}

	///Returns derivative of F coefficient with respect to x
	virtual T GetdFdX(const T& x) override
	{
		return dPhi(x);
	}
};
