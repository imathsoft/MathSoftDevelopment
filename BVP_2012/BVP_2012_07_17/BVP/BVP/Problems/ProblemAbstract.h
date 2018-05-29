#ifndef GUARD_TROESCH_PROBLEM_ABSTRACT
#define GUARD_TROESCH_PROBLEM_ABSTRACT

#include <array>
#include <functional>

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAbstract
{
	public:
	///A method to return std:function wrapper of the A coefficient
	virtual T GetACoeff(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::array<T, 3> GetACoeffGradient(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the B coefficient
	virtual T GetBCoeff(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::array<T, 3> GetBCoeffGradient(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual T GetACoeffInverse(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	virtual std::array<T, 3> GetACoeffInverseGradient(const T& derivative, const T& value, const T& argument ) = 0;

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual T GetBCoeffInverse(const T& derivative, const T& value, const T& argument) = 0;

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::array<T, 3> GetBCoeffInverseGradient(const T& derivative, const T& value, const T& argument) = 0;
};
#endif 