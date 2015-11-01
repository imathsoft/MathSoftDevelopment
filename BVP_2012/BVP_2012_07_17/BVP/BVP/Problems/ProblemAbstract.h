#ifndef GUARD_TROESCH_PROBLEM_ABSTRACT
#define GUARD_TROESCH_PROBLEM_ABSTRACT

#include <array>
#include <functional>

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAbstract
{
	public:
	///A method to return std:function wrapper of Nonlinearity
	virtual std::function<T(const T&)> GetNonLin() = 0;

	///A method to return std:function wrapper of derivative
	virtual std::function<T(const T&)> GetDerivNonLin() = 0;

	///A method to return std:function wrapper of the second derivative
	virtual std::function<T(const T&)> GetSecondDerivNonLin() = 0;

	///A method to return std:function wrapper of the A coefficient
	virtual std::function<T(const T&, const T&, const T&)> GetACoeff() = 0;

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffGradient() = 0;

	///A method to return std:function wrapper of the B coefficient
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeff() = 0;

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffGradient() = 0;

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetACoeffInverse() = 0;

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffInverseGradient() = 0;

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeffInverse() = 0;

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffInverseGradient() = 0;
};
#endif 