#ifndef GUARD_TROESCH_PROBLEM_ABSTRACT
#define GUARD_TROESCH_PROBLEM_ABSTRACT

#include "..\FunctionApproximation\X_Function.h"

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAbstract
{
	public:
	///A method to return std:function wrapper of Nonlinearity
	virtual std::function<T(const T&)> GetNonLin() = 0;

	///A method to return std:function wrapper of derivative
	virtual std::function<T(const T&)> GetDerivNonLin() = 0;

	///A method to return std:function wrapper of step function
	virtual std::function<T(const int, const int)> GetStepFunc() = 0;

	///A method to return std:function wrapper of step function
	virtual std::function<bool(InitCondition<T>&)> GetCheckFunc() = 0;
};
#endif 