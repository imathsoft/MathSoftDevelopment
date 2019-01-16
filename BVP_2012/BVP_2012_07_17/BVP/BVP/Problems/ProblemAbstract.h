#ifndef GUARD_TROESCH_PROBLEM_ABSTRACT
#define GUARD_TROESCH_PROBLEM_ABSTRACT

#include <array>
#include "../FunctionApproximation/InitialCondition.h"
#include "../FunctionApproximation/X_Function.h"

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAbstract
{
	public:
	///A method to return std:function wrapper of the A coefficient
	virtual T GetACoeff(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::array<T, 3> GetACoeffGradient(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the B coefficient
	virtual T GetBCoeff(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::array<T, 3> GetBCoeffGradient(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual T GetACoeffInverse(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	virtual std::array<T, 3> GetACoeffInverseGradient(const T& derivative, const T& value, const T& argument ) const = 0;

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual T GetBCoeffInverse(const T& derivative, const T& value, const T& argument) const = 0;

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::array<T, 3> GetBCoeffInverseGradient(const T& derivative, const T& value, const T& argument) const = 0;

	///Returns value of E coefficient at the given point x
	virtual T GetECoeff(const T& x) const = 0;

	///Returns value of F coefficient at the given point x
	virtual T GetFCoeff(const T& x) const = 0;

	///Returns derivative of E coefficient with respect to x
	virtual T GetdEdX(const T& x) const = 0;

	///Returns derivative of F coefficient with respect to x
	virtual T GetdFdX(const T& x) const = 0;

	///Returns value of the Streight step functions, calculated with the given set of parameters
	virtual InitCondition<T> step_streight(const T& A, const T& B, const T& C, const T& D, const T& E, const T& F, const T& h, const T& precision) const
	{
		return X4_Func(A, B, C, D, E, F, h, precision);
	}

	///Specific version of the "treight" step function which simultaneously computes gradient of the function with respect to its parameters
	virtual X_Func_Gradient<T> step_streight_gradient(const T& A, const T& B, const T& C, const T& D, const T& E, const T& F, const T& h, const T& precision) const
	{
		return X_Func_Gradient<T>::X3_Func_Gradient(A, B, C, D, E, F, h, precision);
	}

	///Returns "optimal" step size for the given set of parameters which define a piece-wise linear approximation of the "streight" problem
	virtual T get_optimal_step_straight(const T& A, const T& B, const T& step_desired) const
	{
		const T test = max(abs(A*step_desired), abs(B));
		T step_optimal = (test*step_desired*step_desired > 0.5) ? 1/auxutils::RoughSqrt(2*test) : abs(step_desired);
		return min(step_optimal, abs(step_desired));
	}

	///Returns value of the "inverse" step function, calculated with the given set of parameters
	virtual InitCondition<T> step_inverse(const T& A, const T& B, const T& C, const T& D, const T& h, const T& precision) const
	{
		return XI_Func(A, B, C, D, h, precision);
	}

	///Specific version of the step_inverse function which simultaneously computes gradient of the function with respect to its parameters
	virtual X_Func_Gradient<T> step_inverse_gradient(const T& A, const T& B, const T& C, const T& D, const T& h, const T& precision) const
	{
		return X_Func_Gradient<T>::XI_Func_Gradient(A, B, C, D, h, precision);
	}

	///Returns "optimal" step size for the given set of parameters which define a piece-wise linear approximation of the "inverse" problem
	virtual T get_optimal_step_inverse(const T& A, const T& B, const T& C, const T& step_desired) const
	{
		T test = max(abs(A*step_desired), abs(B));
		T step_optimal = (test*step_desired > 0.5) ? 1/abs(2*test) : abs(step_desired);
		return min(step_optimal, abs(step_desired));
	}
};
#endif 