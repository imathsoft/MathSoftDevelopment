#pragma once
#include <functional>
#include "JordanBlock.h"

template <class T, typename Func>
class DifferentiableFunction
{
private :

	Func _function;

public:
	//Constructor
	DifferentiableFunction(Func function)
	{
		_function = function;
	}

	T Evaluate(const T& argument)
	{
		return _function(argument);
	}

	T Evaluate(const T& argument, T& firstDerivative)
	{
		std::array<T, 2> source = {argument, T(1)};
		JordanBlock<T, 2> result = _function(JordanBlock<T, 2>(source));

		firstDerivative = result[1];

		return result[0];
	}

	T Evaluate(const T& argument, T& firstDerivative, T& secondDerivative)
	{
		std::array<T, 3> source = {argument, T(1), T(0)};
		JordanBlock<T, 3> result = _function(JordanBlock<T, 3>(source));

		firstDerivative = result[1];
		secondDerivative = 2*result[2];

		return result[0];
	}

	T Evaluate(const T& argument, T& firstDerivative, T& secondDerivative, T& thirdDerivative)
	{
		std::array<T, 4> source = {argument, T(1), T(0), T(0)};
		JordanBlock<T, 4> result = _function(JordanBlock<T, 4>(source));

		firstDerivative = result[1];
		secondDerivative = 2*result[2];
		thirdDerivative = 6*result[3];

		return result[0];
	}
};

template<class T, typename F>
DifferentiableFunction<T, F> DifferentiableFunctionFactory(const F& func)
{
	return DifferentiableFunction<T, F>(func);
}
