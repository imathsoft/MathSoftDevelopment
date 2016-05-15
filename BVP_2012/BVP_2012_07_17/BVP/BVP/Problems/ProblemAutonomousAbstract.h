#ifndef GUARD_TROESCH_PROBLEM_AUTONOMOUS_ABSTRACT
#define GUARD_TROESCH_PROBLEM_AUTONOMOUS_ABSTRACT

#include <array>
#include <functional>
#include "ProblemAbstract.h"

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAutonomousAbstract : public ProblemAbstract<T>
{
	protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) = 0;

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) = 0;

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) = 0;

	public:
	///A method to return std:function wrapper of local method
	std::function<T(const T&)> GetNonLin()
	{
		return [=](const T& u){ return Nonlin(u); };
	}

	///A method to return std:function wrapper of derivative
	std::function<T(const T&)> GetDerivNonLin()
	{
		return [=](const T& u){ return dNonlin(u); };
	}

	///A method to return std:function wrapper of the second derivative
	std::function<T(const T&)> GetSecondDerivNonLin()
	{
		return [=](const T& u){ return ddNonlin(u); };
	}

	///A method to return std:function wrapper of the A coefficient
	virtual std::function<T(const T&, const T&, const T&)> GetACoeff() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{ return dNonlin(value)*derivative; };
	}

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffGradient() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			std::array<T, 3> result;
			result[0] = dNonlin(value);
			result[1] = ddNonlin(value)*derivative;
			result[2] = 0; 
			return result; 
		};
	}

	///A method to return std:function wrapper of the B coefficient
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeff() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{ return Nonlin(value); };
	}

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffGradient() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			std::array<T, 3> result;
			result[0] = 0;
			result[1] = dNonlin(value);
			result[2] = 0; 
			return result; 
		};
	}

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetACoeffInverse() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			T N = Nonlin(argument);
			T dN = dNonlin(argument);
			T derivSquared = auxutils::sqr(derivative);

			/// A_{i} = - (N^{'}(u_{i})*u_{i} + N(u_{i}) - 2(N(u_{i})u_{i})^{2}(x_{i}^{'})^{2})(x_{i}^{'})^{2}
			return - (dN * argument + N - 2*auxutils::sqr(N*argument)*derivSquared)*derivSquared;
		};
	}

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffInverseGradient() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			std::array<T, 3> result;
			T N  = Nonlin(argument);
			T dN  = dNonlin(argument);
			T ddN = ddNonlin(argument);
			T dM = dN * argument + N;
			T derivSquared = auxutils::sqr(derivative);

			/// \frac{\partial A_{i}}{\partial x_{i}^{'}} = -2(N^{'}(u_{i})u_{i} + N(u_{i}))x_{i}^{'} + 
			/// 8(N(u_{i})u_{i})^{2}(x_{i}^{'})^{3}
			result[0] = - 2 * dM * derivative + 8 * auxutils::sqr(N * argument) * derivSquared * derivative;
			result[1] = 0;
			/// \frac{\partial A_{i}}{\partial u_{i}} = - (N^{''}(u_{i})u_{i}+2N^{'}(u_{i})) + 
			/// 4(N^{'}(u_{i})u_{i}+N(u_{i}))N(u_{i})u_{i}(x_{i}^{'})^{4}
			result[2] = - (ddN * argument + 2 * dN) * derivSquared + 4 * dM * N * argument * auxutils::sqr(derivSquared); 
			return result; 
		};
	}

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeffInverse() override
	{
		return [=](const T& derivative, const T& value, const T& argument )			
		{
			///B_{i} = - N(u_{i})u_{i}(x_{i}^{'})^{2}
		    return - Nonlin(argument)*argument*auxutils::sqr(derivative);
		};
	}

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffInverseGradient() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			std::array<T, 3> result;
			T N  = Nonlin(argument);
			T dN  = dNonlin(argument);

			///\frac{\partial B_{i}}{\partial x_{i}^{'}} = -2N(u_{i})u_{i}x_{i}^{'}
			result[0] = - 2 * N * argument * derivative;
			result[1] = 0;
			///\frac{\partial B_{i}}{\partial u_{i}} = -(N^{'}(u_{i})u_{i} + N(u_{i}))(x_{i}^{'})^{2}
			result[2] = - (dN*argument + N) * auxutils::sqr(derivative); 

			return result;
		};
	}
};
#endif