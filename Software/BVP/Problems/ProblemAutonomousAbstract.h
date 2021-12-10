#ifndef GUARD_TROESCH_PROBLEM_AUTONOMOUS_ABSTRACT
#define GUARD_TROESCH_PROBLEM_AUTONOMOUS_ABSTRACT

#include <array>
#include <functional>
#include "ProblemNonUniformAbstract.h"

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemAutonomousAbstract : public ProblemNonUniformAbstract<T>
{
	protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) const = 0;

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) const = 0;

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) const = 0;

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
	virtual T GetACoeff(const T& derivative, const T& value, const T& argument) const override
	{
		return dNonlin(value)*derivative; 
	}

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::array<T, 3> GetACoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		result[0] = dNonlin(value);
		result[1] = ddNonlin(value)*derivative;
		result[2] = 0; 
		return result; 
	}

	///A method to return std:function wrapper of the B coefficient
	virtual T GetBCoeff(const T& derivative, const T& value, const T& argument) const override
	{
		return Nonlin(value);
	}

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::array<T, 3> GetBCoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		result[0] = 0;
		result[1] = dNonlin(value);
		result[2] = 0; 
		return result; 
	}

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual T GetACoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
		T N = Nonlin(argument);
		T dN = dNonlin(argument);
		T derivSquared = auxutils::sqr(derivative);

		/// A_{i} = - (N^{'}(u_{i})*u_{i} + N(u_{i}) + \Phi^{'}(x_{i})*x^{'}_{i} - 2(N(u_{i})u_{i} + \Phi(x_{i}))^{2}(x_{i}^{'})^{2})(x_{i}^{'})^{2}
		return - (dN * argument + N + dPhi(value)*derivative - 2*auxutils::sqr(N*argument + Phi(value))*derivSquared)*derivSquared;
	}

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	virtual std::array<T, 3> GetACoeffInverseGradient(const T& derivative, const T& value, const T& argument ) const override
	{
		std::array<T, 3> result;
		T N  = Nonlin(argument);
		T dN  = dNonlin(argument);
		T ddN = ddNonlin(argument);
		T phi = Phi(value);
		T dphi = dPhi(value);
		T dM = dN * argument + N;
		T rhs = N * argument + phi;
		T derivSquared = auxutils::sqr(derivative);
		T derivQuad = auxutils::sqr(derivSquared);

		/// \frac{\partial A_{i}}{\partial x_{i}^{'}} = -2(N^{'}(u_{i})u_{i} + N(u_{i}))x_{i}^{'} - 3*\Phi^{'}(x_{i})*(x_{i}^{'})^{2} +
		/// 8(N(u_{i})u_{i})^{2}(x_{i}^{'})^{3}
		result[0] = - 2 * dM * derivative - 3 * dphi * derivSquared + 8 * auxutils::sqr(rhs) * derivSquared * derivative;
		///\frac{\partial A_{i}}{\partial x_{i}} = - \Phi^{''}(x_{i})*(x_{i}^{'})^{3} + 4 * \Phi^{'}(x_{i})*(N(u_{i})*u_{i} + \Phi(x_{i}))*(x_{i}^{'})^{4}
		result[1] = - ddPhi(value) * derivative * derivSquared + 4 * dphi * rhs * derivQuad;
		/// \frac{\partial A_{i}}{\partial u_{i}} = - (N^{''}(u_{i})u_{i}+2N^{'}(u_{i})) + 
		/// 4(N^{'}(u_{i})u_{i}+N(u_{i}))(N(u_{i})u_{i} + \Phi(x_{i}))(x_{i}^{'})^{4}
		result[2] = - (ddN * argument + 2 * dN) * derivSquared + 4 * dM * rhs * derivQuad; 
		return result; 
	}

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual T GetBCoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
			///B_{i} = - (N(u_{i})u_{i} + \Phi(u_{i}))(x_{i}^{'})^{2}
		    return - (Nonlin(argument)*argument + Phi(value))*auxutils::sqr(derivative);
	}

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::array<T, 3> GetBCoeffInverseGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		T N  = Nonlin(argument);
		T dN  = dNonlin(argument);
		T derivSquared = auxutils::sqr(derivative);

		///\frac{\partial B_{i}}{\partial x_{i}^{'}} = -2(N(u_{i})u_{i} + \Phi(x))x_{i}^{'}
		result[0] = - 2 * (N * argument + Phi(value)) * derivative;
		///\frac{\partial B_{i}}{\partial x_{i}} = -\Phi^{'}(x)(x_{i}^{'})^{2}
		result[1] = - dPhi(value) * derivSquared;
		///\frac{\partial B_{i}}{\partial u_{i}} = -(N^{'}(u_{i})u_{i} + N(u_{i}))(x_{i}^{'})^{2}
		result[2] = - (dN*argument + N) * derivSquared; 

		return result;
	}
};
#endif