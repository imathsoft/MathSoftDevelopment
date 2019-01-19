#pragma once
#include "ProblemAbstract.h"
#include "../Utils/AuxUtils.h"

template <class T>
class ProblemAutonomousGeneralAbstract : public ProblemAbstract<T>
{
protected:
	///Nonlineariti
	virtual T Nonlin(const T& u_deriv, const T& u) const = 0;

	///Derivative of nonlinearity with respect to u
	virtual T dNonlin_du(const T& u_deriv, const T& u) const = 0;

	///Derivative of nonlinearity with respect to derivative of u_deriv
	virtual T dNonlin_du_deriv(const T& u_deriv, const T& u) const = 0;

	///Second derivative of nonlinearity with respect to u
	virtual T ddNonlin_ddu(const T& u_deriv, const T& u) const = 0;

	///Second derivative of nonlinearity with respect to u_deriv
	virtual T ddNonlin_ddu_deriv(const T& u_deriv, const T& u) const = 0;

	///Second derivative of nonlinearity with respect to u and u_deriv
	virtual T ddNonlin_dudu_deriv(const T& u_deriv, const T& u) const = 0;

public:
	///A method to return std:function wrapper of the A coefficient
	T GetACoeff(const T& derivative, const T& value, const T& argument) const override
	{
		return dNonlin_du(derivative, value)*derivative + dNonlin_du_deriv(derivative, value)*Nonlin(derivative, value)*value;
	}

	///A method to return std:function wrapper of the gradient of A coefficient
	std::array<T, 3> GetACoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;

		const auto dN_du_deriv = dNonlin_du_deriv(derivative, value);
		const auto dN_du = dNonlin_du(derivative, value);
		const auto N = Nonlin(derivative, value);
		const auto ddN_dudu_deriv = ddNonlin_dudu_deriv(derivative, value);
		const auto ddN_ddu_deriv = ddNonlin_ddu_deriv(derivative, value);

		result[0] = ddN_dudu_deriv*derivative + dN_du + ddN_ddu_deriv*N*value + auxutils::sqr(dN_du_deriv)*value;
		result[1] = ddNonlin_ddu(derivative, value)*derivative + ddN_dudu_deriv*N*value + dN_du_deriv*(dN_du*value + N);
		result[2] = T(0);

		return result;
	}

	///A method to return std:function wrapper of the B coefficient
	T GetBCoeff(const T& derivative, const T& value, const T& argument) const override
	{
		return Nonlin(derivative, value);
	}

	///A method to return std:function wrapper of the gradient of B coefficient
	std::array<T, 3> GetBCoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;

		result[0] = dNonlin_du_deriv(derivative, value);
		result[1] = dNonlin_du(derivative, value);
		result[2] = T(0);

		return result;
	}

	///A method to return std:function wrapper of the A coefficient for inverse problem
	T GetACoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
		const auto derivative_inverse = T(1)/derivative;
		const auto N = Nonlin(derivative_inverse, argument);
		const auto dN_du = dNonlin_du(derivative_inverse, argument);
		const auto dN_du_deriv = dNonlin_du_deriv(derivative_inverse, argument);
		return -dN_du*argument - N - dN_du_deriv*N*auxutils::sqr(argument)*derivative;
	}

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	std::array<T, 3> GetACoeffInverseGradient(const T& derivative, const T& value, const T& argument ) const override
	{
		std::array<T, 3> result;
		const auto derivative_inverse = T(1)/derivative;
		const auto N = Nonlin(derivative_inverse, argument);
		const auto dN_du_deriv = dNonlin_du_deriv(derivative_inverse, argument);
		const auto ddN_dudu_deriv = ddNonlin_dudu_deriv(derivative_inverse, argument);
		const auto ddN_ddu_deriv = ddNonlin_ddu_deriv(derivative_inverse, argument);
		const auto dN_du = dNonlin_du(derivative_inverse, argument);
		const auto ddN_ddu = ddNonlin_ddu(derivative_inverse, argument);
		const auto sqr_argument = auxutils::sqr(argument);

		result[0] = (ddN_dudu_deriv*argument + dN_du_deriv + (ddN_ddu_deriv*N + auxutils::sqr(dN_du_deriv))*sqr_argument*derivative)*auxutils::sqr(derivative_inverse) -
			dN_du_deriv*N*sqr_argument;
		result[1] = T(0);
		result[2] = -ddN_ddu*argument - 2*dN_du - (ddN_dudu_deriv*N + dN_du_deriv*dN_du)*sqr_argument*derivative - 2*dN_du_deriv*N*argument*derivative;
		return result;
	}

	///A method to return std:function wrapper of the B coefficient for inverse problem
	T GetBCoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
		return -Nonlin(1/derivative, argument)*argument;
	}

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	std::array<T, 3> GetBCoeffInverseGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		const auto derivative_inverse = T(1)/derivative;
		result[0] = dNonlin_du_deriv(derivative_inverse, argument)*argument*auxutils::sqr(derivative_inverse);
		result[1] = T(0);
		result[2] = -dNonlin_du(derivative_inverse, argument)*argument - Nonlin(derivative_inverse, argument);
		return result;
	}

	///Returns value of E coefficient at the given point x
	T GetECoeff(const T& x) const override {return T(0); };

	///Returns value of F coefficient at the given point x
	T GetFCoeff(const T& x) const override { return T(0); };

	///Returns derivative of E coefficient with respect to x
	T GetdEdX(const T& x) const override {return T(0); };

	///Returns derivative of F coefficient with respect to x
	T GetdFdX(const T& x) const override { return T(0); };

	///Returns value of the Inverse step function, calculated with the given set of parameters
	InitCondition<T> step_inverse(const T& A, const T& B, const T& C, const T& D, const T& h, const T& precision) const override
	{
		return XI_Bernoulli(A, B, C, D, h, precision);
	}

	///Specific version of the step_inverse function which simultaneously computes gradient of the function with respect to its parameters
	X_Func_Gradient<T> step_inverse_gradient(const T& A, const T& B, const T& C, const T& D, const T& h, const T& precision) const override
	{
		return X_Func_Gradient<T>::XI_Func_Bernoulli_Gradient(A, B, C, D, h, precision);
	}

	T get_optimal_step_inverse(const T& A, const T& B, const T& C, const T& step_desired) const override
	{
		T test = abs(step_desired*C*C*(A*step_desired + 2*B));
		T step_optimal = (test > 0.5) ? 1/abs(2*test) : abs(step_desired);
		return min(step_optimal, abs(step_desired));
	}
};