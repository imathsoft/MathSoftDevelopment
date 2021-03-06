#ifndef GUARD_TROESCH_PROBLEM_NON_AUTONOMOUS_ABSTRACT
#define GUARD_TROESCH_PROBLEM_NON_AUTONOMOUS_ABSTRACT

#include <array>
#include <functional>
#include "ProblemNonUniformAbstract.h"

///Abstract class to represent problems with different nonlinearities
template <class T>
class ProblemNonAutonomousAbstract : public ProblemNonUniformAbstract<T>
{
	protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u, const T& x) const = 0;

	///Derivative of nonlinearity (with respect to u)
	virtual T dNonlinDu(const T& u, const T& x) const = 0;
	
	///Derivative of nonlinearity (with respect to x)
	virtual T dNonlinDx(const T& u, const T& x) const = 0;

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDu(const T& u, const T& x) const = 0;

	///Second derivative of nonlinearity
	virtual T ddNonlinDuDx(const T& u, const T& x) const = 0;

	///Second derivative of nonlinearity
	virtual T ddNonlinDxDx(const T& u, const T& x) const = 0;

	public:
	///A method to return std:function wrapper of the A coefficient
	virtual T GetACoeff(const T& derivative, const T& value, const T& argument) const override
	{
		// N^{'}_{u}(u(x), x)u^{'}(x) + N^{'}_{x}(u(x), x)
		return dNonlinDu(value, argument)*derivative + dNonlinDx(value, argument); 
	}

	///A method to return std:function wrapper of the gradient of A coefficient
	virtual std::array<T, 3> GetACoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		T ddNDuDx = ddNonlinDuDx(value, argument);
		/// N^{'}_{u}(u, x)
		result[0] = dNonlinDu(value, argument);
		/// N^{''}_{u^{2}}(u, x)*u^{'} + N^{''}_{ux}(u,x)
		result[1] = ddNonlinDuDu(value, argument)*derivative + ddNDuDx;
		/// N^{''}(ux)(u, x)u^{'} + N^{''}_{x^{2}}(u, x)
		result[2] = ddNDuDx*derivative + ddNonlinDxDx(value, argument); 
		return result; 
	}

	///A method to return std:function wrapper of the B coefficient
	virtual T GetBCoeff(const T& derivative, const T& value, const T& argument) const override
	{
		///N(u,x)
		return Nonlin(value, argument); 
	}

	///A method to return std:function wrapper of the gradient of B coefficient
	virtual std::array<T, 3> GetBCoeffGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		result[0] = 0;
		///N^{'}_{u}(u, x)
		result[1] = dNonlinDu(value, argument);
		///N^{'}_{x}(u, x)
		result[2] = dNonlinDx(value, argument); 
		return result; 
	}

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual T GetACoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
		T derivSquared = auxutils::sqr(derivative);
		T N = Nonlin(argument, value);
		// A_{i} = - (N^{'}_{u}(u, x)*u + N^{'}_{x}(u, x)*u*x' + N(u, x) + \Phi(x)*x' - 
		// - 2(N(u, x)*u + \Phi(x))^{2}(x')^{2})(x')^{2}
		return - ((dNonlinDu(argument, value) + dNonlinDx(argument, value)*derivative) * argument 
			+ N + dPhi(value)*derivative - 2*auxutils::sqr(N*argument + Phi(value))*derivSquared)*derivSquared;
	}

	///A method to return std:function wrapper of the gradient of A coefficient for inverse problem
	virtual std::array<T, 3> GetACoeffInverseGradient(const T& derivative, const T& value, const T& argument ) const override
	{
		std::array<T, 3> result;
		T N  = Nonlin(argument, value);
		T dNdx  = dNonlinDx(argument, value);
		T dNdu  = dNonlinDu(argument, value);
		T ddNddu = ddNonlinDuDu(argument, value);
		T ddNddx = ddNonlinDxDx(argument, value);
		T ddNdudx = ddNonlinDuDx(argument, value);
		T dM = dNdu * argument + N;
		T rhs = N * argument + Phi(value);
		T dphi = dPhi(value);
		T derivSquared = auxutils::sqr(derivative);

		/// - 3*(N'_{x} + \Phi^{'}(x))*u*(x')^{2} - 2(N^{'}(u_{i})u_{i} + N(u_{i}))x_{i}^{'} + 
		/// 8(N(u_{i})u_{i})^{2}(x_{i}^{'})^{3}
		result[0] = - 3 * (dNdx * argument + dphi) * derivSquared 
					- 2 * dM * derivative + 8 * auxutils::sqr(rhs) * derivSquared * derivative;
		/// - ((N''_{ux}(u,x)+N''_{xx}(u,x)*x')u + N'_{x} + \Phi^{''}(x)*x')*(x')^{2} + 4*(N'_{x}(u,x)*u + \Phi^{'}(x))*(N(u,x)*u + \Phi(x))*(x')^{4}
		result[1] = - ((ddNdudx + ddNddx * derivative)*argument + dNdx + ddPhi(value)*derivative)*derivSquared
			+ 4*(dNdx*argument + dphi)*rhs*auxutils::sqr(derivSquared);
		/// - ((N''_{u}(u, x)u + N''_{ux}(u, x)*x')*u + 2*N'_{u}(u, x) + N'_{x}(u,x)*x')*(x')^{2} + 
		/// 4(N'_{u}(u, x)u + N(u, x))*(N(u, x)*u + \Phi(x))*(x')^{4}
		result[2] = - ((ddNddu + ddNdudx * derivative) * argument + 2 * dNdu + dNdx*derivative) * derivSquared 
					+ 4 * dM * rhs * auxutils::sqr(derivSquared); 
		return result; 
	}

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual T GetBCoeffInverse(const T& derivative, const T& value, const T& argument) const override
	{
		///- (N(u, x) + \Phi(x))*u*(x')^{2}
		return - (Nonlin(argument, value)*argument + Phi(value))*auxutils::sqr(derivative);
	}

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::array<T, 3> GetBCoeffInverseGradient(const T& derivative, const T& value, const T& argument) const override
	{
		std::array<T, 3> result;
		T N  = Nonlin(argument, value);
		T dNdu  = dNonlinDu(argument, value);
		T dNdx  = dNonlinDx(argument, value);
		T derivativeSquared = auxutils::sqr(derivative);

		/// - 2 * (N(u, x) * u + \Phi(x)) * x'
		result[0] = - 2 * (N * argument + Phi(value)) * derivative;
		/// - (N'_{x}(u, x) * u + \Phi^{'}(x)) * (x')^{2}
		result[1] = - (dNdx * argument + dPhi(value)) * derivativeSquared;
		/// - (N'_{u}(u, x) * u + N(u, x))*(x')^{2}
		result[2] = - (dNdu * argument + N) * derivativeSquared; 

		return result;
	}
};
#endif