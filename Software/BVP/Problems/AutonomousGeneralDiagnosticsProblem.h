#pragma once

#include "ProblemGeneralAbstract.h"

template <class T>
class AutonomousGeneralDiagnosticsProblem : public ProblemGeneralAbstract<T>
{
	T lambda;
public:
	AutonomousGeneralDiagnosticsProblem(const T lambda_)
	{
		lambda = lambda_;
	}

	///Returns pointer to a deep copy of the current instance of the class
	std::unique_ptr<ProblemAbstract<T>> copy() const override
	{
		return std::unique_ptr<ProblemAbstract<T>>(new AutonomousGeneralDiagnosticsProblem<T>(lambda));
	}

protected:
	///Nonlineariti
	T Nonlin(const T& u_deriv, const T& u, const T& x) const override
	{
		return (1-sin(u + 2*u_deriv))/lambda;
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u, const T& x) const override { 
		return -cos(u + 2*u_deriv)/lambda; 
	}

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u, const T& x) const override
	{
		return -2*cos(u + 2*u_deriv)/lambda; 
	}

	///Derivative of nonlinearity with respect to x
	T dNonlin_dx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u, const T& x) const override { 
		return sin(u + 2*u_deriv)/lambda;  
	}

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u, const T& x) const override { 
		return 4*sin(u + 2*u_deriv)/lambda;  
	};

	///Second derivative of nonlinearity with respect to x
	T ddNonlin_ddx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u, const T& x) const override { 
		return 2*sin(u + 2*u_deriv)/lambda;  
	}

	///Second derivative of nonlinearity with respect to u and x
	T ddNonlin_dudx(const T& u_deriv, const T& u, const T& x) const override { return T(0); }

	///Second derivative of nonlinearity with respect to u_deriv and x
	T ddNonlin_dxdu_deriv(const T& u_deriv, const T& u, const T& x) const override { return T(0); }
};