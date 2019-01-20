#pragma once

#include "ProblemAutnomousGeneralAbstract.h"

template <class T>
class AutonomousGeneralDiagnosticsProblem : public ProblemAutonomousGeneralAbstract<T>
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
	T Nonlin(const T& u_deriv, const T& u) const override
	{
		return (1-sin(u + 2*u_deriv))/lambda;
	}

	///Derivative of nonlinearity with respect to u
	T dNonlin_du(const T& u_deriv, const T& u) const override { 
		return -cos(u + 2*u_deriv)/lambda; 
	}

	///Derivative of nonlinearity with respect to derivative of u_deriv
	T dNonlin_du_deriv(const T& u_deriv, const T& u) const override
	{
		return -2*cos(u + 2*u_deriv)/lambda; 
	}

	///Second derivative of nonlinearity with respect to u
	T ddNonlin_ddu(const T& u_deriv, const T& u) const override { 
		return sin(u + 2*u_deriv)/lambda;  
	}

	///Second derivative of nonlinearity with respect to u_deriv
	T ddNonlin_ddu_deriv(const T& u_deriv, const T& u) const override { 
		return 4*sin(u + 2*u_deriv)/lambda;  
	};

	///Second derivative of nonlinearity with respect to u and u_deriv
	T ddNonlin_dudu_deriv(const T& u_deriv, const T& u) const override { 
		return 2*sin(u + 2*u_deriv)/lambda;  
	}
};