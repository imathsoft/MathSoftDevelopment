#ifndef GUARD_TROESCH_PROBLEM
#define GUARD_TROESCH_PROBLEM

#include "..\Utils\AuxUtils.h"
#include "ProblemAbstract.h"
#include "..\Utils\Exceptions.h"

template <class T>
class TroeschProblem : public ProblemAbstract<T>
{
private:
	T l;
public:
	///Constructor
	TroeschProblem(T lambda)
	{
		l = lambda;
	}

	///Nonlineariti for the Troesch problem
	T inline Nonlin(const T& u)
	{
		T old_sh, sh, t, slu;
		if (abs(u)<0.1)
		{
			int i=1;
			slu=auxutils::sqr(l*u);
			sh=1;
			old_sh=0;
			t=1;
			while (old_sh!=sh) 
			{
				t=t*slu/(2*i*(2*i+1));
			    old_sh=sh;
			    sh=sh+t;
			    i++;
			}
			sh=sh*auxutils::sqr(l);
		} else 
		{
			sh=l*sinh(l*u)/u;
		}
		return sh;
	}

	///A method to return std:function wrapper of local method
	virtual std::function<T(const T&)> GetNonLin() override 
	{
		return [=](const T& u){ return Nonlin(u); };
	}

	///Derivative of nonlinearity
	T inline dNonlin(const T& u)
	{
		T old_sh, sh, t, slu;
		if (abs(u)<0.1)
		{
			int i=1;
   			slu=auxutils::sqr(l*u);
			sh=1; sh=sh/3;
			old_sh=0;
			t=sh;
			while (old_sh!=sh) 
			{
				t=t*slu/(2*i*(2*i+3));
			    old_sh=sh;
			    sh=sh+t;
			    i++;
			}
			sh=u*sh*auxutils::sqr(auxutils::sqr(l));
		} else 
		{
			sh=l*(l*cosh(l*u)/u-sinh(l*u)/auxutils::sqr(u));
		}
		return sh;
	}

	///A method to return std:function wrapper of derivative
	virtual std::function<T(const T&)> GetDerivNonLin() override
	{
		return [=](const T& u){ return dNonlin(u); };
	}

	///Derivative of nonlinearity
	T inline ddNonlin(const T& u)
	{
		T old_sh, sh, t, slu;
		if (abs(u)<0.1)
		{
			int i=1;
   			slu=auxutils::sqr(l*u);
			sh = 1; sh = sh/3;
			old_sh = 0;
			t=sh;
			while (old_sh != sh) 
			{
				t = (2*i + 1)*t*slu/((2*i + 3)*2*i*(2*i - 1));
			    old_sh=sh;
			    sh=sh+t;
			    i++;
			}
			sh = sh*auxutils::sqr(auxutils::sqr(l));
		} else 
		{
			sh=l*(sinh(l*u)*(auxutils::sqr(l) + 2/auxutils::sqr(u)) - 
				2*l*cosh(l*u)/u )/u;
		}
		return sh;
	}

	///A method to return std:function wrapper of the second derivative
	virtual std::function<T(const T&)> GetSecondDerivNonLin() override
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
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffGradient()
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
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffGradient()
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
	virtual std::function<T(const T&, const T&, const T&)> GetACoeffInverse()
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
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetACoeffInverseGradient()
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
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeffInverse()
	{
		return [=](const T& derivative, const T& value, const T& argument )			
		{
			///B_{i} = - N(u_{i})u_{i}(x_{i}^{'})^{2}
		    return - Nonlin(argument)*argument*auxutils::sqr(derivative);
		};
	}

	///A method to return std:function wrapper of the gradient of B coefficient for inverse problem
	virtual std::function<std::array<T, 3>(const T&, const T&, const T&)> GetBCoeffInverseGradient()
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