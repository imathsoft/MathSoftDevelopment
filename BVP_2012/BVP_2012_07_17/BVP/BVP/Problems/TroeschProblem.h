#ifndef GUARD_TROESCH_PROBLEM
#define GUARD_TROESCH_PROBLEM

#include "..\Utils\AuxUtils.h"
#include "ProblemAbstract.h"

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

	///A method to return std:function wrapper of the B coefficient
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeff() override
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{ return Nonlin(value); };
	}

	///A method to return std:function wrapper of the A coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetACoeffInverse()
	{
		return [=](const T& derivative, const T& value, const T& argument ) 
		{
			T F = Nonlin(argument);
			T derivSquared = auxutils::sqr(derivative);
			return - (dNonlin(argument)*argument + F - 2*auxutils::sqr(F*argument)*derivSquared)*derivSquared;
		};
	}

	///A method to return std:function wrapper of the B coefficient for inverse problem
	virtual std::function<T(const T&, const T&, const T&)> GetBCoeffInverse()
	{
		return [=](const T& derivative, const T& value, const T& argument )			
		{
		    return - Nonlin(argument)*argument*auxutils::sqr(derivative);
		};
	}
};
#endif 