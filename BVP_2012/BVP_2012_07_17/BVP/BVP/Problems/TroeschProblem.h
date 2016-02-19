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

protected:
	///Nonlineariti for the Troesch problem
	virtual T Nonlin(const T& u) override
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

	///Derivative of nonlinearity
	virtual T dNonlin(const T& u) override
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

	///Derivative of nonlinearity
	virtual T ddNonlin(const T& u) override
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
};
#endif 