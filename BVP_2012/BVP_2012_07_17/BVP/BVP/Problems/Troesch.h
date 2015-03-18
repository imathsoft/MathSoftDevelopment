#ifndef GUARD_TROESCH
#define GUARD_TROESCH

#include "../AuxUtils.h"
#include "../X_function.h"

/// Class that represents Troesch problem (template version)
/// TO-DO: get rid of this class
template <class T, int l>
class Troesch
{
public:
	static T inline Nonlin(const T& u){
	
		T old_sh, sh, t, slu;
		if (abs(l*u)<0.1){
			int i=1;
			slu=auxutils::sqr(l*u);
			sh=1;
			old_sh=0;
			t=1;
			while (old_sh!=sh) 
			{t=t*slu/(2*i*(2*i+1));
			old_sh=sh;
			sh=sh+t;
			i++;
			}
			sh=sh*auxutils::sqr(l);
		} else {
			sh=l*sinh(l*u)/u;
		}
		return sh;
	
		//return sin(u);
	}
	static T inline dNonlin(const T& u){
	
		T old_sh, sh, t, slu;
		if (abs(l*u)<0.1){
			int i=1;
   			slu=auxutils::sqr(l*u);
			sh=1; sh=sh/3;
			old_sh=0;
			t=sh;
			while (old_sh!=sh) 
			{t=t*slu/(2*i*(2*i+3));
			old_sh=sh;
			sh=sh+t;
			i++;
			}
			int lsqr = l*l;
			sh=u*sh*auxutils::sqr(auxutils::sqr(l));
		} else {
			sh=l*(l*cosh(l*u)/u-sinh(l*u)/auxutils::sqr(u));
		}
		return sh;
	
		//return cos(u);
	}

	static T inline StepFunc(const int stepIndex, const int maStepNumber){// This is a function which can define a non-uniform mesh for the FD-method
		return 1;
	}
	
	static bool CheckFunc(InitCondition<T>& ic)
	{
		T predictedValue = ic.Value + ( "1" - ic.Argument )*ic.Derivative;
		if (abs(predictedValue) > 2)
		{
			ic.Value = predictedValue;
			ic.Argument = "1";
			return false;
		}
		else
			return true;
	}
};
#endif 