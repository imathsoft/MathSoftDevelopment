#ifndef GUARD_SINHC
#define GUARD_SINHC

namespace SpecianFunctions
{
	template<class T>
	class Sinhc
	{
	public:
		///Hyperbolic sinus cardinal
		static T Func(const T& u, const T& l)
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

		///Derivative of the Sinhc(.) function
		static T Deriv(const T& u, const T& l)
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

		///The second derivative
		static T DDeriv(const T& u, const T& l)
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
}

#endif