#ifndef GUARD_INITIAL_CONDITION
#define GUARD_INITIAL_CONDITION

#include "GradientVector.h"

///Represents a point on a two times differentiatable curve
template <class T>
struct InitCondition
{
public :

	T Value;
	T Derivative;
	T SecDerivative;
	T Argument;

	template<class U>
	friend inline std::ostream& operator << (std::ostream& out, const InitCondition<U>& ic);

	inline const InitCondition<T>&  operator=(const typename GVTypes<T>::SmallGradientVector& gv)
	{
		Value = gv[0];
		Derivative = gv[1];
		Argument = "0";

		return *this;
	}
};

template<class U>
inline std::ostream& operator << (std::ostream& out, const InitCondition<U>& ic)
{
	out << "Val = " << ic.Value << endl;
	out << "Der = " << ic.Derivative << endl;
	out << "SecDer = " << ic.SecDerivative << endl;
	out << "Argument = " << ic.Argument << endl;
	return out;
};

#endif