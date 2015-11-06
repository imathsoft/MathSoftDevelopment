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

	//Assignment operator
	inline InitCondition<T>& operator=(const InitCondition<T>& val)
	{
		Value = val.Value;
		Derivative = val.Derivative;
		SecDerivative = val.SecDerivative;
		Argument = val.Argument;
		return *this;
	}

	inline InitCondition<T> operator+(const InitCondition<T>& val) const
	{
		InitCondition<T> result;

		result.Value = Value + val.Value;
		result.Derivative = Derivative + val.Derivative;
		result.SecDerivative = SecDerivative + val.SecDerivative;
		result.Argument = Argument + val.Argument;

		return result;
	}

	inline InitCondition<T> operator-(const InitCondition<T>& val) const
	{
		InitCondition<T> result;

		result.Value = Value - val.Value;
		result.Derivative = Derivative - val.Derivative;
		result.SecDerivative = SecDerivative - val.SecDerivative;
		result.Argument = Argument - val.Argument;

		return result;
	}

	inline InitCondition<T> operator*(const T val) const
	{
		InitCondition<T> result;

		result.Value = Value*val;
		result.Derivative = Derivative*val;
		result.SecDerivative = SecDerivative*val;
		result.Argument = Argument*val;

		return result;
	}

	inline InitCondition<T> operator/(const T val) const
	{
		InitCondition<T> result;

		result.Value = Value/val;
		result.Derivative = Derivative/val;
		result.SecDerivative = SecDerivative/val;
		result.Argument = Argument/val;

		return result;
	}

	inline T NormSquared()
	{
		return Value * Value + Derivative * Derivative +
			SecDerivative * SecDerivative + Argument * Argument;
	}

	inline T NormSquaredNaive()
	{
		return Value * Value + Argument * Argument;
	}

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

namespace std {

	template<class T> class numeric_limits<InitCondition<T>> {
        public:
			static const int digits10 = numeric_limits<T>::digits10;
    };
} 

#endif