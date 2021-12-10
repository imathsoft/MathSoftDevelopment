#ifndef GUARD_INITIAL_CONDITION
#define GUARD_INITIAL_CONDITION

#include "GradientVector.h"

///Represents a point on a two times differentiatable curve
template <class T>
struct InitCondition
{
public :

	using float_t = T;

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

	bool operator<(const InitCondition<T> &ic) const 
	{ 
		return Argument < ic.Argument; 
	}

	inline T NormSquared() const
	{
		return Value * Value + Derivative * Derivative +
			SecDerivative * SecDerivative + Argument * Argument;
	}

	template <bool Deriv>
	inline T NormSquaredPartial() const
	{
		if (Deriv)
			return Derivative * Derivative + Argument * Argument;

		return Value * Value + Argument * Argument;
	}

	template<class U>
	friend inline std::ostream& operator << (std::ostream& out, const InitCondition<U>& ic);

	template<class U>
	friend inline std::istream& operator >> (std::istream& in, const InitCondition<U>& ic);

	inline const InitCondition<T>&  operator=(const typename GVTypes<T>::SmallGradientVector& gv)
	{
		Value = gv[0];
		Derivative = gv[1];
		Argument = "0";

		return *this;
	}
};

template<class U>
inline std::istream& operator >> (std::istream& in, InitCondition<U>& ic)
{
	in >> ic.Value >> ic.Derivative >> ic.SecDerivative >> ic.Argument;

	return in;
}

template<class U>
inline std::ostream& operator << (std::ostream& out, const InitCondition<U>& ic)
{
	out << ic.Value << " " << ic.Derivative << " " << ic.SecDerivative << " " << ic.Argument;
	return out;
};

namespace std {

	template<class T> class numeric_limits<InitCondition<T>> {
        public:
			static const int digits10 = numeric_limits<T>::digits10;
			static const bool is_integer = numeric_limits<T>::is_integer;
			static const int  min_exponent = numeric_limits<T>::min_exponent;
			static const int  max_exponent = numeric_limits<T>::max_exponent;
    };
} 

#endif