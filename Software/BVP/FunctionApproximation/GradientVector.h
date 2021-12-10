#ifndef GUARD_GRADIENT_VECTOR
#define GUARD_GRADIENT_VECTOR

#include <ostream>
#include <array>

template <class T, int _Cols>
class GradientVector
{
private:
	std::array<T, _Cols> _body;
public:
	//Constructor
    GradientVector()
	{};

	GradientVector(const T& initVal)
	{
		_body[0] = initVal;
		for (int i = 1; i < _Cols; ++i)
			 _body[i] = T(0);
	}

	GradientVector(const T val, const int deriv_place_holder_id)
	{
		for (int i = 1; i < _Cols; i++)
			_body[i] = T(0);

		_body[deriv_place_holder_id] = T(1);
		_body[0] = val;
	}

	//Destructor
	~GradientVector() 
	{};

	//Subscript operator
	inline T operator[](const int index) const
	{
		return _body.at(index);
	}

	//Subscript operator
	inline T& operator[](const int index)
	{
		return _body.at(index);
	}

	///Multipication operator
	inline GradientVector<T, _Cols>  operator*=(const GradientVector<T, _Cols>& gv)
	{
		for (int i = 1; i< _Cols; ++i)
			_body[i] = _body[0]*gv._body[i] + _body[i]*gv._body[0];

		_body[0] = _body[0]*gv._body[0];
		return *this;
	}

	///Multipication operator
	inline GradientVector<T, _Cols>  operator*(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result = *this;
		return result *= gv;
	}

	inline GradientVector<T, _Cols>&  operator+=(const GradientVector<T, _Cols>& gv) 
	{
		for (int i = 0; i < _Cols; ++i)
			_body[i] += gv._body[i];
		return *this;
	}

	///Addition operator
	inline GradientVector<T, _Cols>  operator+(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result = *this;
		return result += gv;
	}

	inline GradientVector<T, _Cols>&  operator-=(const GradientVector<T, _Cols>& gv) 
	{
		for (int i = 0; i < _Cols; ++i)
			_body[i] -= gv._body[i];
		return *this;
	}

	///Subtraction operator
	inline GradientVector<T, _Cols>  operator-(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result = *this;
		return result -= gv;
	}

	inline GradientVector<T, _Cols>  operator*=(const T d)
	{
		for (int i = 0; i < _Cols; ++i)
			_body[i] *= d;
		return *this;
	}

	inline GradientVector<T, _Cols>  operator*(const T d) const
	{
		auto result = *this;
		return result *= d;
	}

	inline GradientVector<T, _Cols>  operator/=(const T d)
	{
		T multiplicant = T(1.0)/d;
		return (*this) *= multiplicant;
	}

	inline GradientVector<T, _Cols>  operator/(const T d) const
	{
		GradientVector<T, _Cols> result = *this;
		return result /= d;
	}

	inline GradientVector<T, _Cols>  operator/=(const GradientVector<T, _Cols>& gv)
	{
		GradientVector<T, _Cols> inverse;
		inverse._body[0] = T(1)/gv._body[0];
		T one_over_denominator_squared = inverse._body[0]*inverse._body[0];

		for (int i = 1; i < _Cols; i++)
			inverse._body[i] = -gv._body[i]*one_over_denominator_squared;

		return *this *= inverse;
	}

	inline GradientVector<T, _Cols>  operator/(const GradientVector<T, _Cols>& gv) const
	{
		auto result = *this;
		return result /= gv;
	}

	inline GradientVector<T, _Cols>& operator<<(const T& value)
	{
		for (int i = 0; i < _Cols - 1; ++i)
			(*this)[i] = (*this)[i + 1];

		(*this)[_Cols - 1] = value;
		return *this;
	}

	template<class U, int Size>
	friend inline std::ostream& operator << (std::ostream& out, const GradientVector<U, Size>& gv);

	template<class U, int Size>
	friend inline U abs(const GradientVector<U, Size>& gv);
};

	template<class T, int Size>
	inline GradientVector<T, Size>  operator/(const T numerator, const GradientVector<T, Size>& denominator)
	{
		return GradientVector<T, Size>(numerator) / denominator;
	}

	template<class U, int Size>
	inline std::ostream& operator << (std::ostream& out, const GradientVector<U, Size>& gv)
	{
		for (int i = 0; i < Size; ++i)
			out << (*this)[i];
		return out;
	}

	template<class U, int Size>
	inline U abs(const GradientVector<U, Size>& gv)
	{
		U a = 0;
		for (int i = 0; i < Size; ++i)
		{
			a = a + abs(gv[i]);
		}
		return a;
	}

	const int FULL_GRAD_VECT_SIZE = 6;
	const int OPTIMAL_GRAD_VECT_SIZE = 5;
	const int OPTIMAL_GRAD_VECT_EXTENDED_SIZE = 7;
	const int SMALL_GRAD_VECT_SIZE = 2;

	template<class T>
	class GVTypes
	{
	public : 
		typedef GradientVector<T, SMALL_GRAD_VECT_SIZE> SmallGradientVector;
		typedef GradientVector<T, FULL_GRAD_VECT_SIZE> FullGradientVector;
		typedef GradientVector<T, OPTIMAL_GRAD_VECT_SIZE> OptimalGradientVector;
		typedef GradientVector<T, OPTIMAL_GRAD_VECT_EXTENDED_SIZE> OptimalGradientVectorExtended;
	};

#endif