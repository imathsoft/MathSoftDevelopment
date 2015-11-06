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
	//Destructor
	~GradientVector()
	{};

	///Assignment operator
	inline GradientVector<T, _Cols>& operator=(const GradientVector<T, _Cols>& gv)
	{
		for (int i = 0; i < _Cols; ++i)
			 (*this)[i] = gv[i];

		return *this;
	}

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
	inline GradientVector<T, _Cols>  operator*(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result;
		result[0] = (*this)[0]*gv[0];
		for (int i = 1; i< _Cols; ++i)
			result[i] = (*this)[0]*gv[i] + (*this)[i]*gv[0];
		return result;
	}

	///Addition operator
	inline GradientVector<T, _Cols>  operator+(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result;
		for (int i = 0; i < _Cols; ++i)
			result[i] = (*this)[i] + gv[i];
		return result;
	}

	inline GradientVector<T, _Cols>&  operator+=(const GradientVector<T, _Cols>& gv) 
	{
		for (int i = 0; i < _Cols; ++i)
			(*this)[i] += gv[i];
		return *this;
	}

	inline GradientVector<T, _Cols>  operator/(const double d) const
	{
		GradientVector<T, _Cols> result;
		for (int i = 0; i < _Cols; ++i)
			result[i] = (*this)[i]/d;
		return result;
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
	const int SMALL_GRAD_VECT_SIZE = 2;

	template<class T>
	class GVTypes
	{
	public : 
		typedef GradientVector<T, SMALL_GRAD_VECT_SIZE> SmallGradientVector;
		typedef GradientVector<T, FULL_GRAD_VECT_SIZE> FullGradientVector;
		typedef GradientVector<T, OPTIMAL_GRAD_VECT_SIZE> OptimalGradientVector;
	};

#endif