#ifndef GUARD_GRADIENT_VECTOR
#define GUARD_GRADIENT_VECTOR

#include <Eigen/Dense>
#include <Eigen/MPRealSupport>

template <class T, int _Cols>
class GradientVector : public Eigen::Matrix<T, 1,_Cols>
{
private:
    typedef Eigen::Matrix<T, 1,_Cols> super;
public:
	//Constructor
    GradientVector()
	{};
	//Destructor
	~GradientVector()
	{};
    // Assignment operator
	inline const GradientVector<T,_Cols>&  operator=(const Eigen::Matrix<T, 1,_Cols>& gv)
	{
		super::operator=(gv);
		return *this;
	};

	inline static int Length()
	{
		return _Cols;
	}

	inline GradientVector<T, _Cols>  operator*(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result;
		result[0,0] = (*this)[0,0]*gv[0,0];
		for (int i = 1; i< _Cols; ++i)
			result[0,i] = (*this)[0,0]*gv[0,i] + (*this)[0,i]*gv[0,0];
		return result;
	}

	inline GradientVector<T, _Cols>  operator+(const GradientVector<T, _Cols>& gv) const
	{
		GradientVector<T, _Cols> result;
		result = super::operator+(gv); 
		return result;
	}

	inline GradientVector<T, _Cols>  operator/(const double d) const
	{
		GradientVector<T, _Cols> result;
		result = super::operator/(d);
		return result;
	}

	template<class U, int Size>
	friend inline std::ostream& operator << (std::ostream& out, const GradientVector<U, Size>& gv);

	template<class U, int Size>
	friend inline U abs(const GradientVector<U, Size>& gv);
};

	template<class U, int Size>
	inline std::ostream& operator << (std::ostream& out, const GradientVector<U, Size>& gv)
	{
		out << gv.transpose();
		return out;
	}

	template<class U, int Size>
	inline U abs(const GradientVector<U, Size>& gv)
	{
		U a = 0;
		for (int i = 0; i < gv.Length(); ++i)
		{
			a = a + abs(gv[i]);
		}
		return a;
	}

	const int FULL_GRAD_VECT_SIZE = 6;
	const int SMALL_GRAD_VECT_SIZE = 2;

	template<class T>
	class GVTypes
	{
	public : 
		typedef GradientVector<T, SMALL_GRAD_VECT_SIZE> SmallGradientVector;
		typedef GradientVector<T, FULL_GRAD_VECT_SIZE> FullGradientVector;
	};

#endif