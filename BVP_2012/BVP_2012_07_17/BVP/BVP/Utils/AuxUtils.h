#ifndef GUARD_AUXUTILS
#define GUARD_AUXUTILS

#include <mpreal.h>
#include <fstream>

#include <boost\multiprecision\cpp_dec_float.hpp>
#include <boost/multiprecision/debug_adaptor.hpp> 

//typedef boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::cpp_dec_float<50>>, boost::multiprecision::et_off> float_50_noet;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30>, boost::multiprecision::et_off> float_50_noet;

namespace auxutils
{
	inline double Sqrt(double d)
	{
		return std::sqrt((double)d);
	}

	inline mpfr::mpreal Sqrt(mpfr::mpreal mpr)
	{
		return mpfr::sqrt(mpr);
	}

	inline float_50_noet Sqrt(float_50_noet mpr)
	{
		return boost::multiprecision::sqrt(mpr);
	}

	inline double sqr(double d)
	{
		return d*d;
	}

	inline float sqr(float f)
	{
		return f*f;
	}

	inline float_50_noet sqr(float_50_noet mpr)
	{
		return mpr*mpr;
	}

	inline mpfr::mpreal sqr(mpfr::mpreal mpr)
	{
		return mpfr::sqr(mpr);
	}

	inline int sqr(int intVal)
	{
		return intVal*intVal;
	}

	inline int sgn(float_50_noet mpr)
	{
		return boost::math::sign(mpr);
	}

	inline int sgn(double x)
	{
		return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
	}

	inline float_50_noet abs(float_50_noet mpr)
	{
		return boost::multiprecision::abs(mpr);
	}

	inline float_50_noet max(float_50_noet mpr1, 
		float_50_noet mpr2)
	{
		return mpr1 > mpr2 ? mpr1 : mpr2;
	}

	inline float_50_noet min(float_50_noet mpr1, 
		float_50_noet mpr2)
	{
		return mpr1 < mpr2 ? mpr1 : mpr2;
	}

	///Square root approximation using Halley's method (that is faster than the Newton's one)
	///Is here for two reasons:
	///1) sqrt() produces ambiguity for doubles and floats 
	///2) In some cases we do not need an exact value of a square root
	template <class T>
	T RoughSqrt(T a)
	{
		T curVal = 1;
		T currValSquared = 1;
		T correction;

		do
		{
			correction = 2*curVal*(currValSquared - a)/(3*currValSquared + a);
			curVal = curVal - correction;
			currValSquared = curVal*curVal;
		} while(abs(correction) > 1);

		return curVal;
	}

    void WriteToStream(std::ofstream& stream, mpfr::mpreal value );

	void WriteToStream(std::ofstream& stream, float_50_noet value );

	void WriteToStream(std::ofstream& stream, float value );

	void WriteToStream(std::ofstream& stream, double value );

	///Method to write sparse matrix into a text file
	template <class T>
	void SaveToFile(std::vector<T> mesh, char* filename, std::streamsize pecision = 15)
	{
		 ofstream file;
		 file.precision(pecision);
		 file.open (filename);
		 for (std::vector<T>::iterator m = mesh.begin(); m != mesh.end(); ++m)
		 {
			 file << (*m);
			 file << endl;
		 }
         file.close();
	}
};

#endif