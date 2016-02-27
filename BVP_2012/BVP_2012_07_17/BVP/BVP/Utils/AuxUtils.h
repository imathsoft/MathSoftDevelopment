#ifndef GUARD_AUXUTILS
#define GUARD_AUXUTILS

#include <mpreal.h>
#include <fstream>
#include <strstream>
#include "..\FunctionApproximation\InitialCondition.h"

#include <boost\multiprecision\cpp_dec_float.hpp>
#include <boost/multiprecision/debug_adaptor.hpp> 

//typedef boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::cpp_dec_float<20>>, boost::multiprecision::et_off> float_50_noet;
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

	template <class T>
	inline std::string ToString( const T val)
	{
		std::stringstream ss;
		ss.precision(std::numeric_limits<T>::digits10);
		ss << val;
		return ss.str();
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

	///Method to write vector of knots into a text file
	template <class T>
	void SaveToFile(const std::vector<T> mesh, const char* filename)
	{
		 std::ofstream file;
		 file.precision(std::numeric_limits<T>::digits10);
		 file.open (filename);
		 for (std::vector<T>::const_iterator m = mesh.begin(); m != mesh.end(); ++m)
		 {
			 file << (*m);
			 file << endl;
		 }
         file.close();
	}

	///Method to write vector of knots into a "Maple"-compatible text file
	template <class T>
	void SaveToMapleFile(const std::vector<InitCondition<T>> mesh, const char* filename)
	{
		 std::ofstream file;
		 file.precision(std::numeric_limits<T>::digits10);
		 file.open (filename);
		 for (std::vector<InitCondition<T>>::const_iterator m = mesh.begin(); m != mesh.end(); ++m)
		 {
			 file << (*m).Argument << " " << (*m).Value;
			 file << endl;
		 }
         file.close();
	}

};

#endif