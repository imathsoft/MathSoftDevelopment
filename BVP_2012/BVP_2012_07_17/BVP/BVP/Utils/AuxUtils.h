#ifndef GUARD_AUXUTILS
#define GUARD_AUXUTILS

#include <fstream>
#include <strstream>
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\GradientVector.h"

#include <boost\multiprecision\cpp_dec_float.hpp>
#include <boost/multiprecision/debug_adaptor.hpp> 

using namespace boost::multiprecision;

//typedef boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::cpp_dec_float<30>>, boost::multiprecision::et_off> float_50_noet;
typedef number<cpp_dec_float<30>, et_off> float_50_noet;

namespace auxutils
{
	inline double Sqrt(const double d)
	{
		return std::sqrt(d);
	}

	inline float Sqrt(const float d)
	{
		return std::sqrtf(d);
	}

	template <unsigned Digits10>
	inline number<cpp_dec_float<Digits10>, et_off> Sqrt(const number<cpp_dec_float<Digits10>, et_off>& val)
	{
		return boost::multiprecision::sqrt(val);
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Sqrt(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Sqrt(gv[0]);
		auto factor = U(1)/(2*result[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i]*factor;

		return result;
	}

	inline double Exp(const double d)
	{
		return std::exp(d);
	}

	inline float Exp(const float d)
	{
		return std::expf(d);
	}

	template <unsigned Digits10>
	inline number<cpp_dec_float<Digits10>, et_off> Exp(const number<cpp_dec_float<Digits10>, et_off>& val)
	{
		return boost::multiprecision::exp(val);
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Exp(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Exp(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i]*result[0];

		return result;
	}

	inline double Sin(const double d)
	{
		return std::sin(d);
	}

	inline float Sin(const float d)
	{
		return std::sinf(d);
	}

	template <unsigned Digits10>
	inline number<cpp_dec_float<Digits10>, et_off> Sin(const number<cpp_dec_float<Digits10>, et_off>& val)
	{
		return boost::multiprecision::sin(val);
	}

	inline double Cos(const double d)
	{
		return std::cos(d);
	}

	inline float Cos(const float d)
	{
		return std::cosf(d);
	}

	template <unsigned Digits10>
	inline number<cpp_dec_float<Digits10>, et_off> Cos(const number<cpp_dec_float<Digits10>, et_off>& val)
	{
		return boost::multiprecision::cos(val);
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Sin(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Sin(gv[0]);
		auto deriv = Cos(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i]*deriv;

		return result;
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Cos(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Cos(gv[0]);
		auto deriv = -Sin(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i]*deriv;

		return result;
	}

	template <class T>
	inline T sqr(const T& d)
	{
		return d*d;
	}

	template <class T>
	inline int sgn(T x)
	{
		return (x > T(0)) ? 1 : ((x < T(0)) ? -1 : 0);
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

	template <unsigned Digits10>
	void WriteToStream(std::ofstream& stream, number<cpp_dec_float<Digits10>, et_off>& value )
	{
		stream << value << " ";
	}

	void WriteToStream(std::ofstream& stream, float value );

	void WriteToStream(std::ofstream& stream, double value );

	///Method to write vector of knots into a text file
	template <class T>
	void SaveToFile(const std::vector<T>& mesh, const char* filename)
	{
		 std::ofstream file;
		 file.open (filename);
		 file.precision(std::numeric_limits<T>::digits10 + 2);
		 for (std::vector<T>::const_iterator m = mesh.begin(); m != mesh.end(); ++m)
		 {
			 file << std::endl << (*m);
		 }
         file.close();
	}

	///Method to write vector of knots into a text file
	template <class T>
	std::vector<T> ReadFromFile(const char* filename)
	{
		std::vector<T> mesh;
		std::ifstream file;
		file.open (filename);
		while (!file.eof())
		{
			T data;
			file >> data;
			mesh.push_back(data);
		}
        file.close();

		return mesh;
	}


	///Method to write vector of knots into a "Maple"-compatible text file
	template <class T>
	void SaveToMapleFile(const std::vector<InitCondition<T>>& mesh, const char* filename, 
		bool saveDerivatives = false)
	{
		 std::ofstream file;
		 file.precision(std::numeric_limits<T>::digits10);
		 file.open (filename);
		 for (std::vector<InitCondition<T>>::const_iterator m = mesh.begin(); m != mesh.end(); ++m)
		 {
			 auto ic = (*m);
			 file << ic.Argument << " " << ic.Value;
			 if (saveDerivatives)
				 file << " " << ic.Derivative;
			 file << endl;
		 }
         file.close();
	}

	template <bool Deriv, class T>
	void SaveFunction(const char* argument_name, const char* function_name, const char* filename, const std::vector<InitCondition<T>>& solution)
	{
		 std::ofstream file;
		 file.precision(std::numeric_limits<T>::digits10);
		 file.open (filename);
		 file << argument_name << " " << function_name << std::endl;
		 for (const auto ic : solution)
			 file << ic.Argument << " " << (Deriv ? ic.Derivative : ic.Value) << std::endl;

		 file.close();
	}

};

#endif