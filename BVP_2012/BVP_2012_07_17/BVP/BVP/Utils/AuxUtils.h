#ifndef GUARD_AUXUTILS
#define GUARD_AUXUTILS

#include <mpreal.h>
#include <fstream>

#include <Eigen/Sparse>
#include <Eigen/MPRealSupport>

namespace auxutils
{
	inline double sqr(double d)
	{
		return d*d;
	}

	inline float sqr(float f)
	{
		return f*f;
	}

	inline mpfr::mpreal sqr(mpfr::mpreal mpr)
	{
		return mpfr::sqr(mpr);
	}

	inline int sqr(int intVal)
	{
		return intVal*intVal;
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

	void WriteToStream(std::ofstream& stream, float value );

	void WriteToStream(std::ofstream& stream, double value );

	///Method to write sparse matrix into a text file
	template <class T>
	void SaveToFile(Eigen::SparseMatrix<T> matrix, char* filename)
	{
		 ofstream file;
		 file.open (filename);
         file << matrix;
         file.close();
	}

	///Method to write sparse matrix into a text file
	template <class T>
	void SaveToFile(Eigen::Matrix<T, Eigen::Dynamic, 1> matrix, char* filename)
	{
		 ofstream file;
		 file.open (filename);
         file << matrix;
         file.close();
	}

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
}

#endif