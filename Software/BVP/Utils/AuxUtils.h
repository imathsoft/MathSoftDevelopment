#pragma once

#include <fstream>
#include <strstream>
#include <vector>
#include <Utils.h>
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\GradientVector.h"

namespace auxutils
{
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

	template<class U, int Size>
	inline GradientVector<U, Size> Exp(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Exp(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i]*result[0];

		return result;
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Sin(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Sin(gv[0]);
		auto deriv = Cos(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i] * deriv;

		return result;
	}

	template<class U, int Size>
	inline GradientVector<U, Size> Cos(const GradientVector<U, Size>& gv)
	{
		GradientVector<U, Size> result;
		result[0] = Cos(gv[0]);
		auto deriv = -Sin(gv[0]);
		for (int i = 1; i < Size; i++)
			result[i] = result[i] * deriv;

		return result;
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