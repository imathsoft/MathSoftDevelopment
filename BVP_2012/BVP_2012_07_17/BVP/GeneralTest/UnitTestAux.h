#ifndef GUARD_UNITTESTAUX
#define GUARD_UNITTESTAUX

#ifdef _DEBUG
	#ifdef WIN32
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/Win32/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/Win32/Release/mpfr.lib")
	#else
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/x64/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/x64/Release/mpfr.lib")
	#endif
#else
	#ifdef WIN32
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/Win32/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/Win32/Release/mpfr.lib")
	#else
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/x64/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/x64/Release/mpfr.lib")
	#endif
#endif

#include <string>
#include <functional>
#include <vector>
#include "../BVP/FunctionApproximation/InitialCondition.h"

namespace UnitTestAux
{
	///Method to calculate the deviation of the given set of knots to the given "exact solution"
	template<class T>
	T CalcDeviationFromExactSolution(std::vector<InitCondition<T>> knots, 
		std::function<T(const T&)> exactSolution)
	{
		T deviation = 0;
		for (std::vector<InitCondition<T>>::const_iterator iter = knots.begin(); iter!= knots.end(); ++iter)
		{
			InitCondition<T> knot = (*iter);
			deviation = max(deviation, abs(exactSolution(knot.Argument) - knot.Value));
		}

		return deviation;
	}

	///Method to calculate maximal squared distance between neighbour knots
	template<class T>
	T CalcMaxSquaredDistanceBetweenNeighbourKnots(std::vector<InitCondition<T>> knots)
	{
		T maxDistance = 0;
		if (knots.empty())
		{
			return maxDistance;
		}

		InitCondition<T> prevKnot, nextKnot;
		prevKnot = knots.front();

		for (std::vector<InitCondition<T>>::const_iterator 
			iter = std::next(knots.begin()); iter!= knots.end(); ++iter)
		{
			nextKnot = (*iter);
			maxDistance = max(maxDistance, (prevKnot-nextKnot).NormSquaredNaive());
			prevKnot = nextKnot;
		}

		return maxDistance;
	}



	static wchar_t* Message(const char* text)
	{
		size_t size = strlen(text) + 1;
		size_t convertedSize;
		wchar_t* wa = new wchar_t[size];
		mbstowcs_s(&convertedSize, wa, size, text, _TRUNCATE);
		return wa;
	}

	static wchar_t* Message(const std::string text)
	{
		return Message(text.c_str());
	}
};

#endif