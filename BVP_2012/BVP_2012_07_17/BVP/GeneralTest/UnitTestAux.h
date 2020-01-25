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
#include "../BVP/Problems/ProblemAbstract.h"

namespace UnitTestAux
{
	template <class T>
	bool CheckQuadraticConvergenceOfNewtonMethd(const std::vector<T>& successiveCorrections)
	{
		int numberOfAcceptableCorrections = 0;
		for (size_t index = successiveCorrections.size() - 1; index > 0; index--)
		{
			T prevCorrectionSquared = auxutils::sqr(successiveCorrections[index - 1]);
			if (3* prevCorrectionSquared >= successiveCorrections[index])
				numberOfAcceptableCorrections++;
			else if (index < successiveCorrections.size() - 1) //we can skip the very last correction because it can be not "clear" enough
				break;
		}

		return numberOfAcceptableCorrections >= 2;
	}

	///Method to calculate the deviation of the given set of knots to the given "exact solution"
	template<class T>
	T CalcDeviationFromExactSolution(std::vector<InitCondition<T>> knots, 
		std::function<T(const T&)> exactSolution, bool checkDerivative = false)
	{
		T deviation = 0;
		for (std::vector<InitCondition<T>>::const_iterator iter = knots.begin(); iter!= knots.end(); ++iter)
		{
			InitCondition<T> knot = (*iter);
			if (checkDerivative)
				deviation = max(deviation, abs(exactSolution(knot.Argument) - knot.Derivative));
			else
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
			maxDistance = max(maxDistance, (prevKnot-nextKnot).NormSquaredPartial<false>());
			prevKnot = nextKnot;
		}

		return maxDistance;
	}

	///A "Standard test" for a problem whose solution is exp(sin(x));
	///The "problem" can be either autonomous or not
	template<class T,class P>
	void StandardOscillatinProblemMultipleShoothingTest(P problem)
	{
		T preH = 0.1;
		T finalH = 0.001;
		T targetValue = 0.5804096620;
		T targetArgument = 10;

		std::function<bool(const InitCondition<T>&)> checkFunc = 
			[=](const InitCondition<T>& ic) { return (abs(ic.Value) <= targetArgument) 
			&& (abs(ic.Argument) <= targetArgument); };

		HybridCannon<T> cannon(problem, preH, preH/10.0, checkFunc);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[=](const InitCondition<T>& ic) { return sgn(ic.Value - targetValue); };
		BisectionComponent<T> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, targetArgument, 1.0, 100, 0.95, 1.01, evalFunc));

		auto knots = cannon.GetKnotVectorStreight();
		knots[knots.size() - 1].Value = targetValue;
		knots[knots.size() - 1].Argument = targetArgument;

		HybridMultipleShootingComponent<T> HMSComp(problem);
		HMSComp.MaxNumberOfNewtonIterations = 20;

		bool succeeded;
		std::vector<InitCondition<T>> solution = HMSComp.Run(knots, finalH, succeeded);

		Assert::IsTrue(succeeded, Message("Multiple shoothing method has failed"));

		T maxDev = CalcDeviationFromExactSolution<T>(solution, 
	    [](const T& u){ return exp(sin(u)); });

		T maxDerivativeDev = CalcDeviationFromExactSolution<T>(solution, 
			[](const T& u){ return cos(u)*exp(sin(u)); }, true);

		T vaxDistanceBetweenKnots = CalcMaxSquaredDistanceBetweenNeighbourKnots(solution);
		Assert::IsTrue(vaxDistanceBetweenKnots < 2*finalH*finalH, 
			Message("Too big maximal distance between neighbour knots " + 
			auxutils::ToString(vaxDistanceBetweenKnots)));
		Assert::IsTrue(maxDev < finalH*finalH, 
			Message("Too big deviation, maxDev= " + auxutils::ToString(maxDev)));
		Assert::IsTrue(maxDerivativeDev < 2*finalH*finalH, 
			Message("Too big deviation between derivatives, maxDev= " + 
			auxutils::ToString(maxDev)));
		Assert::IsTrue(solution.size() > targetArgument/finalH, 
			Message("Too few knot in the solution vector, " + 
			auxutils::ToString(solution.size())));
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