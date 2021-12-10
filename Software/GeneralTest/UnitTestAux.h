#ifndef GUARD_UNITTESTAUX
#define GUARD_UNITTESTAUX

#include <string>
#include <functional>
#include <vector>
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Problems/ProblemAbstract.h"
#include "../BVP/LinearAlgebra/Matrix.h"
#include <stdlib.h>

namespace UnitTestAux
{
	/// <summary>
	/// Returns pseudorandom double from [0, 1]
	/// </summary>
	inline double Random()
	{
		return (double)std::rand() / RAND_MAX;
	}

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

	/*
	* Assuming that each element of the given sequence can be approximately represented as M*q^{2^{n}}, 
	* the method calculates the "M" and "q" parameters via the least squares fitting approach as well as the maximal 
	* relative error of the approximation 
	*/
	template <class R>
	void ApproximateParametersOfQuadraticallyDecayingSequence(const std::vector<R>& sequence, R& M, R& q, R& max_rel_error)
	{
		LinAlg::Matrix<R, 2, 2> A{}; //Matrix of the linear system to solve
		LinAlg::Matrix<R, 2, 1> b{}; //The right hand side of the system to solve

		R power_of_two = R(1);
		for (auto val_id = 0; val_id < sequence.size(); val_id++)
		{
			A[0][0] += R(1);
			A[0][1] += power_of_two;

			A[1][1] += power_of_two * power_of_two;

			const auto log_val = auxutils::Log(sequence[val_id]);

			b[0][0] += log_val;
			b[1][0] += log_val * power_of_two;

			power_of_two *= R(2);
		}

		A[1][0] = A[0][1];

		const auto x = A.Inverse() * b;

		M = auxutils::Exp(x[0][0]);
		q = auxutils::Exp(x[1][0]);

		R power_of_q = q;
		max_rel_error = R(0);
		for (auto val_id = 0; val_id < sequence.size(); val_id++)
		{
			max_rel_error = std::max<R>(max_rel_error, auxutils::Abs((sequence[val_id] - M * power_of_q)/sequence[val_id]));
			power_of_q *= power_of_q;
		}
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