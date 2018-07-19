#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include "UnitTestAux.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Problems/AutonomousNonUniformProblem.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef double numTypeMp;

namespace GeneralTest
{
	TEST_CLASS(NonUniformProblemTest)
	{
	private:
		std::vector<InitCondition<double>> GetGybridShootingResultForward(const double alpha, const double h)
		{
			 AutonomousNonUniformProblem<double> problem(alpha);
			 HybridCannon<double> thc(problem, h, 10*std::numeric_limits<double>::epsilon());

			 auto result = thc.Shoot(1, 4, sin(alpha), 10, alpha*cos(alpha));
			 return thc.GetKnotVectorStreight();
		}

		std::vector<InitCondition<double>> GetGybridShootingResultBackward(const double alpha, const double h)
		{
			 AutonomousNonUniformProblem<double> problem(alpha);
			 HybridCannon<double> thc(problem, h, 10*std::numeric_limits<double>::epsilon());

			 auto result = thc.Shoot(4, 1, sin(4*alpha), 10, alpha*cos(4*alpha));
			 return thc.GetKnotVectorStreight();
		}

		template <class T>
		void GetMaxValueAndDerivativeErrors(const std::vector<InitCondition<T>>& result, const T alpha, T& maxValueError, T& maxDerivError)
		{
			 std::vector<InitCondition<T>> errorVector(result.size()); 

			 std::transform(result.begin(), result.end(), errorVector.begin(), [=](InitCondition<T> ic) -> InitCondition<T> 
			 {
				 InitCondition<T> result;

				 result.Argument = ic.Argument;
				 result.Value = ic.Value - sin(alpha*ic.Argument);
				 result.Derivative = ic.Derivative - alpha*cos(alpha*ic.Argument);
				 result.SecDerivative = ic.SecDerivative + alpha*alpha*sin(alpha*ic.Argument);
				 
				 return result;
			 });

			 maxValueError = std::max_element(errorVector.begin(), errorVector.end(), 
				 [](const InitCondition<T>& ic1, const InitCondition<T>& ic2) {return abs(ic1.Value) < abs(ic2.Value); })->Value;

			 maxDerivError = std::max_element(errorVector.begin(), errorVector.end(), 
				 [](const InitCondition<T>& ic1, const InitCondition<T>& ic2) {return abs(ic1.Derivative) < abs(ic2.Derivative); })->Derivative;
		}

		template <class T>
		std::vector<InitCondition<T>> GetHybridBisectionResult(const T alpha, const T h)
		{
			 T argStart = T(1);
			 T uStart = sin(alpha*argStart);
			 T argTarget = T(4);
			 T uTarget = sin(argTarget*alpha);
			 T duLeft = alpha*cos(alpha*argStart) - T(0.1);
			 T duRight = alpha*cos(alpha*argStart) + T(0.1);
			 AutonomousNonUniformProblem<T> problem(alpha);

			 std::function<bool(const InitCondition<T>&)> checkFunc = 
				 [](const InitCondition<T>& ic) { return (abs(ic.Value) < 10) && (abs(ic.Argument) < 4); };

			 HybridCannon<T> cannon(problem, h, 10*std::numeric_limits<T>::epsilon(), checkFunc);

			 std::function<int(const InitCondition<T>&)> evalFunc = 
				 [=](const InitCondition<T>& ic) { return sgn(ic.Value - uTarget); };
			 BisectionComponent<T> bc(cannon);
			 Assert::IsTrue(bc.DerivativeBisectionGen(argStart, argTarget, uStart, uTarget, duLeft, duRight, evalFunc));

			 return cannon.GetKnotVectorStreight();
		}

		template <class T>
		bool CheckQuadraticConvergenceOfNewtonMethd(const std::vector<T>& successiveCorrections)
		{
			int numberOfAcceptableCorrections = 0;
			for (size_t index = successiveCorrections.size() - 1; index > 0; index--)
			{
				T prevCorrectionSquared = auxutils::sqr(successiveCorrections[index - 1]);
				if (2* prevCorrectionSquared >= successiveCorrections[index])
					numberOfAcceptableCorrections++;
				else if (index < successiveCorrections.size() - 1) //we can skip the very last correction because it can be not "clear" enough
					break;
			}

			return numberOfAcceptableCorrections >= 2;
		}

	public:
		
		TEST_METHOD(NonUniformHybridShootingTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;
			 auto result_h1 = GetGybridShootingResultForward(alpha, h1);
			 auto result_h2 = GetGybridShootingResultForward(alpha, h2);

			 double maxValueError_h1, maxDerivError_h1;

			 GetMaxValueAndDerivativeErrors(result_h1, alpha, maxValueError_h1, maxDerivError_h1);

			 double maxValueError_h2, maxDerivError_h2;

			 GetMaxValueAndDerivativeErrors(result_h2, alpha, maxValueError_h2, maxDerivError_h2);

			 double valueRatio = maxValueError_h2/maxValueError_h1;
			 double derivativeRatio = maxDerivError_h2/maxDerivError_h1;


			 Assert::IsTrue(valueRatio >= auxutils::sqr(h2/h1)*(1 - 0.002), Message("Too low ratio for values: " + auxutils::ToString(valueRatio)));
			 Assert::IsTrue(derivativeRatio >= auxutils::sqr(h2/h1)*(1 - 0.002), Message("Too low ratio for derivatives: " + auxutils::ToString(derivativeRatio)));
		}

		TEST_METHOD(NonUniformHybridShootingBackwardTest)
		{
			 double h1 = -0.0001;
			 double h2 = 10*h1;
			 double alpha = 2;
			 auto result_h1 = GetGybridShootingResultBackward(alpha, h1);
			 auto result_h2 = GetGybridShootingResultBackward(alpha, h2);

			 double maxValueError_h1, maxDerivError_h1;

			 GetMaxValueAndDerivativeErrors(result_h1, alpha, maxValueError_h1, maxDerivError_h1);

			 double maxValueError_h2, maxDerivError_h2;

			 GetMaxValueAndDerivativeErrors(result_h2, alpha, maxValueError_h2, maxDerivError_h2);

			 double valueRatio = maxValueError_h2/maxValueError_h1;
			 double derivativeRatio = maxDerivError_h2/maxDerivError_h1;

			 Assert::IsTrue(valueRatio >= auxutils::sqr(h2/h1), Message("Too low ratio for values: " + auxutils::ToString(valueRatio)));
			 Assert::IsTrue(derivativeRatio >= auxutils::sqr(h2/h1), Message("Too low ratio for derivatives: " + auxutils::ToString(derivativeRatio)));
		}

		TEST_METHOD(NonUniformHybridBisectionTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;

			 auto result_h1 = GetHybridBisectionResult(alpha, h1);
			 auto result_h2 = GetHybridBisectionResult(alpha, h2);

			 double maxValueError_h1, maxDerivError_h1;

			 GetMaxValueAndDerivativeErrors(result_h1, alpha, maxValueError_h1, maxDerivError_h1);

			 double maxValueError_h2, maxDerivError_h2;

			 GetMaxValueAndDerivativeErrors(result_h2, alpha, maxValueError_h2, maxDerivError_h2);

			 double valueRatio = maxValueError_h2/maxValueError_h1;
			 double derivativeRatio = maxDerivError_h2/maxDerivError_h1;

			 Assert::IsTrue(valueRatio >= auxutils::sqr(h2/h1), Message("Too low ratio for values: " + auxutils::ToString(valueRatio)));
			 Assert::IsTrue(derivativeRatio >= auxutils::sqr(h2/h1)*(1-0.002), Message("Too low ratio for derivatives: " + auxutils::ToString(derivativeRatio)));
		}

		TEST_METHOD(NonUniformHybridMultipleShootingTest)
		{
			 numTypeMp alpha = numTypeMp(2)/1;

			 auto init_guess = auxutils::ReadFromFile<InitCondition<numTypeMp>>("f:\\NonUniformAutonomousProblemInitGuess.txt");

			 AutonomousNonUniformProblem<numTypeMp> problem(alpha);
			 HybridMultipleShootingComponent<numTypeMp> HMSComp(problem);

			 bool succeeded;
			 numTypeMp h1 = numTypeMp(1)/10000;
			 std::vector<InitCondition<numTypeMp>> result_h1 = HMSComp.Run(init_guess, h1, succeeded);

			 Assert::IsTrue(succeeded, Message("Hybrid multiple shooting failed for h1"));

			 Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));

			 numTypeMp h2 = 10*h1;
			 std::vector<InitCondition<numTypeMp>> result_h2 = HMSComp.Run(init_guess, h2, succeeded);

			 Assert::IsTrue(succeeded, Message("Hybrid multiple shooting failed for h2"));

			 Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));

			 numTypeMp maxValueError_h1, maxDerivError_h1;

			 GetMaxValueAndDerivativeErrors(result_h1, alpha, maxValueError_h1, maxDerivError_h1);

			 numTypeMp maxValueError_h2, maxDerivError_h2;

			 GetMaxValueAndDerivativeErrors(result_h2, alpha, maxValueError_h2, maxDerivError_h2);

			 numTypeMp valueRatio = maxValueError_h2/maxValueError_h1;
			 numTypeMp derivativeRatio = maxDerivError_h2/maxDerivError_h1;

			 Assert::IsTrue(valueRatio >= auxutils::sqr(h2/h1)*(1-0.012), Message("Too low ratio for values: " + auxutils::ToString(valueRatio)));
			 Assert::IsTrue(derivativeRatio >= auxutils::sqr(h2/h1)*(1-0.012), Message("Too low ratio for derivatives: " + auxutils::ToString(derivativeRatio)));
    	}

	};
}