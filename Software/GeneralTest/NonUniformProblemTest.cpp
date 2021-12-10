#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include "UnitTestAux.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Problems/AutonomousNonUniformProblem.h"
#include "../BVP/Problems/NonAutonomousNonUniformProblem.h"
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
		template <class T>
		std::vector<InitCondition<T>> GetGybridShootingResultForward(const T alpha, const T h, const bool nonAutonomousProblem)
		{
			ProblemAbstract<T>* problem; 
			if (nonAutonomousProblem)  
				problem = new NonAutonomousNonUniformProblem<T>(alpha); 
			else 
				problem = new AutonomousNonUniformProblem<T>(alpha);

			HybridCannon<T> thc(*problem, h, 10*std::numeric_limits<T>::epsilon());

			auto result = thc.Shoot(1, 4, sin(alpha), 10, alpha*cos(alpha));

			delete problem;
			return thc.GetKnotVectorStreight();
		}

		template <class T>
		std::vector<InitCondition<T>> GetGybridShootingResultBackward(const T alpha, const T h, const bool nonAutonomousProblem)
		{
			ProblemAbstract<T>* problem; 
			if (nonAutonomousProblem)  
				problem = new NonAutonomousNonUniformProblem<T>(alpha); 
			else 
				problem = new AutonomousNonUniformProblem<T>(alpha);

			 HybridCannon<T> thc(*problem, h, 10*std::numeric_limits<T>::epsilon());

			 auto result = thc.Shoot(4, 1, sin(4*alpha), 10, alpha*cos(4*alpha));
			 delete problem;
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
		std::vector<InitCondition<T>> GetHybridBisectionResult(const T alpha, const T h, const bool nonAutonomousProblem)
		{
			 T argStart = T(1);
			 T uStart = sin(alpha*argStart);
			 T argTarget = T(4);
			 T uTarget = sin(argTarget*alpha);
			 T duLeft = alpha*cos(alpha*argStart) - T(0.1);
			 T duRight = alpha*cos(alpha*argStart) + T(0.1);

			ProblemAbstract<T>* problem; 
			if (nonAutonomousProblem)  
				problem = new NonAutonomousNonUniformProblem<T>(alpha); 
			else 
				problem = new AutonomousNonUniformProblem<T>(alpha);

			 std::function<bool(const InitCondition<T>&)> checkFunc = 
				 [](const InitCondition<T>& ic) { return (abs(ic.Value) < 10) && (abs(ic.Argument) < 4); };

			 HybridCannon<T> cannon(*problem, h, 10*std::numeric_limits<T>::epsilon(), checkFunc);

			 std::function<int(const InitCondition<T>&)> evalFunc = 
				 [=](const InitCondition<T>& ic) { return sgn(ic.Value - uTarget); };
			 BisectionComponent<T> bc(cannon);
			 Assert::IsTrue(bc.DerivativeBisectionGen(argStart, argTarget, uStart, uTarget, duLeft, duRight, evalFunc));

			 delete problem;
			 return cannon.GetKnotVectorStreight();
		}

		template <class T>
		std::vector<InitCondition<T>> GetHybridMultipleShootingResult(const T alpha, const T h, const std::vector<InitCondition<T>> init_guess,
			const bool nonAutonomousProblem) 
		{
			ProblemAbstract<T>* problem; 
			if (nonAutonomousProblem)  
				problem = new NonAutonomousNonUniformProblem<T>(alpha); 
			else 
				problem = new AutonomousNonUniformProblem<T>(alpha);

			 HybridMultipleShootingComponent<numTypeMp> HMSComp(*problem);

			 bool succeeded;
			 std::vector<InitCondition<numTypeMp>> result = HMSComp.Run(init_guess, h, succeeded);

			 Assert::IsTrue(succeeded, Message("Hybrid multiple shooting failed."));

			 Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));

			 delete problem;
			 return result;
		}

		template <class T>
		void AssertQuadraticOrderOfApproximation(const std::vector<InitCondition<T>>& approximation_h1, const T h1, 
			const std::vector<InitCondition<T>>& approximation_h2, const T h2, const T alpha, 
			const T percentToleranceFunc, const T percentToleranceDeriv)
		{
			 T maxValueError_h1, maxDerivError_h1;

			 GetMaxValueAndDerivativeErrors(approximation_h1, alpha, maxValueError_h1, maxDerivError_h1);

			 T maxValueError_h2, maxDerivError_h2;

			 GetMaxValueAndDerivativeErrors(approximation_h2, alpha, maxValueError_h2, maxDerivError_h2);

			 T valueRatio = maxValueError_h2/maxValueError_h1;
			 T derivativeRatio = maxDerivError_h2/maxDerivError_h1;

			 T h2_to_h1_ratio_squared = auxutils::sqr(h2/h1);

			 Assert::IsTrue(valueRatio >= h2_to_h1_ratio_squared*(1 - percentToleranceFunc), Message("Too low ratio for values: " + auxutils::ToString(valueRatio)));
			 Assert::IsTrue(derivativeRatio >= h2_to_h1_ratio_squared*(1 - percentToleranceDeriv), Message("Too low ratio for derivatives: " + auxutils::ToString(derivativeRatio)));
		}

	public:
		
		TEST_METHOD(NonUniformHybridShootingTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;
			 auto result_h1 = GetGybridShootingResultForward(alpha, h1, false);
			 auto result_h2 = GetGybridShootingResultForward(alpha, h2, false);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.002, 0.002);
		}

		TEST_METHOD(NonAutoNonUniformHybridShootingTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;
			 auto result_h1 = GetGybridShootingResultForward(alpha, h1, true);
			 auto result_h2 = GetGybridShootingResultForward(alpha, h2, true);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.001, 0.001);
		}


		TEST_METHOD(NonUniformHybridShootingBackwardTest)
		{
			 double h1 = -0.0001;
			 double h2 = 10*h1;
			 double alpha = 2;
			 auto result_h1 = GetGybridShootingResultBackward(alpha, h1, false);
			 auto result_h2 = GetGybridShootingResultBackward(alpha, h2, false);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.0, 0.0);
		}

		TEST_METHOD(NonAutoNonUniformHybridShootingBackwardTest)
		{
			 double h1 = -0.0001;
			 double h2 = 10*h1;
			 double alpha = 2;
			 auto result_h1 = GetGybridShootingResultBackward(alpha, h1, true);
			 auto result_h2 = GetGybridShootingResultBackward(alpha, h2, true);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.0, 0.009);
		}

		TEST_METHOD(NonUniformHybridBisectionTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;

			 auto result_h1 = GetHybridBisectionResult(alpha, h1, false);
			 auto result_h2 = GetHybridBisectionResult(alpha, h2, false);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.0, 0.002);
		}

		TEST_METHOD(NonAutoNonUniformHybridBisectionTest)
		{
			 double h1 = 0.0001;
			 double alpha = 2;
			 double h2 = h1*10;

			 auto result_h1 = GetHybridBisectionResult(alpha, h1, true);
			 auto result_h2 = GetHybridBisectionResult(alpha, h2, true);

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.0, 0.0003);
		}

		TEST_METHOD(NonUniformHybridMultipleShootingTest)
		{
			 auto init_guess = auxutils::ReadFromFile<InitCondition<numTypeMp>>("TestData\\NonUniformProblemInitGuess.txt");

			 numTypeMp alpha = numTypeMp(2)/1;
			 numTypeMp h1 = numTypeMp(1)/10000;
			 numTypeMp h2 = 10*h1;

			 auto result_h1 = GetHybridMultipleShootingResult(alpha, h1, init_guess, false); 
			 auto result_h2 = GetHybridMultipleShootingResult(alpha, h2, init_guess, false); 

			 AssertQuadraticOrderOfApproximation<numTypeMp>(result_h1, h1, result_h2, h2, alpha, 0.012, 0.012);
    	}

		TEST_METHOD(NonAutoNonUniformHybridMultipleShootingTest)
		{
			 auto init_guess = auxutils::ReadFromFile<InitCondition<numTypeMp>>("TestData\\NonUniformProblemInitGuess.txt");

			 numTypeMp alpha = numTypeMp(2)/1;
			 numTypeMp h1 = numTypeMp(1)/5000;
			 numTypeMp h2 = 10*h1;

			 auto result_h1 = GetHybridMultipleShootingResult(alpha, h1, init_guess, true); 
			 auto result_h2 = GetHybridMultipleShootingResult(alpha, h2, init_guess, true); 

			 AssertQuadraticOrderOfApproximation(result_h1, h1, result_h2, h2, alpha, 0.0, 0.0);
    	}
	};
}