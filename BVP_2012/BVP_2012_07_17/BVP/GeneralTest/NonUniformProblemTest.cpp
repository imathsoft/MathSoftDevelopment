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

typedef float_50_noet numTypeMp;

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

		void GetMaxValueAndDerivativeErrors(const std::vector<InitCondition<double>>& result, const double alpha, double& maxValueError, double& maxDerivError)
		{
			 std::vector<InitCondition<double>> errorVector(result.size()); 

			 std::transform(result.begin(), result.end(), errorVector.begin(), [=](InitCondition<double> ic) -> InitCondition<double> 
			 {
				 InitCondition<double> result;

				 result.Argument = ic.Argument;
				 result.Value = ic.Value - sin(alpha*ic.Argument);
				 result.Derivative = ic.Derivative - alpha*cos(alpha*ic.Argument);
				 result.SecDerivative = ic.SecDerivative + alpha*alpha*sin(alpha*ic.Argument);
				 
				 return result;
			 });

			 maxValueError = std::max_element(errorVector.begin(), errorVector.end(), 
				 [](const InitCondition<double>& ic1, const InitCondition<double>& ic2) {return abs(ic1.Value) < abs(ic2.Value); })->Value;

			 maxDerivError = std::max_element(errorVector.begin(), errorVector.end(), 
				 [](const InitCondition<double>& ic1, const InitCondition<double>& ic2) {return abs(ic1.Derivative) < abs(ic2.Derivative); })->Derivative;
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

	};
}