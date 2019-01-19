#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include "UnitTestAux.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Problems/bvp_t_30.h"
#include "../BVP/Problems/TroeschProblem.h"
#include "../BVP/Problems/NonAutonomousNonUniformProblem.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef double numTypeMp;

#define XI_BERNOULLI

namespace GeneralTest
{
	TEST_CLASS(bvp_t_30_test)
	{
public:
		
		TEST_METHOD(bvp_t_30_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
				[](const InitCondition<numTypeMp>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

			numTypeMp h = 0.001;

			HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

			std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
				[](const InitCondition<numTypeMp>& ic) { return sgn(ic.Argument*4 +  3*(ic.Value - 1)); };
			BisectionComponent<numTypeMp> bc(cannon);
			Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 1.0, -1.0/3, -9, -10, evalFunc), L"Bisection failed");

			const auto knots = cannon.GetKnotVectorStreight();

			const auto diff_deriv_initial = abs(knots.begin()->Derivative + 9.35131588439187988);
			const auto diff_deriv_final = abs(knots.rbegin()->Derivative + 1.90926225716916762);

			Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
			Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");
		}

		TEST_METHOD(bvp_t_30_hybrid_multiple_shooting_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
				[](const InitCondition<numTypeMp>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

			numTypeMp h = 0.1;

			HybridCannon<numTypeMp> cannon(problem, h, 0.001, checkFunc);

			std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
				[](const InitCondition<numTypeMp>& ic) { return sgn(ic.Argument*4 +  3*(ic.Value - 1)); };
			BisectionComponent<numTypeMp> bc(cannon);
			Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 1.0, -1.0/3, -9, -10, evalFunc), L"Bisection failed");

			auto init_guess = cannon.GetKnotVectorStreight();

			init_guess.rbegin()->Argument = 1.0;
			init_guess.rbegin()->Value = -1.0/3;

			HybridMultipleShootingComponent<numTypeMp> HMSComp(problem);
			h = 0.0001;
			bool succeeded;
			std::vector<InitCondition<numTypeMp>> result = HMSComp.Run(init_guess, h, succeeded);

			const auto diff_deriv_initial = abs(result.begin()->Derivative + 9.35131588439187988);
			const auto diff_deriv_final = abs(result.rbegin()->Derivative + 1.90926225716916762);

			Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
			Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");

			Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));
		}

	};

}