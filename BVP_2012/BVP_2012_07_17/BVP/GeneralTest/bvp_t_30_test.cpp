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
#include "../BVP/Problems/AutonomousGeneralDiagnosticsProblem.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef double numTypeMp;

namespace GeneralTest
{
	void run_hybrid_cannon_test(const ProblemAbstract<numTypeMp>& problem, 
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		numTypeMp h = 0.001;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[](const InitCondition<numTypeMp>& ic) { return sgn(ic.Argument*4 +  3*(ic.Value - 1)); };
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 1.0, -1.0/3, -2, -10, evalFunc), L"Bisection failed");

		const auto knots = cannon.GetKnotVectorStreight();

		const auto diff_deriv_initial = abs(knots.begin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(knots.rbegin()->Derivative - deriv_final_reference);

		Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
		Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");
	}

	template<class T>
	void run_hybrid_multiple_shooting_test(const ProblemAbstract<T>& problem, 
		const T deriv_initial_reference, const T deriv_final_reference)
	{
		std::function<bool(const InitCondition<T>&)> checkFunc = 
			[](const InitCondition<T>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		T h = 0.01;

		HybridCannon<T> cannon(problem, h, 0.001, checkFunc);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[](const InitCondition<T>& ic) { return sgn(ic.Argument*4 +  3*(ic.Value - 1)); };
		BisectionComponent<T> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 1.0, -1.0/3, -2, -10, evalFunc), L"Bisection failed");

		auto init_guess = cannon.GetKnotVectorStreight();

		init_guess.rbegin()->Argument = 1.0;
		init_guess.rbegin()->Value = -1.0/3;

		HybridMultipleShootingComponent<T> HMSComp(problem);
		h = 0.0001;
		bool succeeded;
		std::vector<InitCondition<T>> result = HMSComp.Run(init_guess, h, succeeded);

		const auto diff_deriv_initial = abs(result.begin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(result.rbegin()->Derivative - deriv_final_reference);

		Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
		Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");

		//auxutils::SaveToFile(HMSComp.GetCorrectionMgnitudes(), "H:\\corrections.txt");

		Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));
	}

	TEST_CLASS(bvp_t_30_test)
	{
public:
		
		TEST_METHOD(bvp_t_30_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_cannon_test(problem, -9.35131588439187988, -1.90926225716916762);
		}

		TEST_METHOD(AutonomousGeneralDiagnosticsProblem_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.04;
			AutonomousGeneralDiagnosticsProblem<numTypeMp> problem(lambda);

			run_hybrid_cannon_test(problem, -2.6726349376177122, -1.7518920136960743);
		}


		TEST_METHOD(bvp_t_30_hybrid_multiple_shooting_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_multiple_shooting_test(problem, -9.35131588439187988, -1.90926225716916762);
		}

		TEST_METHOD(AutonimousGeneralDiagnosticsProblem_hybrid_multiple_shooting_test)
		{
			double lambda = 0.04;
			AutonomousGeneralDiagnosticsProblem<double> problem(lambda);

			run_hybrid_multiple_shooting_test(problem, double(-2.6726349376177122), double(-1.7518920136960743));
		}
	};
}