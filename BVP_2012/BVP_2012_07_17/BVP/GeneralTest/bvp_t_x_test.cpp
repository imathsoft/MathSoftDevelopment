#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include "UnitTestAux.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Problems/bvp_t_x.h"
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
		const numTypeMp a, const numTypeMp b,
		const numTypeMp u_a, const numTypeMp u_b,
		const numTypeMp deriv_min, const numTypeMp deriv_max,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		numTypeMp h = 0.001;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[=](const InitCondition<numTypeMp>& ic) { return sgn((u_a - u_b)*(a - ic.Argument) - (u_a - ic.Value)*(a-b)); };
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(a, b, u_a, u_b, deriv_min, deriv_max, evalFunc), L"Bisection failed");

		const auto knots = cannon.GetKnotVectorStreight();

		const auto diff_deriv_initial = abs(knots.begin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(knots.rbegin()->Derivative - deriv_final_reference);

		Assert::IsTrue(diff_deriv_initial < h*h, L"Too big deviation from the reference");
		Assert::IsTrue(diff_deriv_final < h*h, L"Too big deviation from the reference");
	}

	template<class T>
	void run_hybrid_multiple_shooting_test(const ProblemAbstract<T>& problem, 
		const numTypeMp a, const numTypeMp b,
		const numTypeMp u_a, const numTypeMp u_b,
		const numTypeMp deriv_min, const numTypeMp deriv_max,
		const T deriv_initial_reference, const T deriv_final_reference)
	{
		std::function<bool(const InitCondition<T>&)> checkFunc = 
			[](const InitCondition<T>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		T h = 0.01;

		HybridCannon<T> cannon(problem, h, 0.001, checkFunc);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[=](const InitCondition<numTypeMp>& ic) { return sgn((u_a - u_b)*(a - ic.Argument) - (u_a - ic.Value)*(a-b)); };
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(a, b, u_a, u_b, deriv_min, deriv_max, evalFunc), L"Bisection failed");

		auto init_guess = cannon.GetKnotVectorStreight();

		init_guess.rbegin()->Argument = b;
		init_guess.rbegin()->Value = u_b;

		HybridMultipleShootingComponent<T> HMSComp(problem);
		h = 0.0001;
		bool succeeded;
		std::vector<InitCondition<T>> result = HMSComp.Run(init_guess, h, succeeded);

		const auto diff_deriv_initial = abs(result.begin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(result.rbegin()->Derivative - deriv_final_reference);

		Assert::IsTrue(diff_deriv_initial < h*h, L"Too big deviation from the reference");
		Assert::IsTrue(diff_deriv_final < h*h, L"Too big deviation from the reference");

		//auxutils::SaveToFile(HMSComp.GetCorrectionMgnitudes(), "H:\\corrections.txt");

		Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));
	}

	TEST_CLASS(bvp_t_x_test)
	{
		TEST_METHOD(bvp_t_x_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.01;
			Bvp_t_x<numTypeMp> problem(lambda);

			numTypeMp u_a = exp(-1*(-1+2)/auxutils::Sqrt(lambda));

			const auto sol = [lambda](const numTypeMp x) { return exp(x*(x+2)/auxutils::Sqrt(lambda)); };
			const auto sol_deriv = [lambda](const numTypeMp x) { return 2*(x+1)*exp(x*(x+2)/auxutils::Sqrt(lambda))/auxutils::Sqrt(lambda); };

			run_hybrid_cannon_test(problem, -1, 0, u_a, 1, 0, 1,  sol_deriv(-1), sol_deriv(0));
		}

		TEST_METHOD(bvp_t_x_hybrid_multiple_shooting_test)
		{
			numTypeMp lambda = 0.05;
			Bvp_t_x<numTypeMp> problem(lambda);

			numTypeMp u_a = exp(-1*(-1+2)/auxutils::Sqrt(lambda));

			const auto sol = [lambda](const numTypeMp x) { return exp(x*(x+2)/auxutils::Sqrt(lambda)); };
			const auto sol_deriv = [lambda](const numTypeMp x) { return 2*(x+1)*exp(x*(x+2)/auxutils::Sqrt(lambda))/auxutils::Sqrt(lambda); };

			run_hybrid_multiple_shooting_test(problem,-1, 0, u_a, 1, 0, 1, sol_deriv(-1), sol_deriv(0));
		}
	};
}