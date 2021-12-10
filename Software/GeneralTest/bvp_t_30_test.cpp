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
#include <string>

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef double numTypeMp;

namespace GeneralTest
{
	void run_hybrid_cannon_test(const ProblemAbstract<numTypeMp>& problem, 
		const numTypeMp u_init, const numTypeMp u_final, 
		const numTypeMp deriv_left, const numTypeMp deriv_right,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference, 
		const numTypeMp critical_d_value = 1)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		numTypeMp h = 0.001;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc, critical_d_value);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[u_init, u_final](const InitCondition<numTypeMp>& ic) { return sgn(ic.Argument - (ic.Value - u_init)/(u_final-u_init)); };
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, u_init, u_final, deriv_left, deriv_right, evalFunc), L"Bisection failed");

		const auto knots = cannon.GetKnotVectorStreight();

		const auto diff_deriv_initial = abs(knots.begin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(knots.rbegin()->Derivative - deriv_final_reference);

		Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
		Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");
	}

	template<class T>
	void run_hybrid_multiple_shooting_test(const ProblemAbstract<T>& problem, 
		const T u_init, const T u_final, 
		const T deriv_left, const T deriv_right,
		const T deriv_initial_reference, const T deriv_final_reference, 
		const T critical_d_value = 1)
	{
		std::function<bool(const InitCondition<T>&)> checkFunc = 
			[](const InitCondition<T>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		T h = 0.01;

		HybridCannon<T> cannon(problem, h, 100*std::numeric_limits<T>::epsilon(), checkFunc, critical_d_value);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[u_init, u_final](const InitCondition<T>& ic) { return sgn(ic.Argument - (ic.Value - u_init)/(u_final-u_init)); };
		BisectionComponent<T> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, u_init, u_final, deriv_left, deriv_right, evalFunc), L"Bisection failed");

		auto init_guess = cannon.GetKnotVectorStreight();

		init_guess.rbegin()->Argument = 1;
		init_guess.rbegin()->Value = u_final;

		HybridMultipleShootingComponent<T> HMSComp(problem, critical_d_value);
		h = 0.001;
		bool succeeded;
		init_guess = HMSComp.Run(init_guess, h, succeeded);

		const auto deriv_final_diff = init_guess.rbegin()->Derivative - deriv_final_reference;
		const auto deriv_initial_diff = init_guess.begin()->Derivative - deriv_initial_reference;

		Logger::WriteMessage(("deriv_initial_diff = " + auxutils::ToString(deriv_initial_diff) + "\n").c_str());
		Logger::WriteMessage(("deriv_final_diff = " + auxutils::ToString(deriv_final_diff) + "\n").c_str());

		Assert::IsTrue(abs(deriv_initial_diff) < 12*h*h, L"Too big difference for the initial value of the derivative");
		Assert::IsTrue(abs(deriv_final_diff) < 12*h*h, L"Too big difference for the initial value of the derivative");

		Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));
	}


	template<class T>
	void run_hybrid_multiple_shooting_test_chasing(const ProblemAbstract<T>& problem, 
		const T u_init, const T u_final, 
		const T deriv_left, const T deriv_right,
		const T deriv_initial_reference, const T deriv_final_reference, 
		const T critical_d_value = 1)
	{
		std::function<bool(const InitCondition<T>&)> checkFunc = 
			[](const InitCondition<T>& ic) { return (abs(ic.Value) < 3) && (abs(ic.Argument) < 3); };

		T h = 0.001;

		HybridCannon<T> cannon(problem, h, 100*std::numeric_limits<T>::epsilon(), checkFunc, critical_d_value);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[u_init, u_final](const InitCondition<T>& ic) { return sgn(ic.Argument - (ic.Value - u_init)/(u_final-u_init)); };
		BisectionComponent<T> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, u_init, u_final, deriv_left, deriv_right, evalFunc), L"Bisection failed");

		auto init_guess = cannon.GetKnotVectorStreight();

		auxutils::SaveToFile(init_guess, "E:\\Research\\init_guess.txt");

		init_guess.rbegin()->Argument = 1;
		init_guess.rbegin()->Value = 1.5;

		std::vector<std::pair<T,T>> params(11);
		params[0] = std::pair<T,T>(0.01, 0.001);
		params[1] = std::pair<T,T>(0.006, 0.001);
		params[2] = std::pair<T,T>(0.005, 0.001);
		params[3] = std::pair<T,T>(0.004, 0.001);
		params[4] = std::pair<T,T>(0.003, 0.001);
		params[5] = std::pair<T,T>(0.002, 0.001);
		params[6] = std::pair<T,T>(0.0015, 0.001);
		params[7] = std::pair<T,T>(0.00125, 0.001);
		params[8] = std::pair<T,T>(0.001125, 0.001);
		params[9] = std::pair<T,T>(0.001, 0.001);
		params[10] = std::pair<T,T>(0.001, 0.0001);

		for (auto param : params )
		{
			const auto lambda = param.first;

			Bvp_t_30<T> problem_ch(lambda);


			//init_guess.rbegin()->Value = u_final;
			HybridMultipleShootingComponent<T> HMSComp(problem_ch, critical_d_value);
			h = param.second;
			bool succeeded;
			init_guess = HMSComp.Run(init_guess, h, succeeded);

			std::string filename = std::string("E:\\Research\\corrections_") + auxutils::ToString(lambda) + "_" + auxutils::ToString(h) + ".txt";
			auxutils::SaveToFile(HMSComp.GetCorrectionMgnitudes(), filename.c_str());
		}

		auxutils::SaveToFile(init_guess, "E:\\Research\\final_mesh.txt");

		//std::vector<InitCondition<double>> solution_double_precision(init_guess.size());

		//for (int point_id = 0; point_id < solution_double_precision.size(); point_id++)
		//{
		//	solution_double_precision[point_id].Argument      = init_guess[point_id].Argument.convert_to<double>();
		//	solution_double_precision[point_id].Value         = init_guess[point_id].Value.convert_to<double>();
		//	solution_double_precision[point_id].Derivative    = init_guess[point_id].Derivative.convert_to<double>();
		//	solution_double_precision[point_id].SecDerivative = init_guess[point_id].SecDerivative.convert_to<double>();
		//}

		//const auto solution_simplified = auxutils::simplify_polyline(init_guess, 0.1);

		//auxutils::SaveToFile(solution_double_precision, "E:\\Research\\solution_double_precision.txt");


		//const auto diff_deriv_initial = abs(result.begin()->Derivative - deriv_initial_reference);
		//const auto diff_deriv_final = abs(result.rbegin()->Derivative - deriv_final_reference);

		//Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
		//Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");

		//auxutils::SaveToFile(HMSComp.GetCorrectionMgnitudes(), "H:\\corrections.txt");

		//Assert::IsTrue(CheckQuadraticConvergenceOfNewtonMethd(HMSComp.GetCorrectionMgnitudes()), Message("Cannot confirm quadratic convergence rate"));
	}

	TEST_CLASS(bvp_t_30_test)
	{
public:
		
		TEST_METHOD(bvp_t_26_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_cannon_test(problem, numTypeMp(1), -numTypeMp(1)/3, -2, -10, -9.35131588439187988, -1.90926225716916762);
		}

		TEST_METHOD(bvp_t_30_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.001;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_cannon_test(problem, -numTypeMp(7)/6, numTypeMp(3)/2, 0.5, 1.5, 1, 1, 2);
		}

		TEST_METHOD(AutonomousGeneralDiagnosticsProblem_hybrid_cannot_test)
		{
			numTypeMp lambda = 0.04;
			AutonomousGeneralDiagnosticsProblem<numTypeMp> problem(lambda);

			run_hybrid_cannon_test(problem, numTypeMp(1), -numTypeMp(1)/3, -2, -10, -2.6726349376177122, -1.7518920136960743);
		}


		TEST_METHOD(bvp_t_26_hybrid_multiple_shooting_test)
		{
			numTypeMp lambda = 0.072;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_multiple_shooting_test<numTypeMp>(problem, numTypeMp(1), -numTypeMp(1)/3, -2, -10, -9.35131588439187988, -1.90926225716916762);
		}

		TEST_METHOD(bvp_t_30_hybrid_multiple_shooting_test)
		{
			numTypeMp lambda = 0.01;
			Bvp_t_30<numTypeMp> problem(lambda);

			run_hybrid_multiple_shooting_test<numTypeMp>(problem, -numTypeMp(7)/6, numTypeMp(3)/2, 0.5, 1.5, 1, 1, 2);
		}


		//TEST_METHOD(bvp_t_30_hybrid_multiple_shooting_chasing_test)
		//{
		//	numTypeMp lambda = 0.01;
		//	Bvp_t_30<numTypeMp> problem(lambda);

		//	run_hybrid_multiple_shooting_test_chasing<numTypeMp>(problem, -numTypeMp(7)/6, numTypeMp(3)/2, 0.5, 1.5, 1, 1, 2);
		//}

		TEST_METHOD(AutonimousGeneralDiagnosticsProblem_hybrid_multiple_shooting_test)
		{
			double lambda = 0.04;
			AutonomousGeneralDiagnosticsProblem<double> problem(lambda);

			run_hybrid_multiple_shooting_test<double>(problem, double(1), -double(1)/3, -2, -10, double(-2.6726349376177122), double(-1.7518920136960743));
		}
	};
}