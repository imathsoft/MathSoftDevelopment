#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include "UnitTestAux.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Problems/bvp_t_21.h"
#include "../BVP/Problems/TroeschProblem.h"
#include "../BVP/Problems/AutonomousGeneralDiagnosticsProblem.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//typedef float_50_noet numTypeMp;
typedef double numTypeMp;

namespace GeneralTest
{
	std::vector<InitCondition<numTypeMp>> get_init_guess(const ProblemAbstract<numTypeMp>& problem, const numTypeMp& u_target,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (ic.Value <= 1) && (ic.Value >= 0) && (ic.Argument >= 0) && (ic.Argument <= 1); };

		numTypeMp h = -0.01;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[u_target](const InitCondition<numTypeMp>& ic) 
		{
			if (ic.Argument > 0 && ic.Value >= 1)
				return 1;

			return sgn(ic.Value -  1.0); 
		};
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(1.0, 0, u_target, 1.0, 0, -0.1, evalFunc), L"Bisection failed");


		auto knots = cannon.GetKnotVectorStreight();

		knots.rbegin()->Argument = 0;
		knots.rbegin()->Value = 1.0;

		return knots;
	}

	std::vector<InitCondition<numTypeMp>> run_multiple_shooting(const ProblemAbstract<numTypeMp>& problem, const std::vector<InitCondition<numTypeMp>>& init_guess, const numTypeMp& u_target,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference, const numTypeMp h)
	{

		HybridMultipleShootingComponent<numTypeMp> HMSComp(problem);
		bool succeeded;
		std::vector<InitCondition<numTypeMp>> result = HMSComp.Run(init_guess, h, succeeded);

		const auto diff_deriv_initial = abs(result.rbegin()->Derivative - deriv_initial_reference);
		const auto diff_deriv_final = abs(result.begin()->Derivative - deriv_final_reference);

		return result;
	}

	void run_hybrid_cannon_test(const ProblemAbstract<numTypeMp>& problem, const numTypeMp& u_target,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (ic.Value <= 1) && (ic.Value >= 0) && (ic.Argument >= 0) && (ic.Argument <= 1); };

		numTypeMp h = -0.01;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

		//std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
		//	[u_target](const InitCondition<numTypeMp>& ic) 
		//{
		//	if(ic.Derivative > 1 && ic.Value > 1)
		//		return 1;

		//	if (ic.Argument < 0.999)
		//		return -1;

		//		return sgn(ic.Value -  u_target); 
		//};
		//BisectionComponent<numTypeMp> bc(cannon);
		//Assert::IsTrue(bc.DerivativeBisectionGen(0, 1.0, 1.0, u_target, -90, -110, evalFunc), L"Bisection failed");

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[u_target](const InitCondition<numTypeMp>& ic) 
		{
			if (ic.Argument > 0 && ic.Value >= 1)
				return 1;

			return sgn(ic.Value -  1.0); 
		};
		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(1.0, 0, u_target, 1.0, 0, -0.1, evalFunc), L"Bisection failed");


		auto knots = cannon.GetKnotVectorStreight();

		auto diff_deriv_initial = abs(knots.rbegin()->Derivative - deriv_initial_reference);
		auto diff_deriv_final = abs(knots.begin()->Derivative - deriv_final_reference);


		//knots.rbegin()->Argument = 1.0;
		//knots.rbegin()->Value = u_target;

		knots.rbegin()->Argument = 0;
		knots.rbegin()->Value = 1.0;

		HybridMultipleShootingComponent<numTypeMp> HMSComp(problem);
		h = 0.0001;
		bool succeeded;
		std::vector<InitCondition<numTypeMp>> result = HMSComp.Run(knots, h, succeeded);

		diff_deriv_initial = abs(result.rbegin()->Derivative - deriv_initial_reference);
		diff_deriv_final = abs(result.begin()->Derivative - deriv_final_reference);
		//Assert::IsTrue(diff_deriv_initial < 20*h*h, L"Too big deviation from the reference");
		//Assert::IsTrue(diff_deriv_final < 4*h*h, L"Too big deviation from the reference");
	}

	std::vector<InitCondition<numTypeMp>> get_init_guess_forward(const ProblemAbstract<numTypeMp>& problem, const numTypeMp& u_target,
		const numTypeMp deriv_initial_reference, const numTypeMp deriv_final_reference)
	{
		std::function<bool(const InitCondition<numTypeMp>&)> checkFunc = 
			[](const InitCondition<numTypeMp>& ic) { return (ic.Value <= 1) && (ic.Value >= -0.1) && (ic.Argument >= 0) && (ic.Argument <= 1); };

		numTypeMp h = 0.001;

		HybridCannon<numTypeMp> cannon(problem, h, 10*std::numeric_limits<numTypeMp>::epsilon(), checkFunc);

		std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
			[u_target](const InitCondition<numTypeMp>& ic) 
		{
			if(ic.Derivative > 1 && ic.Value >= 1)
				return 1;

			if (ic.Argument < 0.999)
				return -1;

				return sgn(ic.Value -  u_target); 
		};

		BisectionComponent<numTypeMp> bc(cannon);
		Assert::IsTrue(bc.DerivativeBisectionGen(0, 1.0, 1.0, u_target, deriv_initial_reference-5 , 
			deriv_initial_reference + 4, evalFunc), L"Bisection failed");

		auto knots = cannon.GetKnotVectorStreight();

		auto diff_deriv_initial = abs(knots.rbegin()->Derivative - deriv_initial_reference);
		auto diff_deriv_final = abs(knots.begin()->Derivative - deriv_final_reference);


		Logger::WriteMessage((auxutils::ToString(knots.rbegin()->Argument) + "\n").c_str());
		Logger::WriteMessage((auxutils::ToString(knots.rbegin()->Value) + "\n").c_str());
		Logger::WriteMessage((auxutils::ToString(knots.rbegin()->Derivative) + "\n").c_str());

		knots.rbegin()->Argument = 1.0;
		knots.rbegin()->Value = u_target;

		return knots;
	}


	TEST_CLASS(bvp_t_21_test)
	{
public:
		
		TEST_METHOD(bvp_t_21_hybrid_cannot_test)
		{
			numTypeMp lambda = numTypeMp(0.00001);
			Bvp_t_21<numTypeMp> problem_init_guess(lambda);

			const auto init_guess = get_init_guess(problem_init_guess, exp(-1/sqrt(lambda)), -1/sqrt(lambda), (-1/sqrt(lambda))*exp(-1/sqrt(lambda)));

			Bvp_t_21<numTypeMp> problem(lambda);

			const auto u = [lambda](const numTypeMp x) { return exp(-x/auxutils::Sqrt(lambda)); };
			const auto u_prime = [lambda](const numTypeMp x) { return -exp(-x/auxutils::Sqrt(lambda))/auxutils::Sqrt(lambda); };
			const auto u_prime_prime = [lambda](const numTypeMp x) { return exp(-x/auxutils::Sqrt(x))/lambda; };

			const auto x = [lambda](const numTypeMp u) { return -auxutils::Sqrt(lambda)*log(u);};
			const auto x_prime = [lambda](const numTypeMp u) { return -auxutils::Sqrt(lambda)/u;};
			const auto x_prime_prime = [lambda](const numTypeMp u) { return auxutils::Sqrt(lambda)/(u*u);};

			const numTypeMp h = 0.0001;
			const auto result = run_multiple_shooting(problem, init_guess, u(numTypeMp(1)), u_prime(numTypeMp(0)), u_prime(numTypeMp(1)), h);

			std::vector<InitCondition<numTypeMp>> diff_S;
			std::vector<InitCondition<numTypeMp>> diff_I;
			int N_I = 0;
			int N_S = 0;
			InitCondition<numTypeMp> c;

			numTypeMp k_s_0 = 0, k_s_1 = 0, k_i_0 = 0, k_i_1 = 0;
			numTypeMp k_s_0_x, k_s_1_x, k_i_0_u, k_i_1_u;

			for (int knot_id = 0; knot_id < result.size(); knot_id++)
			{
				const auto& data = result[knot_id];

				if (abs(data.Derivative) > 1)
				{
					if ((k_i_0 < abs(data.Argument - x(data.Value))/(h*h)))
					{
						k_i_0 = abs(data.Argument - x(data.Value))/(h*h);
						k_i_0_u = data.Value;
					}

					if (k_i_1 < abs(1/data.Derivative - x_prime(data.Value))/(h*h))
					{
						k_i_1 = abs(1/data.Derivative - x_prime(data.Value))/(h*h);
						k_i_1_u = data.Value;
					}

					N_I++;
					if(knot_id == 0 || abs(result[knot_id-1].Derivative) <= 1 || 
						knot_id == result.size() - 1 || abs(result[knot_id+1].Derivative) <= 1 || 
						floor(10*result[knot_id-1].Value) != floor(10*result[knot_id].Value))
					{

						if (knot_id > 0 && abs(result[knot_id-1].Derivative) <= 1)
						{
							c = data;
						}

						const auto u = data.Value;
						InitCondition<numTypeMp> diff;
						diff.Value = u;
						diff.Derivative = abs(1/data.Derivative - x_prime(u))/(h*h);
						diff.SecDerivative = 0;
						diff.Argument = abs(data.Argument - x(u))/(h*h);
						diff_I.push_back(diff);
					}
				} else
				{
					if (k_s_0 < abs(data.Value - u(data.Argument))/(h*h))
					{
						k_s_0 = abs(data.Value - u(data.Argument))/(h*h);
						k_s_0_x = data.Argument;
					}

					if (k_s_1 < abs(data.Derivative - u_prime(data.Argument))/(h*h))
					{
						k_s_1 = abs(data.Derivative - u_prime(data.Argument))/(h*h);
						k_s_1_x = data.Argument;
					}

					N_S++;
					if(knot_id == 0 || abs(result[knot_id-1].Derivative) > 1 || 
						knot_id == result.size() - 1 || abs(result[knot_id+1].Derivative) > 1 || 
						floor(10*result[knot_id-1].Argument) != floor(10*result[knot_id].Argument))
					{
						const auto x = data.Argument;
						InitCondition<numTypeMp> diff;
						diff.Value = abs(data.Value - u(x))/(h*h);
						diff.Derivative = abs(data.Derivative - u_prime(x))/(h*h);
						diff.SecDerivative = 0;
						diff.Argument = x;
						diff_S.push_back(diff);
					}
				}
			}

			Logger::WriteMessage(("lambda = " + auxutils::ToString(lambda) + "\n").c_str());
			Logger::WriteMessage(("N_S = " + auxutils::ToString(N_S) + "\n").c_str());
			Logger::WriteMessage(("N_I = " + auxutils::ToString(N_I) + "\n").c_str());
			Logger::WriteMessage(("c = " + auxutils::ToString(c.Argument) + "\n").c_str());
			Logger::WriteMessage(("u(c) = " + auxutils::ToString(c.Value) + "\n").c_str());
			Logger::WriteMessage(("u'(c) = " + auxutils::ToString(c.Derivative) + "\n").c_str());
			Logger::WriteMessage(("k_s_0 = " + auxutils::ToString(k_s_0) + " k_s_0_x = " + auxutils::ToString(k_s_0_x) + "\n").c_str());
			Logger::WriteMessage(("k_s_1 = " + auxutils::ToString(k_s_1) + " k_s_1_x = " + auxutils::ToString(k_s_1_x) + "\n").c_str());
			Logger::WriteMessage(("k_i_0 = " + auxutils::ToString(k_i_0) + " k_i_0_u = " + auxutils::ToString(k_i_0_u) + "\n").c_str());
			Logger::WriteMessage(("k_i_1 = " + auxutils::ToString(k_i_1) + " k_i_1_u = " + auxutils::ToString(k_i_1_u) + "\n").c_str());

			Assert::IsTrue(k_s_0 < 20, L"Too big value of the kappa coefficient for the solution of the straight problem");
			Assert::IsTrue(k_i_0 < 20, L"Too big value of the rappa coefficient for the solution of the inverse problem");
			Assert::IsTrue(k_s_1 < 6000, L"Too big value of the kappa coefficient for the derivative of the straight problem");
			Assert::IsTrue(k_i_1 < 3000, L"Too big value of the kappa coefficient for the derivative of the iverse problem");

			//auxutils::SaveToFile(result, "E:\\Research\\bvp_t_21.txt");
			//auxutils::SaveToFile(diff_I, "E:\\Research\\bvp_t_21_diff_I.txt");
			//auxutils::SaveToFile(diff_S, "E:\\Research\\bvp_t_21_diff_S.txt");
		}
	};
}