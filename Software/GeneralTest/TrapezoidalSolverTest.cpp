#include "CppUnitTest.h"
#include <vector>
#include <string>
#include <exception>
#include <numeric>
#include "../BVP/LinearAlgebra/Matrix.h"
#include "../BVP/Systems/bvp_sys_factory.h"
#include "../BVP/Problems/References/TroeschReferences.h"
#include "../BVP/Systems/FiniteDifferenceSolver/trapezoidal_solver.h"
#include "../BVP/Systems/FiniteDifferenceSolver/ts_diagnostics.h"
#include "../BVP/Utils/AuxUtils.h"
#include "UnitTestAux.h"
#include "../BVP/Utils/LoggingUtils.h"
#include <map>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace LinAlg;

namespace GeneralTest
{
	TEST_CLASS(TrapezoidalSolverTest)
	{
		/// <summary>
		/// Method to generate initial guess for an iterative process of solving boundary value problem
		/// </summary>
		template <class R, int varCnt>
		std::vector<mesh_point<R, varCnt>> GenerateInitialGuess(const R& t0, const R& t1, std::function<mesh_point<R, varCnt>(const R& t)> generator, const int intervalCnt)
		{
			std::vector<mesh_point<R, varCnt>> result;
			const auto step = (t1 - t0) / intervalCnt;

			for (int interval_Id = 0; interval_Id <= intervalCnt; interval_Id++)
				result.push_back(generator(t0 + step*interval_Id));

			return result;
		}

		/// <summary>
		/// Returns randomly a value equals to either magnitude or -magnitude
		/// </summary>
		template <class R>
		static R rand_perturbation(const R& magnitude)
		{
			return std::rand() % 2 == 0 ? magnitude : -magnitude;
		}

		/// <summary>
		/// Introduces random perturbation into the given initial guess
		/// </summary>
		template <class R, int varCnt>
		std::vector<mesh_point<R, varCnt>> perturbate(const std::vector<mesh_point<R, varCnt>>& init_guess, const R max_perturbation_magnitude, const bool first_func_bc)
		{
			std::srand(1);

			std::vector<mesh_point<R, varCnt>> result = init_guess;

			for (int pt_id = 0; pt_id < result.size(); pt_id++)
			{
				result[pt_id][0] += R(rand_perturbation(max_perturbation_magnitude));
				result[pt_id][1] += R(rand_perturbation(max_perturbation_magnitude));
			}

			for (int iter = 0; iter < 5; iter++)
				for (int pt_id = 1; pt_id < result.size() - 1; pt_id++)
				{
					result[pt_id][0] = (result[pt_id - 1][0] + result[pt_id][0] + result[pt_id + 1][0])/R(3);
					result[pt_id][1] = (result[pt_id - 1][1] + result[pt_id][1] + result[pt_id + 1][1]) / R(3);
				}

			if (first_func_bc)
			{
				result[0][0] = init_guess[0][0];
				result[result.size() - 1][0] = init_guess[result.size() - 1][0];
			}
			else
			{
				result[0][1] = init_guess[0][1];
				result[result.size() - 1][1] = init_guess[result.size() - 1][1];
			}

			return result;
		}

		/// <summary>
		/// Method to check that the convergence rate of the iteration process is quadratic
		/// </summary>
		template <class R, class S>
		void CheckQuadraticConvergence(const trapezoidal_solver<R, S>& solver, const R& acceptable_approximation_error,
			const bool skip_silently_if_not_enough_corrections = false)
		{
			const auto corrections = solver.get_correcion_magnitudes();

			std::vector<R> chosen_corrections;

			for (int corr_id = corrections.size() - 2; corr_id >= 1 &&
				!corrections[corr_id].RefinementApplied &&
				corrections[corr_id].CorrectionMagnitude < R(0.1); corr_id--)
				chosen_corrections.insert(chosen_corrections.begin(), corrections[corr_id].CorrectionMagnitude);

			const auto enought_corrections_available = chosen_corrections.size() >= 3;

			if (!enought_corrections_available && skip_silently_if_not_enough_corrections)
				return;

			Assert::IsTrue(enought_corrections_available, L"Too few corrections to detect quadratic convergence rate");

			R M, q, max_rel_error;

			UnitTestAux::ApproximateParametersOfQuadraticallyDecayingSequence(chosen_corrections, M, q, max_rel_error);

			Logger::WriteMessage((std::string("Actual diviation of the convergence rate from the quadratic one is ") +
				auxutils::ToString(max_rel_error) + "\n").c_str());

			Assert::IsTrue(max_rel_error <= acceptable_approximation_error, L"Too high approximation error");
		}

		/// <summary>
		/// Returns maximal deviations between the aproximated solution and the reference fucntions
		/// </summary>
		template <class R, int eqCnt>
		std::array<R, eqCnt> get_max_deviations(const std::array<std::function<R(R)>, eqCnt>& reference,
			const std::vector<mesh_point<R, eqCnt + 1>>& solution)
		{
			std::array<R, eqCnt> result{};

			for (int func_id = 0; func_id < eqCnt; func_id++)
			{
				for (int pt_id = 0; pt_id < solution.size(); pt_id++)
					result[func_id] = std::max<R>(result[func_id], auxutils::Abs(reference[func_id](solution[pt_id][eqCnt]) - solution[pt_id][func_id]));
			}

			return result;
		}

		/// <summary>
		/// Returns average deviations between the aproximated solution and the reference fucntions
		/// </summary>
		template <class R, int eqCnt>
		std::array<R, eqCnt> get_average_deviations(const std::array<std::function<R(R)>, eqCnt>& reference,
			const std::vector<mesh_point<R, eqCnt + 1>>& solution)
		{
			std::array<R, eqCnt> result{};

			for (int func_id = 0; func_id < eqCnt; func_id++)
			{
				result[func_id] = std::accumulate(solution.begin(), solution.end(), R(0.0), [func_id, &reference](const R sum, const auto pt)
					{
						return auxutils::Abs(reference[func_id](pt[eqCnt]) - pt[func_id]) + sum;
					});

				result[func_id] /= solution.size();
			}

			return result;
		}

		/// <summary>
		/// Generic mathod to perform test with the Troesch's problem
		/// </summary>
		template <class R = double>
		void perform_Troesch_test(const bool use_reparametrization, const bool use_inversion, const int discretization = 1000,
			const simple_bvp<R, 2>& pr = bvp_sys_factory<R>::Troesch(3.0), const R& quadratic_convergence_tolerance = R(0.2))
		{
			const auto init_guess = GenerateInitialGuess<R, 3>(R(0.0), R(1.0), [](const auto& t)
				{
					return mesh_point<R, 3>{ t, 1, t };
				}, discretization);

			const auto step = 1.0 / discretization;

			trapezoidal_solver<R> solver{};
			solver.DoMeshRefinement = use_reparametrization;

			const auto solution = solver.solve(pr.get_system(), init_guess, { {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
				step, use_reparametrization, use_inversion);

			Assert::IsTrue(solver.success(), L"Failed to achieve desired precision or iteration procedure is divergent.");

			CheckQuadraticConvergence(solver, quadratic_convergence_tolerance);

			const auto init_slope_diff = auxutils::Abs(solution[0][1] - 0.25560421);
			//check initial slope
			Assert::IsTrue(init_slope_diff < 2 * step * step, L"Too big deviation from the referance initial slope value");

			const auto final_slope_diff =auxutils::Abs((*solution.rbegin())[1] - 4.266223);
			//check final slope
			Assert::IsTrue(final_slope_diff < 15 * step * step, L"Too big deviation from the referance final slope value");
		}

		TEST_METHOD(TroeschProblemTest)
		{
			perform_Troesch_test(false, false);
		}

		TEST_METHOD(TroeschProblemReparamTest)
		{
			perform_Troesch_test<number<cpp_dec_float<50>, et_off>>(true, true);
		}

		/// <summary>
		/// General method to perform test with arbitrary problem
		/// </summary>
		template <class R>
		void perform_test(const simple_bvp<R, 2>& pr, const bool first_func_bc, const R tolerance_factor,
			const bool use_reparametrization = false, const bool use_inversion = false,
			const int discretization = 100, const R t0 = R(0), const R t1 = R(3),
			const R& quadratic_convergence_tolerance = R(0.5), const R& perturbationMagnitude = R(0.2))
		{
			const auto reference_solution = pr.get_solution();

			Assert::IsTrue(std::all_of(reference_solution.begin(), reference_solution.end(), [](const auto f) { return f != nullptr; }), L"Exact solution is unknown");

			const auto reference = GenerateInitialGuess<R, 3>(t0, t1, [t0, t1, &reference_solution, first_func_bc](const auto& t)
				{
					return mesh_point<R, 3>{reference_solution[0](t),	reference_solution[1](t), t};
				}, discretization);

			const auto init_guess = perturbate(reference, perturbationMagnitude, first_func_bc);

			const auto init_guess_average_deviation = get_average_deviations(reference_solution, init_guess);

			Assert::IsTrue(std::all_of(init_guess_average_deviation.begin(), init_guess_average_deviation.end(),
				[perturbationMagnitude](const auto x) { return x > perturbationMagnitude/R(4); }),
				L"Unexpectedly low deviation for the initial guess");

			const auto step = (t1 - t0) / discretization;

			trapezoidal_solver<R> solver{};
			solver.DoMeshRefinement = use_reparametrization;

			const auto solution = solver.solve(pr.get_system(), init_guess, { {first_func_bc, !first_func_bc},{ first_func_bc, !first_func_bc} },
				100*std::numeric_limits<R>::epsilon(), step, use_reparametrization, use_inversion);

			Assert::IsTrue(solver.success(), L"Failed to achieve decired precision");

			const auto correction_magnitudes = solver.get_correcion_magnitudes();
			Assert::IsTrue(use_reparametrization || std::all_of(correction_magnitudes.begin(), correction_magnitudes.end(),
				[](const auto m) { return !m.RefinementApplied; }), L"Mesh refinement was applied unexpectedly");

			CheckQuadraticConvergence(solver, quadratic_convergence_tolerance);

			const auto diff = auxutils::Abs((*solution.rbegin())[1] - (*reference.rbegin())[1]);

			if (first_func_bc)
				Assert::IsTrue(auxutils::Abs(solution[0][0] - reference[0][0]) < 4e4 * std::numeric_limits<R>::epsilon() &&
					auxutils::Abs((*solution.rbegin())[0] - (*reference.rbegin())[0]) < 4e4 * std::numeric_limits<R>::epsilon(),
					L"The first function violates boundary conditions.");
			else
				Assert::IsTrue(auxutils::Abs(solution[0][1] - reference[0][1]) < 4e4 * std::numeric_limits<R>::epsilon() &&
					auxutils::Abs((*solution.rbegin())[1] - (*reference.rbegin())[1]) < 4e4 * std::numeric_limits<R>::epsilon(),
					L"The second function violates boundary conditions.");

			const auto max_deviations = get_max_deviations(reference_solution, solution);

			for (int func_id = 0; func_id < max_deviations.size(); func_id++)
				Logger::WriteMessage((std::string("Function ") + auxutils::ToString(func_id) + " max. deviation : " + auxutils::ToString(max_deviations[func_id]) + "\n").c_str());

			const auto max_allowed_deviation = step * step * tolerance_factor;

			Assert::IsTrue(std::all_of(max_deviations.begin(), max_deviations.end(), [max_allowed_deviation](const auto x) { return x < max_allowed_deviation; }), L"Unexpectedly high deviation");
		}

		TEST_METHOD(BVP_1_FirstFuncBCTest)
		{
			perform_test(bvp_sys_factory<double>::BVP_1(), true, 4.0, false, false, 100, 0.0, 3.0, 0.05, 1.0);
		}

		TEST_METHOD(BVP_1_SecondFuncBCTest)
		{
			perform_test(bvp_sys_factory<double>::BVP_1(), false, 1.0, false, false, 100, 0.0, 3.0, 0.05, 1.0);
		}

		TEST_METHOD(BVP_2_ReparamAndInversionFirstFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.55, true, true, 100, -3, 3, 0.64, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamAndInversionSecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), false, 0.5, true, true, 100, -3, 3, 0.48, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamOnlyFirtsFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.5, true, false, 100, -3, 3, 0.03, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamOnlySecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), false, 0.5, true, false, 100, -3, 3, 1.4, 0.2);
		}

		TEST_METHOD(BVP_2_InversionOnlyFirstFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 3.0, false, true, 300, -3, 3, 0.09, 0.2);
		}

		TEST_METHOD(BVP_2_InversionOnlySecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), false, 3.1, false, true, 300, -3, 3, 0.08, 0.2);
		}

		TEST_METHOD(BVP_2_PureTrapezoidalFirstFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.25, false, false, 100, -3, 3, 0.35, 0.2);
		}

		TEST_METHOD(BVP_2_PureTrapezoidalSecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.25, false, false, 100, -3, 3, 0.35, 0.2);
		}

		template <class R>
		void perform_BVP_T30_1_test(const R lambda, const int discretization = 100,
			const bool use_reparametrization = false, const bool use_inversion = false,
			const R& derivative_threshold = R(1), const R& quadratic_convergence_tolerance = 1.75)
		{
			const auto pr = bvp_sys_factory<R>::BVP_T30_1(lambda);

			const auto t0 = pr.get_bc().left.t;
			const auto u0 = pr.get_bc().left.u;

			const auto t1 = pr.get_bc().right.t;
			const auto u1 = pr.get_bc().right.u;
			const auto tangent = (u1 - u0) / (t1 - t0);

			auto init_guess = GenerateInitialGuess<R, 3>(t0, t1, [u0, u1, t0, t1, tangent](const auto& t)
				{
					return mesh_point<R, 3>{ u0 + (t - t0) * (u1 - u0) / (t1 - t0), tangent, t };
				}, discretization);

			const auto lambdas = std::vector<R>({R(2), R(5),  R(10), R(20), R(40), R(60), R(80), R(100), R(140),
				R(180), R(200), R(220), R(260), R(300), R(350), R(400), R(450), R(500), R(550), R(600), R(700), R(800), R(900), R(1000)
				, R(1100) , R(1200) , R(1300) , R(1400) , R(1500) , R(1600) /*, R(1700)*/ , R(1800) , R(1900), R(2000)
				, R(2100) , R(2200) , R(2300) , R(2400) , R(2500) , R(2600) , R(2700) , R(2800) , R(2900), R(3000) 
				, R(3100) , R(3200) , R(3300) , R(3400) , R(3500) , R(3600) , R(3700) , R(3800) , R(3900), R(4000)
				, R(4100) , R(4200) , R(4300) , R(4400) , R(4500), R(4600) , R(4700) , R(4800) , R(4900), R(5000)
				, R(5100) , R(5200) , R(5300) , R(5400) , R(5500) , R(5600) , R(5700) , R(5800) , R(5900), R(6000) 
				, R(6100) , R(6200) , R(6300) , R(6400) , R(6500) , R(6600) , R(6700) , R(6800) , R(6900), R(7000)
				, R(7100) , R(7200) , R(7300) , R(7400) /*, R(7500)*/ , R(7600) , R(7700) , R(7800) , R(7900), R(8000)
				, R(8100) , R(8200) , R(8300) , R(8400) , R(8500) , R(8600) , R(8700) , R(8800) , R(8900), R(9000)
				, R(9100) , R(9200) , R(9300) , R(9400) , R(9500) , R(9600) , R(9700) , R(9800) , R(9900), R(10000)
				, R(10500) , R(11000)
				, R(11500) , R(12000)
				, R(12500) , R(13000)
				, R(13500) , R(14000)
				, R(14500) , R(15000)
				, R(15500) , R(16000)
				, R(16500) , R(17000)
				, R(17500) , R(18000)
				, R(18500) , R(19000)
				, R(19500) , R(20000)
				});

			const auto second_deriv_refinement_threshold = [](const auto lambda)
			{
				return R(20);
			};

			const auto min_h_threshold = [](const auto lambda)
			{
				if (lambda < 220)
					return R(5e-4);
				
				if (lambda < 4500)
					return R(5e-5);

				if (lambda < 11500)
					return R(2e-5);

				return R(1e-5);
			};

			const auto step = [=](const auto lambda)
			{
				if (lambda < 7500)
					return R(0.01);

				if (lambda < 10500)
				return R(0.007);

				if (lambda < 15000)
					return R(0.004);

				if (lambda < 19000)
					return R(0.003);

				return R(0.0024);
			};

			trapezoidal_solver<R> solver{};
			solver.DoMeshRefinement = use_reparametrization;

			for (const auto l : lambdas)
			{
				Logger::WriteMessage((std::string("Lambda = ") + std::to_string(l) + "\n").c_str());

				init_guess = solver.solve(bvp_sys_factory<R>::BVP_T30_1(l).get_system(), init_guess,
					{ {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
					step(l), use_reparametrization, use_inversion, derivative_threshold,
					second_deriv_refinement_threshold(l), min_h_threshold(l),
					std::nullopt, //"swap" map
					std::array<bool, 2>{false, true} //"flip" map
				);

				Assert::IsTrue(solver.success(), (std::wstring(L"Failed to achieve desired precision or iteration procedure is divergent. L = ") +
					std::to_wstring(l)).c_str());

				CheckQuadraticConvergence(solver, R(100), true);

				const auto diff_left = auxutils::Abs((*init_guess.begin())[0] - u0);
				Logger::WriteMessage((std::string("diff_left =") + auxutils::ToString(diff_left) + std::string("\n")).c_str());
				Assert::IsTrue(auxutils::Abs((*init_guess.begin())[0] - u0) < std::numeric_limits<R>::epsilon(),
					L"Too big deviation form the left boundary condition");

				const auto diff_right = auxutils::Abs((*init_guess.rbegin())[0] - u1);
				Logger::WriteMessage((std::string("diff_right =") + auxutils::ToString(diff_right) + std::string("\n")).c_str());
				Assert::IsTrue(auxutils::Abs((*init_guess.rbegin())[0] - u1) < std::numeric_limits<R>::epsilon(),
					L"Too big deviation form the right boundary condition");

				const auto start_slope = auxutils::Abs((*init_guess.begin())[1]);
				Logger::WriteMessage((std::string("start_slope =") + auxutils::ToString(start_slope) + std::string("\n")).c_str());
				Assert::IsTrue(l < 110 || start_slope < 300 * std::numeric_limits<R>::epsilon(),
					L"Unexpected slope at the left boundary point");

				const auto end_slope = auxutils::Abs((*init_guess.rbegin())[1]);
				Logger::WriteMessage((std::string("end_slope =") + auxutils::ToString(end_slope) + std::string("\n")).c_str());
				Assert::IsTrue(l < 110 || end_slope < 25 *std::numeric_limits<R>::epsilon(),
					L"Unexpected slope at the right boundary point");

				(*init_guess.rbegin())[0] = u1;//Ensure that the initial guess has "perfect" values of boundary conditions
			}
		}

		TEST_METHOD(BVP_T30_1_Test)
		{
			typedef double R;

			perform_BVP_T30_1_test<R>(R(100.0), 100, true, true, R(1.0));
		}

		TEST_METHOD(BVP_T30_full_reparametrization)
		{
			using R = double;

			const auto l_0 = 100;
			auto init_guess = auxutils::ReadFromBinaryFile<mesh_point<R, 3>>((std::string("TestData\\BVP_T30_solution_") + std::to_string((int)l_0) + std::string(".dat")).c_str());

			const auto lambdas = std::vector<R>({R(140),
				R(180), R(200), R(220), R(260), R(300), R(350), R(400), R(450), R(500), R(550), R(600), R(700), R(800), R(900), R(1000)
				, R(1100) , R(1200) , R(1300) , R(1400) , R(1500) , R(1600) , R(1700) , R(1800) , R(1900), R(2000)
				, R(2100) , R(2200) , R(2300) , R(2400) , R(2500) , R(2600) , R(2700) , R(2800) , R(2900), R(3000)
				, R(3100) , R(3200) , R(3300) , R(3400) , R(3500) , R(3600) , R(3700) , R(3800) , R(3900), R(4000)
				, R(4100) , R(4200) , R(4300) , R(4400) , R(4500), R(4600) , R(4700) , R(4800) , R(4900), R(5000)
				, R(5100) , R(5200) , R(5300) , R(5400) , R(5500) , R(5600) , R(5700) , R(5800) , R(5900), R(6000)
				, R(6100) , R(6200) , R(6300) , R(6400) , R(6500) , R(6600) , R(6700) , R(6800) , R(6900), R(7000)
				, R(7100) , R(7200) , R(7300) , R(7400) , R(7500) , R(7600) , R(7700) , R(7800) , R(7900), R(8000)
				, R(8100) , R(8200) , R(8300) , R(8400) , R(8500) , R(8600) , R(8700) , R(8800) , R(8900), R(9000)
				, R(9100) , R(9200) , R(9300) , R(9400) , R(9500) , R(9600) , R(9700) , R(9800) , R(9900), R(10000)
				, R(10500) , R(11000)
				, R(11500) , R(12000)
				, R(12500) , R(13000)
				, R(13500) , R(14000)
				, R(14500) , R(15000)
				, R(15500) , R(16000)
				, R(16500) , R(17000)
				, R(17500) , R(18000)
				, R(18500) , R(19000)
				, R(19500) , R(20000)
				});

			const auto use_reparametrization = true;
			const auto use_inversion = false;
			const auto derivative_threshold = 1.0;

			const auto step = [=](const auto lambda)
			{
				if (lambda < 600)
					return R(0.1);

				if (lambda < 5700)
					return R(0.03);

				return R(0.01);

			};

			const auto second_deriv_refinement_threshold = [](const auto lambda)
			{
				if (lambda < 3600)
					return R(1);

				if (lambda < 19000)
					return R(0.1);

				return R(0.03);
			};

			trapezoidal_solver<R, ts_experimental> solver{};
			solver.DoPreRefinement = true;
			solver.OptimizeStepSize = true;
			solver.RefinementRemoveFactor = R(0.01);
			solver.DoMeshRefinement = true;

			for (const auto l : lambdas)
			{
				if ((int)l == 5900)
					int i = 0;

				Logger::WriteMessage((std::string("Lambda = ") + std::to_string(l) + "\n").c_str());

				init_guess = solver.solve(bvp_sys_factory<R>::BVP_T30_1(l).get_system(), init_guess,
					{ {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
					step(l), use_reparametrization, use_inversion, derivative_threshold,
					second_deriv_refinement_threshold(l), R(1e-5),
					std::nullopt, //"swap" map
					std::nullopt //"flip" map
				);

				Assert::IsTrue(solver.success(), (std::wstring(L"Failed to achieve desired precision or iteration procedure is divergent. L = ") +
					std::to_wstring(l)).c_str());

				CheckQuadraticConvergence(solver, R(100), true);
			}
		}

		template <class R>
		void perform_BVP_T20_test(const R lambda, const int discretization = 100,
			const bool use_reparametrization = false, const bool use_inversion = false,
			const R& derivative_threshold = R(1), const R& quadratic_convergence_tolerance = 1.75)
		{
			const auto alpha = R(0.5);

			const auto lambdas = std::vector<R>({ R(0.5), R(1), R(2), R(5),  R(10), R(20), R(40), R(60), R(80), R(100), R(140) });

			for (const auto l : lambdas)
			{

				const auto pr = bvp_sys_factory<R>::BVP_T20(l, alpha);

				const auto t0 = pr.get_bc().left.t;
				const auto u0 = pr.get_bc().left.u;

				const auto t1 = pr.get_bc().right.t;
				const auto u1 = pr.get_bc().right.u;
				const auto tangent0 = (R(1) - u0) / (alpha - t0);
				const auto tangent1 = (u1 - R(1)) / (t1 - alpha);

				auto init_guess = GenerateInitialGuess<R, 3>(t0, t1, [u0, u1, t0, t1, tangent0, tangent1, alpha](const auto& t)
					{
						if (t < alpha)
							return mesh_point<R, 3>{ u0 + (t - t0) * tangent0, tangent0, t };
						if (t > alpha)
							return mesh_point<R, 3>{ R(1) + (t - alpha) * tangent1, tangent1, t };

						return mesh_point<R, 3>{ R(1), R(0), t };

					}, discretization);

				auxutils::SaveToTextFile(init_guess, (std::string("D:\\Development\\SandBox\\init_guess_iter_") + std::to_string((int)l) + std::string(".txt")).c_str());


				trapezoidal_solver<R/*, ts_experimental*/> solver{};

				const auto pr_work = bvp_sys_factory<R>::BVP_T20(l, alpha);

				init_guess[0][0] = pr_work.get_bc().left.u;
				(*init_guess.rbegin())[0] = pr_work.get_bc().right.u;

				init_guess = solver.solve(pr_work.get_system(), init_guess,
					{ {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
					R(0.1), use_reparametrization, use_inversion, derivative_threshold,
					R(20), R(1e-4),
					std::nullopt, //"swap" map
					std::nullopt //"flip" map
				);

				auxutils::SaveToTextFile(init_guess, (std::string("D:\\Development\\SandBox\\solution_") + std::to_string((int)l) + std::string(".txt")).c_str());
			}
		}

		///Below you can fing a functionality that was used to generate results of numerical experiments in arXiv:2106.09928

		/// <summary>
		/// Describes 4 modes the trapezoidal scheme can be used for the case of Troesch's problem
		/// </summary>
		enum TrapezoidalMode
		{
			/// <summary>
			/// Regular trapezoidal scheme on an uniform mesh
			/// </summary>
			Regular_uniform,
			/// <summary>
			/// Regular trapezoidal scheme on a non-uniform mesh
			/// </summary>
			Regular_nonuniform,
			/// <summary>
			/// I-SP_1FP_2 transformation strategy (see arXiv:2106.09928)
			/// </summary>
			I_SP1FP2,
			/// <summary>
			/// I-SP_2-SP_1FP_2 transformation strategy (see arXiv:2106.09928)
			/// </summary>
			I_SP2_SP1FP2,
		};

		/// <summary>
		/// Converts trapezoidal mode into a string
		/// </summary>
		std::string ModeToString(const TrapezoidalMode& mode)
		{
			switch (mode)
			{
			case Regular_uniform: return "Regular_uniform";
			case Regular_nonuniform: return "Regular_nonuniform";
			case I_SP1FP2: return "I_SP1FP2";
			case I_SP2_SP1FP2: return "I_SP2_SP1FP2";
			default: throw std::exception("Unknown mode");
			}
		}

		/// <summary>
		/// This method was used to generate numerical data published in arXiv:2106.09928
		/// </summary>
		/// <param name="mode"> Mode of the Trapezoidal solved we want to use </param>
		/// <param name="M"> "M" parameter of the mesh refinement procedure (see arXiv:2106.09928), ignored if `mode` == `Regular_uniform`</param>
		/// <param name="h_min">h_{min} parameter of the mesh refinement procedure (see arXiv:2106.09928), ignored if `mode` == `Regular_uniform`</param>
		/// <param name="h_max">h_{max} parameter of the mesh refinement procedure (see arXiv:2106.09928), ignored if `mode` == `Regular_uniform`</param>
		template <class R>
		void arXiv_2106_09928_experiment_data_generator(const TrapezoidalMode& mode, const R& M, const R& h_min, const R& h_max, const std::string folder)
		{
			TroeschReferences reference; //access to the reference data for the Troesch's problem

			//Check the machine epsilon for the current floating point type
			Logger::WriteMessage((auxutils::ToString(auxutils::estimate_epsilon<R>()) + "\n").c_str());

			const bool use_reparametrization = ((mode == TrapezoidalMode::I_SP1FP2) || (mode == TrapezoidalMode::I_SP2_SP1FP2));
			const bool use_inversion = use_reparametrization;
			const int discretization = auxutils::ToDouble(R(1) / h_max);

			auto init_guess = GenerateInitialGuess<R, 3>(R(0.0), R(1.0), [](const auto& t)
				{
					return mesh_point<R, 3>{ t, 1, t };
				}, discretization);

			trapezoidal_solver<R> solver{};
			solver.DoMeshRefinement = mode != TrapezoidalMode::Regular_uniform;
			solver.set_max_iterations(30);
			solver.RefinementRemoveFactor = R(0.01);

			trapezoidal_solver<R, ts_experimental_troesch> solver_I_SP2_SP1FP2{};
			solver_I_SP2_SP1FP2.DoMeshRefinement = false;
			solver_I_SP2_SP1FP2.set_max_iterations(30);
			solver_I_SP2_SP1FP2.RefinementRemoveFactor = R(0.0001);

			auto file_name_base = std::string("Report_") + ModeToString(mode);

			if (mode != TrapezoidalMode::Regular_uniform)
				file_name_base += std::string("_M_") + auxutils::ToString(auxutils::ToDouble(M)) +
				std::string("_h_min_") + auxutils::ToString(auxutils::ToDouble(h_min)) + std::string("_h_max_") + auxutils::ToString(auxutils::ToDouble(h_max));

			const auto path_base = folder + file_name_base;

			table_logger<10> logger(path_base + std::string(".txt"), { "Lambda", "Tangent left", "Tangent right", "Value right", "Final correction", "Points count",
				"Left rel. diff (Maple)", "Right rel. diff (Maple)", "Left rel. diff (Vazquez)", "Right rel. diff (Vazquez)" });

			for (int l = 2; l < 501; l += 1)
			{
				const auto sys = bvp_sys_factory<R>::Troesch(l).get_system();

				init_guess = solver.solve(sys, init_guess,
					{ {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
					h_max, use_reparametrization, use_inversion, R(1), M, h_min);

				const auto final_correction = solver.get_correcion_magnitudes()[solver.get_correcion_magnitudes().size() - 1].CorrectionMagnitude;

				const auto left_tangent = init_guess[0][1];
				const auto right_tangent = init_guess[init_guess.size() - 1][1];
				const auto value_right = init_guess[init_guess.size() - 1][0];

				const auto rel_diff_maple = reference.get_rel_diff(l, auxutils::ToDouble(left_tangent), auxutils::ToDouble(right_tangent), true /*Maple*/);
				const auto rel_diff_Vazquez = reference.get_rel_diff(l, auxutils::ToDouble(left_tangent), auxutils::ToDouble(right_tangent), false /*Maple*/);

				logger.write_line<R>({ R(l), left_tangent, right_tangent, value_right, final_correction, R(init_guess.size()),
					rel_diff_maple.LeftPointTangent, rel_diff_maple.RightPointTangent, rel_diff_Vazquez.LeftPointTangent, rel_diff_Vazquez.RightPointTangent });

				Logger::WriteMessage((std::string("Lambda = ") + std::to_string(l) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Tangent left = ") + auxutils::ToString(init_guess[0][1]) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Tangent right = ") + auxutils::ToString(init_guess[init_guess.size() - 1][1]) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Value right = ") + auxutils::ToString(init_guess[init_guess.size() - 1][0]) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Argument right = ") + auxutils::ToString(init_guess[init_guess.size() - 1][2]) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Points count = ") + std::to_string(init_guess.size()) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("final_correction = ") + auxutils::ToString(final_correction) + std::string("\n")).c_str());
				Logger::WriteMessage((std::string("Total iterations = ") + auxutils::ToString(solver.get_correcion_magnitudes().size()) + std::string("\n")).c_str());
				Logger::WriteMessage("\n");

				if (l > 5 && (mode == TrapezoidalMode::I_SP2_SP1FP2))
				{
					const auto init_guess_I_SP2_SP1FP2 = solver_I_SP2_SP1FP2.solve(sys, init_guess,
						{ {true, false},{ true, false} }, 1000 * std::numeric_limits<R>::epsilon(),
						h_max, use_reparametrization, use_inversion, R(1), M, h_min);

					const auto final_correction = solver_I_SP2_SP1FP2.get_correcion_magnitudes()[solver_I_SP2_SP1FP2.get_correcion_magnitudes().size() - 1].CorrectionMagnitude;

					const auto left_tangent = init_guess_I_SP2_SP1FP2[0][1];
					const auto right_tangent = init_guess_I_SP2_SP1FP2[init_guess_I_SP2_SP1FP2.size() - 1][1];
					const auto value_right = init_guess_I_SP2_SP1FP2[init_guess_I_SP2_SP1FP2.size() - 1][0];

					const auto rel_diff_maple = reference.get_rel_diff(l, auxutils::ToDouble(left_tangent), auxutils::ToDouble(right_tangent), true /*Maple*/);
					const auto rel_diff_Vazquez = reference.get_rel_diff(l, auxutils::ToDouble(left_tangent), auxutils::ToDouble(right_tangent), false /*Maple*/);

					logger.write_line<R>({ R(l), left_tangent, right_tangent, value_right, final_correction, R(init_guess_I_SP2_SP1FP2.size()),
						rel_diff_maple.LeftPointTangent, rel_diff_maple.RightPointTangent, rel_diff_Vazquez.LeftPointTangent, rel_diff_Vazquez.RightPointTangent });
				}

				Assert::IsTrue(final_correction < R(1));
				init_guess[init_guess.size() - 1][0] = bvp_sys_factory<R>::Troesch(l).get_bc().right.u;
				Assert::IsTrue(solver.success(), L"Failed to achieve desired precision or iteration procedure is divergent.");
			}
		}

		TEST_METHOD(arXiv_2106_09928_Troesch_regular_uniform_experiments)
		{
			return; // remove to enable
			arXiv_2106_09928_experiment_data_generator<double>(TrapezoidalMode::Regular_uniform, 0.1, 1e-3, 1e-2, "D:\\Development\\SandBox\\Test_0\\");
		}
		TEST_METHOD(arXiv_2106_09928_Troesch_regular_nonuniform_experiments)
		{
			return; // remove to enable
			arXiv_2106_09928_experiment_data_generator<double>(TrapezoidalMode::Regular_nonuniform, 0.1, 1e-3, 1e-2, "D:\\Development\\SandBox\\Test_1\\");
		}

		TEST_METHOD(arXiv_2106_09928_Troesch_I_SP1FP2_experiments)
		{
			return; // remove to enable
			arXiv_2106_09928_experiment_data_generator<double>(TrapezoidalMode::I_SP1FP2, 0.1, 1e-3, 1e-2, "D:\\Development\\SandBox\\Test_2\\");
		}

		TEST_METHOD(arXiv_2106_09928_Troesch_I_SP2_SP1FP2_experiments)
		{
			return; // remove to enable
			arXiv_2106_09928_experiment_data_generator<double>(TrapezoidalMode::I_SP2_SP1FP2, 0.1, 1e-3, 1e-2, "D:\\Development\\SandBox\\Test_3\\");
		}

		TEST_METHOD(arXiv_2106_09928_Troesch_I_SP2_SP1FP2_multiprecision_experiments)
		{
			return; // remove to enable
			arXiv_2106_09928_experiment_data_generator<number<cpp_dec_float<110>, et_off>>(TrapezoidalMode::I_SP2_SP1FP2, 0.1, 1e-3, 1e-2, "D:\\Development\\SandBox\\Test_4\\");
		}
	};
}