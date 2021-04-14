#include "CppUnitTest.h"
#include <vector>
#include <numeric>
#include "../BVP/LinearAlgebra/Matrix.h"
#include "../BVP/Systems/bvp_sys_factory.h"
#include "../BVP/Systems/FiniteDifferenceSolver/trapezoidal_solver.h"
#include "../BVP/Utils/AuxUtils.h"
#include "UnitTestAux.h"

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
		template <class R>
		void CheckQuadraticConvergence(const trapezoidal_solver<R>& solver, const R& acceptable_approximation_error)
		{
			const auto corrections = solver.get_correcion_magnitudes();

			std::vector<R> chosen_corrections;

			for (int corr_id = corrections.size() - 2; corr_id >= 1 &&
				!corrections[corr_id].RefinementApplied &&
				corrections[corr_id].CorrectionMagnitude < R(0.1); corr_id--)
				chosen_corrections.insert(chosen_corrections.begin(), corrections[corr_id].CorrectionMagnitude);

			Assert::IsTrue(chosen_corrections.size() >= 3, L"Too few corrections to detect quadratic convergence rate");

			R M, q, max_rel_error;

			UnitTestAux::ApproximateParametersOfQuadraticallyDecayingSequence(chosen_corrections, M, q, max_rel_error);

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
			const simple_bvp<R, 2>& pr = bvp_sys_factory<R>::Troesch(3.0), const R& quadratic_convergence_tolerance = R(0.3))
		{
			const auto init_guess = GenerateInitialGuess<R, 3>(R(0.0), R(1.0), [](const auto& t)
				{
					return mesh_point<R, 3>{ t, 1, t };
				}, discretization);

			const auto step = 1.0 / discretization;

			trapezoidal_solver<R> solver{};

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
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.5, true, true, 100, -3, 3, 0.6, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamAndInversionSecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), false, 0.5, true, true, 100, -3, 3, 0.25, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamOnlyFirtsFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 0.5, true, false, 100, -3, 3, 0.3, 0.2);
		}

		TEST_METHOD(BVP_2_ReparamOnlySecondFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), false, 0.5, true, false, 100, -3, 3, 0.9, 0.2);
		}

		TEST_METHOD(BVP_2_InversionOnlyFirstFuncBCTest)
		{
			perform_test<double>(bvp_sys_factory<double>::BVP_2(), true, 3.0, false, true, 300, -3, 3, 0.08, 0.2);
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
	};
}