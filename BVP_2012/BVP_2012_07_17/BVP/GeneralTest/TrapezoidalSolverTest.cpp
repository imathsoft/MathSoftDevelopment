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
		/// Method to check that the convergence rate of the iteration process is quadratic
		/// </summary>
		template <class R>
		void CheckQuadraticConvergence(const trapezoidal_solver<R>& solver, const R& acceptable_approximation_error)
		{
			const auto corrections = solver.get_correcion_magnitudes();

			std::vector<R> chosen_corrections;

			for (int corr_id = corrections.size() - 2; corr_id >= 1 && !corrections[corr_id].RefinementApplied; corr_id--)
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
			const int discretization = 100, const R t0 = R(0), const R t1 = R(3), const R& quadratic_convergence_tolerance = R(0.2))
		{
			const auto reference_solution = pr.get_solution();

			Assert::IsTrue(std::all_of(reference_solution.begin(), reference_solution.end(), [](const auto f) { return f != nullptr; }), L"Exact solution is unknown");

			const auto init_guess = GenerateInitialGuess<double, 3>(t0, t1, [t0, t1, &reference_solution, first_func_bc](const auto& t)
				{
					return mesh_point<double, 3>{ 
						first_func_bc*(reference_solution[0](t0) + (t - t0) * (reference_solution[0](t1) - reference_solution[0](t0)) / (t1 - t0)) + 2.0* !first_func_bc,
						(!first_func_bc)* (reference_solution[1](t0) + (t - t0) * (reference_solution[1](t1) - reference_solution[1](t0)) / (t1 - t0)), t };
				}, discretization);

			const auto init_guess_average_deviation = get_average_deviations(reference_solution, init_guess);

			Assert::IsTrue(std::any_of(init_guess_average_deviation.begin(), init_guess_average_deviation.end(), [](const auto x) { return x > 0.6; }),
				L"Unexpectedly low deviation for the initial guess");

			const auto step = (t1 - t0) / discretization;

			trapezoidal_solver<double> solver{};

			const auto solution = solver.solve(pr.get_system(), init_guess, { {first_func_bc, !first_func_bc},{ first_func_bc, !first_func_bc} }, 1e-12, 0.01, false, false);

			Assert::IsTrue(solver.success(), L"Failed to achieve decired precision");

			CheckQuadraticConvergence(solver, quadratic_convergence_tolerance);

			const auto max_deviations = get_max_deviations(reference_solution, solution);

			for (int func_id = 0; func_id < max_deviations.size(); func_id++)
				Logger::WriteMessage((std::string("Function ") + auxutils::ToString(func_id) + " max. deviation : " + auxutils::ToString(max_deviations[func_id]) + "\n").c_str());

			const auto max_allowed_deviation = step * step * tolerance_factor;

			Assert::IsTrue(std::all_of(max_deviations.begin(), max_deviations.end(), [max_allowed_deviation](const auto x) { return x < max_allowed_deviation; }), L"Unexpectedly high deviation");
		}

		TEST_METHOD(BVP_1_FirstFuncBCTest)
		{
			perform_test(bvp_sys_factory<double>::BVP_1(), true, 1e-1);
		}

		TEST_METHOD(BVP_1_SecondFuncBCTest)
		{
			perform_test(bvp_sys_factory<double>::BVP_1(), false, 1e-1);
		}
	};
}