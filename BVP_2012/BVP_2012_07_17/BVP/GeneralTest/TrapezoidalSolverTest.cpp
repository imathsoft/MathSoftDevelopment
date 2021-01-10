#include "CppUnitTest.h"
#include <vector>
#include "../BVP/LinearAlgebra/Matrix.h"
#include "../BVP/Systems/bvp_sys_factory.h"
#include "../BVP/Systems/FiniteDifferenceSolver/trapezoidal_solver.h"

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
		void CheckQuadraticConvergence(const trapezoidal_solver<R>& solver)
		{
			const auto corrections = solver.get_correcion_magnitudes();

			Assert::IsTrue(corrections.size() >= 4, L"Too few iterations to detect quadratic convergence rate");

			for (int corr_id = 0; corr_id < corrections.size() - 2; corr_id++)
			{
				Assert::IsTrue(corrections[corr_id] * corrections[corr_id] > corrections[corr_id + 1], L"The convergence rate is not quadratic");
			}
		}

		TEST_METHOD(TroeschProblemTest)
		{
			const int discretization = 1000;

			const auto lambda = 3.0;
			const auto pr = bvp_sys_factory<double>::Troesch(lambda);
			const auto init_guess = GenerateInitialGuess<double, 3>(0.0, 1.0, [](const auto& t)
				{
					return mesh_point<double, 3>{ t, 1, t };
				}, discretization);

			const auto step = 1.0 / discretization;

			trapezoidal_solver<double> solver;

			const auto solution = solver.solve(pr.get_system(), init_guess, { {true, false},{ true, false} }, 1e-12);

			Assert::IsTrue(solver.success(), L"Failed to achieve decired precision");

			CheckQuadraticConvergence(solver);

			const auto init_slope_diff = std::abs(solution[0][1] - 0.25560421);
			//check initial slope
			Assert::IsTrue(init_slope_diff < 2*step * step, L"Too big deviation from the referance initial slope value");

			const auto final_slope_diff = std::abs((*solution.rbegin())[1] - 4.266223);
			//check final slope
			Assert::IsTrue(final_slope_diff < 15 * step * step, L"Too big deviation from the referance final slope value");
		}
	};
}