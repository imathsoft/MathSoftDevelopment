#include "stdafx.h"
#include "CppUnitTest.h"

#include <vector>
#include <algorithm>
#include <boost\timer.hpp>
#include <CppUnitTestLogger.h>
#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Utils/SolutionUtils.h"
#include <string>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
namespace GeneralTest
{
TEST_CLASS(simplificatio_test)
{
public:
	template <bool Deriv>
	void Run_simplification_test()
	{
		const auto solution = auxutils::ReadFromFile<InitCondition<double>>("TestData\\solution_double_precision.txt");
		const auto length_init = solutionutils::length<Deriv>(solution);

		boost::timer t;

		double percentage_of_points_to_keep = 0.01;
		const auto solution_simplified = solutionutils::simplify_polyline<Deriv>(solution, percentage_of_points_to_keep);

		Logger::WriteMessage((std::string("Simplification time : ") + auxutils::ToString(t.elapsed())+"\n").c_str());

		const auto length_simplified = solutionutils::length<Deriv>(solution_simplified);

		auto diff = std::abs(length_init - length_simplified);

		Logger::WriteMessage((std::string("Length before simplification : ") + auxutils::ToString(length_init)+"\n").c_str());
		Logger::WriteMessage((std::string("Length after simplification : ") + auxutils::ToString(length_simplified) + "\n").c_str());
		Logger::WriteMessage((std::string("Length difference : ") + auxutils::ToString(diff)+"\n").c_str());

		Assert::IsTrue(diff < 2e-7, L"Too big difference in th elength before and after simplification");
		Assert::IsTrue(solution_simplified.size()/solution.size() <= percentage_of_points_to_keep, L"Unexpectedly many points have been left");

		//auxutils::SaveToFile(solution_simplified, "E:\\Research\\Derivative_simplified_double_precision.txt");
		//auxutils::SaveFunction<Deriv>("x", "u", "E:\\Research\\derivative_to_plot.txt", solution_simplified);
	}

	TEST_METHOD(Function_simplification_test)
	{
		Run_simplification_test<false>();
	}

	TEST_METHOD(Derivative_simplification_test)
	{
		Run_simplification_test<true>();
	}
};
}