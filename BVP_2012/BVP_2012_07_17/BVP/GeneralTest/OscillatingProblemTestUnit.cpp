#include "CppUnitTest.h"

#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "..\BVP\Problems\OscillatingTestProblem.h"
#include "../BVP/Problems/AutonomousOscillatingProblem.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef float_50_noet numTypeMp;

namespace GeneralTest
{
	TEST_CLASS(OscillatingProblemTestUnit)
	{
	public:
		
		TEST_METHOD(OscillatingProblemTestMethod)
		{
			OscillatingTestProblem<numTypeMp> problem;

			auto N = problem.GetNonLin();
			Assert::IsTrue(abs(N((numTypeMp)1) - cos((numTypeMp)1)) < std::numeric_limits<numTypeMp>::epsilon(), Message("Function mismatch"));
			Assert::IsTrue(abs(N((numTypeMp)3.14) - cos((numTypeMp)3.14)) < std::numeric_limits<numTypeMp>::epsilon(), Message("Function mismatch"));

			auto dN = problem.GetDerivNonLin();
			Assert::IsTrue(abs(dN((numTypeMp)1) + sin((numTypeMp)1)) < std::numeric_limits<numTypeMp>::epsilon(), Message("First derivative mismatch"));
			Assert::IsTrue(abs(dN((numTypeMp)3.14) + sin((numTypeMp)3.14)) < std::numeric_limits<numTypeMp>::epsilon(), Message("First derivative mismatch"));

			auto ddN = problem.GetSecondDerivNonLin();
			Assert::IsTrue(abs(ddN((numTypeMp)1) + N((numTypeMp)1)) < std::numeric_limits<numTypeMp>::epsilon(), Message("Second derivative mismatch"));
			Assert::IsTrue(abs(ddN((numTypeMp)3.14) + N((numTypeMp)3.14)) < std::numeric_limits<numTypeMp>::epsilon(), Message("SEcond derivative mismatch"));
		}

		TEST_METHOD(OscillatingProblemCanonHybridDoubleTestMethod)
		{
			 double h = 0.01;
			 OscillatingTestProblem<double> problem;
			 HybridCannon<double> thc(problem, h, 10*std::numeric_limits<double>::epsilon());

			 auto result = thc.Shoot(0, 10, 1, 1e6, 1);
			 auto knots = thc.GetKnotVectorStreight();

			 Assert::IsTrue(abs(knots[knots.size() - 1].Value + 0.50873075567966530) < std::numeric_limits<double>::epsilon(), Message("Function mismatch"));
			 Assert::IsTrue(abs(knots[knots.size() - 1].Derivative - 0.69192015740737922) < std::numeric_limits<double>::epsilon(), Message("Derivative mismatch"));
		}

		TEST_METHOD(OscillatingProblemHybridBisectionDoubleTestMethod)
		{
			 double h = 0.01;
			 OscillatingTestProblem<double> problem;

			 std::function<bool(const InitCondition<double>&)> checkFunc = 
				 [](const InitCondition<double>& ic) { return (abs(ic.Value) < 10) && (abs(ic.Argument) < 10); };

			 HybridCannon<double> cannon(problem, h, 10*std::numeric_limits<double>::epsilon(), checkFunc);

			 std::function<int(const InitCondition<double>&)> evalFunc = 
				 [](const InitCondition<double>& ic) { return sgn(ic.Value - 1.1); };
			 BisectionComponent<double> bc(cannon);
			 Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 10.0, 1.0, 100, 0.0, 1.25, evalFunc));

			 auto knots = cannon.GetKnotVectorStreight();

			 Assert::IsTrue(abs(knots[knots.size() - 1].Value - 1.1) < 10*std::numeric_limits<double>::epsilon(), Message("Function mismatch"));
			 Assert::IsTrue(abs(knots[knots.size() - 1].Derivative + 0.32755849998534936) < std::numeric_limits<double>::epsilon(), Message("Derivative mismatch"));
			 Assert::IsTrue(abs(knots[0].Derivative - 0.054683724794614211) < std::numeric_limits<double>::epsilon(), Message("Derivative mismatch"));
		}

		TEST_METHOD(OscillatingProblemMultiShootingHybridDoubleTestMethod)
		{
			 OscillatingTestProblem<double> problem;

			 auto knots = auxutils::ReadFromFile<InitCondition<double>>("TestData\\OscillatingTestProblemMSInitGuess.txt");

			 double targetValue = knots[knots.size() - 1].Value;

			 HybridMultipleShootingComponent<double> HMSComp(problem);

			 bool succeeded;
			 std::vector<InitCondition<double>> solution = HMSComp.Run(knots, 0.002, succeeded);

			 Assert::IsTrue(succeeded, Message("Algorithm has not succeeded"));
			 Assert::IsTrue(abs(solution[solution.size() - 1].Value - targetValue) < 10*std::numeric_limits<double>::epsilon(), 
				 Message("Function mismatch " + auxutils::ToString(solution[solution.size() - 1].Value - targetValue)));
			 Assert::IsTrue(abs(solution[solution.size() - 1].Derivative - 1.31696084860275) < 100*std::numeric_limits<double>::epsilon(), 
				 Message("Derivative mismatch " + auxutils::ToString(solution[solution.size() - 1].Derivative)));
			 Assert::IsTrue(abs(solution[0].Derivative - 1.2576169315833) < 100*std::numeric_limits<double>::epsilon(), 
				 Message("Derivative mismatch " + auxutils::ToString(solution[0].Derivative)));
		}

		TEST_METHOD(StandardAutonomousOscilatingProblemMultimpeShootingDouble)
		{
			 AutonomousOscillatingProblem<double> problem;
   			 StandardOscillatinProblemMultipleShoothingTest<double, 
				 AutonomousOscillatingProblem<double>>(problem);
		}

	};
}