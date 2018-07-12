#include <mpreal.h>
#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "..\BVP\Problems\TroeschProblem.h"
#include "..\BVP\FunctionApproximation\PointSimple.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/FunctionApproximation/InitialCondition.h"

using namespace auxutils;

using namespace mpfr;
using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef float_50_noet numTypeMp;
typedef double numType;

namespace GeneralTest
{
	TEST_CLASS(MultipleShootingHybridTestUnit)
	{
	public:
		
		TEST_METHOD(MultipleShootingHybridTest)
		{
			TroeschProblem<numType> tp(20);
			try
			{
				PointSimple<numType> ptLeft;
				ptLeft.Argument  = 0;
				ptLeft.Value  = 0;

				PointSimple<numType> ptRight;
				ptRight.Argument  = 1;
				ptRight.Value  = 1;

				HybridMultipleShootingComponent<numType> HMSComp(tp);

				bool succeeded;
				std::vector<InitCondition<numType>> solution = HMSComp.Run(ptLeft, ptRight, 0.0001, succeeded);

			    Assert::IsTrue(succeeded, Message("Algorithm has not succeeded"));
				auto startDerivativeDiff = abs(solution[0].Derivative - 1.64877350732915e-008); 
				Assert::IsTrue(startDerivativeDiff <= 5e-20, Message("du(0) is different" + auxutils::ToString(startDerivativeDiff)));
				auto finalDerivativeDiff = abs(solution[solution.size() - 1].Derivative - 22026.4657494062);
				Assert::IsTrue(finalDerivativeDiff <= 1e-10, Message("du(1) is different" + auxutils::ToString(finalDerivativeDiff)));
			}
			catch (exception e)
			{
				Assert::IsTrue(false, Message(e.what()));

			}
			// TODO: Your test code here
		}

		TEST_METHOD(MultipleShootingHybridTestMultiPrec)
		{
			TroeschProblem<numTypeMp> tp(20);
			try
			{
				PointSimple<numTypeMp> ptLeft;
				ptLeft.Argument  = 0;
				ptLeft.Value  = 0;

				PointSimple<numTypeMp> ptRight;
				ptRight.Argument  = 1;
				ptRight.Value  = 1;

				HybridMultipleShootingComponent<numTypeMp> HMSComp(tp);

				bool succeeded;
				std::vector<InitCondition<numTypeMp>> solution = HMSComp.Run(ptLeft, ptRight, (numTypeMp)1/50, succeeded);

			    Assert::IsTrue(succeeded, Message("Algorithm has not succeeded"));

				auto startDerivativeDiff = abs(solution[0].Derivative - (numTypeMp)"1.65428160979175801124447457829e-08");
				Assert::IsTrue(startDerivativeDiff <= std::numeric_limits<numTypeMp>::epsilon(), 
					Message("du(0) is different " + auxutils::ToString(solution[0].Derivative)+ " " + auxutils::ToString(startDerivativeDiff)));

				auto finalDerivativeDiff = abs(solution[solution.size() - 1].Derivative - (numTypeMp)"22026.4657402382174165421786901");
				Assert::IsTrue( finalDerivativeDiff <= 1e-24, 
					Message("du(1) is different " + auxutils::ToString(solution[solution.size() - 1].Derivative) + " " + auxutils::ToString(finalDerivativeDiff)));
			}
			catch (exception e)
			{
				Assert::IsTrue(false, Message(e.what()));

			}
			// TODO: Your test code here
		}

	};
}