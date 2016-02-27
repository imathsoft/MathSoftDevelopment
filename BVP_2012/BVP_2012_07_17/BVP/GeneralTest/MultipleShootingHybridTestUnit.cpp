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

				std::vector<InitCondition<numType>> solution = HMSComp.Run(ptLeft, ptRight, 0.0001);

				Assert::IsTrue(abs(solution[0].Derivative - 1.64877364654916e-008) <= 1e-21, 
					Message("du(0) is different" + auxutils::ToString(solution[0].Derivative)));
				Assert::IsTrue(abs(solution[solution.size() - 1].Derivative - 22026.4657494062) <= 1e-10, 
					Message("du(1) is different" + auxutils::ToString(solution[solution.size() - 1].Derivative)));
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

				std::vector<InitCondition<numTypeMp>> solution = HMSComp.Run(ptLeft, ptRight, (numTypeMp)1/50);

				Assert::IsTrue(abs(solution[0].Derivative - 1.6570519017559527928e-08)<= 1e-25, 
					Message("du(0) is different " + auxutils::ToString(solution[0].Derivative)+ " " +
					auxutils::ToString(abs(solution[0].Derivative - 1.6570519017559527928e-08))));

				Assert::IsTrue(abs(solution[solution.size() - 1].Derivative - 22026.465708960743113) <= 1e-12, 
					Message("du(1) is different " + auxutils::ToString(solution[solution.size() - 1].Derivative) + " " +
					auxutils::ToString(abs(solution[solution.size() - 1].Derivative - 22026.465708960743113))));
			}
			catch (exception e)
			{
				Assert::IsTrue(false, Message(e.what()));

			}
			// TODO: Your test code here
		}

	};
}