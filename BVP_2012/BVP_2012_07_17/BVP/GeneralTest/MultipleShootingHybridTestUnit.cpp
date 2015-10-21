#include <mpreal.h>
#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "..\BVP\Problems\TroeschProblem.h"
#include "..\BVP\FunctionApproximation\PointSimple.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"

using namespace mpfr;
using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{
	TEST_CLASS(MultipleShootingHybridTestUnit)
	{
	public:
		
		TEST_METHOD(MultipleShootingHybridTest)
		{
			TroeschProblem<double> tp(30);
			try
			{
				mpfr::mpreal::set_default_prec(128);

				PointSimple<double> ptLeft;
				ptLeft.Argument  = 0;
				ptLeft.Value  = 0;

				PointSimple<double> ptRight;
				ptRight.Argument  = 1;
				ptRight.Value  = 1;

				HybridMultipleShootingComponent<double> HMSComp(tp, ptLeft, ptRight, 1e-14);
			}
			catch (exception e)
			{
				Assert::IsTrue(false, Message(e.what()));

			}
			// TODO: Your test code here
		}

	};
}