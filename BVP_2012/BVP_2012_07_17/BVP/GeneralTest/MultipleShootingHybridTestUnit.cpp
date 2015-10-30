#include <mpreal.h>
#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "..\BVP\Problems\TroeschProblem.h"
#include "..\BVP\FunctionApproximation\PointSimple.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../BVP/Utils/AuxUtils.h"

using namespace auxutils;

using namespace mpfr;
using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef float_50_noet numType;
//typedef double numType;

namespace GeneralTest
{
	TEST_CLASS(MultipleShootingHybridTestUnit)
	{
	public:
		
		TEST_METHOD(MultipleShootingHybridTest)
		{
			TroeschProblem<numType> tp(2);
			try
			{
				mpfr::mpreal::set_default_prec(128);

				PointSimple<numType> ptLeft;
				ptLeft.Argument  = 0;
				ptLeft.Value  = 0;

				PointSimple<numType> ptRight;
				ptRight.Argument  = 1;
				ptRight.Value  = 1;

				HybridMultipleShootingComponent<numType> HMSComp(tp, ptLeft, ptRight, 1e-14);
			}
			catch (exception e)
			{
				Assert::IsTrue(false, Message(e.what()));

			}
			// TODO: Your test code here
		}

	};
}