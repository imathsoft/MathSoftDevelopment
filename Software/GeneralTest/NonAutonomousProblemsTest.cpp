#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "..\BVP\Problems\NonAutonomousTroeschProblem.h"
#include "..\BVP\Problems\NonAutonomousOscillatingProblem.h"
#include "..\BVP\FunctionApproximation\PointSimple.h"
#include "..\BVP\MultipleShooting\HybridMultipleShootingComponent.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../BVP/FunctionApproximation/InitialCondition.h"

using namespace auxutils;

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef float_50_noet numTypeMp;
typedef double numType;

namespace GeneralTest
{
	TEST_CLASS(NonAutonomousProblemsTest)
	{
	public:

		TEST_METHOD(StandardNonAutonomousOscilatingProblemMultimpeShootingDouble)
		{
			 NonAutonomousOscillatingProblem<numType> problem;
   			 StandardOscillatinProblemMultipleShoothingTest<numType, NonAutonomousOscillatingProblem<numType>>(problem);
		}
	};
}