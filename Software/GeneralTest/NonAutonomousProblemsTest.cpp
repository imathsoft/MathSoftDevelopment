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
		
		//TEST_METHOD(NonAutonomousTroeschProblemDoubleMultimpeShooting)
		//{
		//	NonAutonomousTroeschProblem<numType> tpf(20);
		//	PointSimple<numType> ptLeft;
		//	ptLeft.Argument  = 0;
		//	ptLeft.Value  = 0;

		//	PointSimple<numType> ptRight;
		//	ptRight.Argument  = 1;
		//	ptRight.Value  = 1;

		//	HybridMultipleShootingComponent<numType> HMSComp(tpf);

		//	bool succeeded;
		//	std::vector<InitCondition<numType>> solution = HMSComp.Run(ptLeft, ptRight, 0.0001, succeeded, 0.1);

		//	auxutils::SaveToMapleFile(solution, "f:\\NonAutoTroeschProblem.txt", true);

		//}

		TEST_METHOD(StandardNonAutonomousOscilatingProblemMultimpeShootingDouble)
		{
			 NonAutonomousOscillatingProblem<numType> problem;
   			 StandardOscillatinProblemMultipleShoothingTest<numType, NonAutonomousOscillatingProblem<numType>>(problem);
		}
	};
}