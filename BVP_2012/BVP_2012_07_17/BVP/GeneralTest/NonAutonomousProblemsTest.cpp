#include <mpreal.h>
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

using namespace mpfr;
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

		TEST_METHOD(NonAutonomousOscilatingProblemMultimpeShootingDouble)
		{
			 numType preH = 0.1;
			 numType finalH = 0.001;
			 numType targetValue = 0.5804096620;
			 numType targetArgument = 10;
			 NonAutonomousOscillatingProblem<numType> problem;

			 std::function<bool(const InitCondition<numType>&)> checkFunc = 
				 [=](const InitCondition<numType>& ic) { return (abs(ic.Value) <= targetArgument) 
				 && (abs(ic.Argument) <= targetArgument); };

			 HybridCannon<numType> cannon(problem, preH, preH/10.0, checkFunc);

			 std::function<int(const InitCondition<numType>&)> evalFunc = 
				 [=](const InitCondition<numType>& ic) { return sgn(ic.Value - targetValue); };
			 BisectionComponent<numType> bc(cannon);
			 Assert::IsTrue(bc.DerivativeBisectionGen(0.0, targetArgument, 1.0, 100, 0.9, 1.05, evalFunc));

			 auto knots = cannon.GetKnotVectorStreight();
			 knots[knots.size() - 1].Value = targetValue;
			 knots[knots.size() - 1].Argument = targetArgument;

			 HybridMultipleShootingComponent<numType> HMSComp(problem);

			 bool succeeded;
			 std::vector<InitCondition<numType>> solution = HMSComp.Run(knots, finalH, succeeded);

			 numType maxDev = CalcDeviationFromExactSolution<numType>(solution, 
			[](const numType& u){ return exp(sin(u)); });

			 numType vaxDistanceBetweenKnots = CalcMaxSquaredDistanceBetweenNeighbourKnots(solution);
			 Assert::IsTrue(vaxDistanceBetweenKnots < 2*finalH*finalH, 
				 Message("Too big maximal distance between neighbour knots " + 
				 auxutils::ToString(vaxDistanceBetweenKnots)));
			 Assert::IsTrue(maxDev < finalH*finalH, 
				 Message("Too big deviation, maxDev= " + auxutils::ToString(maxDev)));
			 Assert::IsTrue(solution.size() > targetArgument/finalH, 
				 Message("Too few knot in the solution vector, " + 
				 auxutils::ToString(solution.size())));
		}

	};
}