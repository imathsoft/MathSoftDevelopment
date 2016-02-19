#include "CppUnitTest.h"

#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "..\BVP\Problems\OscillatingTestProblem.h"

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

	};
}