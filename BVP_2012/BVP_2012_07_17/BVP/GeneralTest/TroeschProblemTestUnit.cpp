
#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "..\BVP\Problems\TroeschProblem.h"

using namespace UnitTestAux;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef number<cpp_dec_float<40>, et_off> multiprec_float;

namespace GeneralTest
{
	TEST_CLASS(TroeschProblemTestUnit)
	{
	public:
		
		TEST_METHOD(TroeschProblemInputDataTest)
		{
			TroeschProblem<multiprec_float> tpf(1);

			auto N = tpf.GetNonLin();
			Assert::IsTrue(abs(N(multiprec_float("1")) - multiprec_float("1.17520119364380145688238185059560081516")) < 1e-38, 
				Message("Function mismatch"));
			Assert::IsTrue(abs(N(multiprec_float("0.01")) - multiprec_float("1.00001666675000019841297398614117380176")) < 1e-38, 
				Message("Function mismatch"));

			auto dN = tpf.GetDerivNonLin();
			Assert::IsTrue(abs(dN(multiprec_float("1")) - multiprec_float("0.36787944117144232159552377016146086744")) < 1e-38, 
				Message("Derivative mismatch"));
			Assert::IsTrue(abs(dN(multiprec_float("0.01")) - multiprec_float("0.003333366666785714506173090027449403")) < 1e-36, 
				Message("Derivative mismatch"));

			auto ddN = tpf.GetSecondDerivNonLin();
			Assert::IsTrue(abs(ddN(multiprec_float("1")) - multiprec_float("0.43944231130091681369133431027267908028")) < 1e-37, 
				Message("Derivative mismatch"));
			Assert::IsTrue(abs(ddN(multiprec_float("0.01")) - multiprec_float("0.3333433333928572971783559806512932")) < 1e-34, 
				Message("Derivative mismatch"));
		}
	};
}