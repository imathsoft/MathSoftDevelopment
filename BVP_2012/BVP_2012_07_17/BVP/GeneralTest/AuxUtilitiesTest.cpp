#include "CppUnitTest.h"
#include "UnitTestAux.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace LinAlg;

namespace GeneralTest
{
	TEST_CLASS(AuxUtilitiesTest)
	{
		TEST_METHOD(ApproximateParametersOfQuadraticallyDecayingSequenceTest)
		{
			double M_orig = 3;
			double q_orig = 0.375;

			int element_count = 7;

			std::vector<double> sequence;

			double power_of_q = q_orig;
			for (auto val_id = 0; val_id < element_count; val_id++)
			{
				sequence.push_back(M_orig * power_of_q);
				power_of_q *= power_of_q;
			}

			double M, q, max_rel_error;

			UnitTestAux::ApproximateParametersOfQuadraticallyDecayingSequence(sequence, M, q, max_rel_error);

			Assert::IsTrue(auxutils::Abs(M - M_orig) < 35 * std::numeric_limits<double>::epsilon(), L"Too big deviation for the M parameter");
			Assert::IsTrue(auxutils::Abs(q - q_orig) < 35 * std::numeric_limits<double>::epsilon(), L"Too big deviation for the q parameter");
			Assert::IsTrue(max_rel_error >= 0.0, L"Relative error can't be negative");
			Assert::IsTrue(max_rel_error < 100 * std::numeric_limits<double>::epsilon(), L"Too big maximal relative error of the approximation");
		}
	};
}