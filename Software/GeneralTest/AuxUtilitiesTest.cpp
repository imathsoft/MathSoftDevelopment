#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/Systems/ode_system.h"
#include <sstream>
#include <numeric>

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

		TEST_METHOD(ColelctionBinaryStreamReadWriteTest)
		{
			const int Dim = 3;

			const auto collection_size = 1000;
			std::vector<mesh_point<double, Dim>> collection_original(collection_size);

			for (auto& pt : collection_original)
				for (int id = 0; id < Dim; id++)
					pt[id] = UnitTestAux::Random();

			const auto mean = std::accumulate(collection_original.begin(), collection_original.end(), 0.0,
				[](const auto sum, const auto pt) { return sum + pt[0]; }) / collection_original.size();;

			//sanity check 
			Assert::IsTrue(std::abs(mean - 0.5) < 0.1, L"Unexpected mean");

			std::stringstream stream;

			stream << collection_original;

			std::vector<mesh_point<double, Dim>> collection_deserialized;

			stream >> collection_deserialized;

			Assert::IsTrue(collection_original.size() == collection_deserialized.size(), L"colelctions have different size");

			for (int pt_id = 0; pt_id < collection_original.size(); pt_id++)
			{
				for (int item_id = 0; item_id < Dim; item_id++)
					Assert::IsTrue(collection_original[pt_id][item_id] == collection_deserialized[pt_id][item_id], L"Items are not the same");
			}
		}
	};
}