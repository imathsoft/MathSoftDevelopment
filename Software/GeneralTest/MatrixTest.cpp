
#include "CppUnitTest.h"
#include <Matrix.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace LinAlg;

namespace GeneralTest
{
	constexpr int RowDim = 5;
	constexpr int ColDim = 7;

	constexpr double TOLL = 20*std::numeric_limits<double>::epsilon();

	TEST_CLASS(MatrixTest)
	{
		template<class R, int RD, int CD>
		static bool MatricesAreEqual(const Matrix<R, RD, CD>& a, const Matrix<R, RD, CD>& b, const R& tolerance = TOLL)
		{
			for (int i = 0; i < RD; i++)
				for (int j = 0; j < CD; j++)
					if (std::abs(a[i][j] - b[i][j]) > tolerance)
						return false;

			return true;
		}

		template<class R, int RD, int CD>
		static bool IsZeroMatrix(const Matrix<R, RD, CD>& m, const R& tolerance = TOLL)
		{
			for (int i = 0; i < RD; i++)
				for (int j = 0; j < CD; j++)
					if (std::abs(m[i][j]) > tolerance)
						return false;

			return true;
		}

		TEST_METHOD(NormOfTransposedMatrixTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto m = Matrix<double, RowDim, ColDim>::Random();

			Assert::IsTrue(std::abs(m.Norm() - m.Transpose().Norm()) < TOLL, L"Norm of transposed matrix is not equal to the notm of the initial matrix");
		}

		TEST_METHOD(ZeroMatrixNormTest)
		{
			Matrix<double, RowDim, ColDim> m;

			Assert::IsTrue(IsZeroMatrix(m, 0.0), L"Matrix is not zero iniialized");

			const auto m_norm = m.Norm();
			Assert::IsTrue(m_norm >= 0.0, L"Matrix notm is negative");
			Assert::IsTrue(m_norm < TOLL, L"Zero matrix has positive norm");
		}

		TEST_METHOD(NormHomogeneousPropertyTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto m = Matrix<double, RowDim, ColDim>::Random();
			const auto scalar = Matrix<double, RowDim, ColDim>::Random()[0][0];

			Assert::IsTrue(MatricesAreEqual(m * scalar, scalar * m, 0.0), L"Matrix-scalar multiplication is broken");

			Assert::IsTrue(std::abs((m * scalar).Norm() - scalar * m.Norm()) < TOLL, L"Matrix norm is not homogeneous with respect to scalar multiplicant");
		}

		TEST_METHOD(NormTriangleInequalityTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, RowDim, ColDim>::Random();

			Assert::IsTrue((a + b).Norm() <= a.Norm() + b.Norm(), L"Triangle inequality does not hold true");
		}

		TEST_METHOD(NormOfProductTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, ColDim, RowDim>::Random();

			Assert::IsTrue((a * b).Norm() <= a.Norm() * b.Norm(), L"Triangle inequality does not hold true");
		}

		TEST_METHOD(ProductAssociativityTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, ColDim, RowDim>::Random();
			const auto c = Matrix<double, RowDim, ColDim>::Random();

			Assert::IsTrue(MatricesAreEqual(a * (b * c), (a * b) * c), L"Product associativity does not hold true");
		}

		TEST_METHOD(ProductDistributivityTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, RowDim, ColDim>::Random();
			const auto c = Matrix<double, ColDim, RowDim>::Random();
			const auto d = Matrix<double, ColDim, RowDim>::Random();

			Assert::IsTrue(MatricesAreEqual((a + b) * c, a * c + b * c), L"Product distributivity does not hold true");
			Assert::IsTrue(MatricesAreEqual(a * (c + d), a * c + a * d), L"Product distributivity does not hold true");
		}

		TEST_METHOD(ProductOfTransposedMatricesTest)
		{
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, ColDim, RowDim>::Random();

			Assert::IsTrue(MatricesAreEqual((a * b).Transpose(), b.Transpose() * a.Transpose()), L"Transposition of a product is not equal to the product of transposed matrices");
		}
		TEST_METHOD(NormNonNegativeTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto m = Matrix<double, RowDim, ColDim>::Random();

			Assert::IsTrue(m.Norm() >= 0, L"Norm can't be negative");
		}

		TEST_METHOD(AdditionCommutativityTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto a = Matrix<double, RowDim, ColDim>::Random();
			const auto b = Matrix<double, RowDim, ColDim>::Random();

			Assert::IsFalse(MatricesAreEqual(a, b), L"Matrices A and B coincide");

			const auto a_plus_b = a + b;
			const auto b_plus_a = b + a;

			Assert::IsTrue(MatricesAreEqual(a + b, b + a), L"Addition commutativity does not hold true");
		}

		TEST_METHOD(AdditionOfZeroMatrixTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto m = Matrix<double, RowDim, ColDim>::Random();
			Matrix<double, RowDim, ColDim> z = Matrix<double, RowDim, ColDim>{};

			Assert::IsTrue(IsZeroMatrix(z), L"Matrix is nonzero");

			Assert::IsTrue(MatricesAreEqual(m, m + z), L"Matrices are not equal");
			Assert::IsTrue(MatricesAreEqual(m, z + m), L"Matrices are not equal");
		}

		TEST_METHOD(IdentityMatrixProductTest)
		{
			Matrix<double, RowDim, ColDim>::Randomize();
			const auto m = Matrix<double, RowDim, RowDim>::Random();
			const auto i = Matrix<double, RowDim, RowDim>::Identity();

			Assert::IsTrue(MatricesAreEqual(m * i, m), L"Matrices are not equal");
			Assert::IsTrue(MatricesAreEqual(i * m, m), L"Matrices are not equal");
		}

		/// <summary>
		/// Generic method to run "standard" test of the matrix inversion functionality
		/// </summary>
		template<int Dim>
		void RunStandardInverseMatrixTest()
		{
			const auto m = Matrix<double, Dim, Dim>::Random() + Matrix<double, Dim, Dim>::Identity();

			const auto  m_inverse = m.Inverse();

			Assert::IsTrue(MatricesAreEqual(Matrix<double, Dim, Dim>::Identity(), m * m_inverse), L"Matrices are not equal");
			Assert::IsTrue(MatricesAreEqual(Matrix<double, Dim, Dim>::Identity(), m_inverse * m), L"Matrices are not equal");
		}

		TEST_METHOD(InverseMatrix1x1Test)
		{
			RunStandardInverseMatrixTest<1>();
		}

		TEST_METHOD(InverseMatrix2x2Test)
		{
			RunStandardInverseMatrixTest<2>();
		}

		TEST_METHOD(Inverse_general_test)
		{
			auto m = Matrix<double, ColDim, ColDim>::Random();

			//make sure that the random matrix is nonsingular
			//by adding ones. On purpose we do that not to the main diagonal,
			//to exercise the row swapping inside the generic inversion implementation
			for (int row_id = 0; row_id < ColDim; row_id++)
			{
				const auto col_id = ColDim - 1 - row_id;
				m[row_id][col_id] += 1.0;
			}

			const auto  inv = m.Inverse();
			const auto check_diff = (inv * m - Matrix<double, ColDim, ColDim>::Identity()).Norm();

			Logger::WriteMessage((std::string("check_diff = ") + auxutils::ToString(check_diff)).c_str());

			Assert::IsTrue(check_diff < 50 * std::numeric_limits<double>::epsilon(), L"Too high deviation");
		}

		TEST_METHOD(UnaryMinusOperatorTest)
		{
			const auto m = Matrix<double, RowDim, ColDim>::Random();
			const auto zero = Matrix<double, RowDim, ColDim>{};

			Assert::IsFalse(IsZeroMatrix(m), L"The matrix is expected to be non-zero");
			Assert::IsTrue(IsZeroMatrix(zero), L"The matrix is expected to be zero");

			const auto minus_m = -m;

			Assert::IsTrue(MatricesAreEqual(m + minus_m, zero), L"The result is expected to be zero");
		}
	};
}