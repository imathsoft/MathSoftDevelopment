#include "CppUnitTest.h"
#include "../BVP/Utils/AuxUtils.h"
#include "../BVP/FunctionApproximation/DerivativeEvaluator/Dual.h"
#include <functional>
#include <random>
#include <ctime>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{
	constexpr int RND_MIN = -5;
	constexpr int RND_MAX = 5;

	TEST_CLASS(DualNumberTest)
	{
		/// <summary>
		/// Returns a "random" double
		/// </summary>
		double get_random_double()
		{
			return RND_MIN + (double)(rand()) / ((double)(RAND_MAX / (RND_MAX - RND_MIN)));
		}

		/// <summary>
		/// Function performing arithmetic operations
		/// </summary>
		template <class R>
		static R arithm_func(const R& x, const R& y, const R& z)
		{
			return (2 * x * x * x * 4 - 3 * x * x * y + 4 * x * y * z - 5 * y * y * z + 6 * y * z * z - 7 * y * y * y + 8 * z * z * z) /
				(2 * x * x + 5 * y * y + 8 * z * z + 1) + x / 3 - 1 / (x * x + 1) + 2;
		}

		/// <summary>
		/// Function performing some trigonometric operations
		/// </summary>
		template <class R>
		static R trigonometric_func(const R& x, const R& y, const R& z)
		{
			return Sin(x) + Cos(y) - Sin(z) + Sin(x) * Cos(y) * Cos(z) + (Sin(x) + Cos(z)) / (2 + Cos(x))
				+ Sin(x * x + y * y + z) - Cos(x * x - y * y - z) + Sin(Cos(x + y + z));
		}

		/// <summary>
		/// Function performing some operations with exponent function
		/// </summary>
		template <class R>
		static R exponent_func(const R& x, const R& y, const R& z)
		{
			return (Exp(1.3 * x) + Exp(y / 2) + Exp(z * 1.5) + 1) / (Exp(x) + Exp(y) + Exp(z)) + (Sinh(x) + Cosh(y))/Sinh(x + y + z);
		}

		/// <summary>
		/// Function performing some operations with square root
		/// </summary>
		template <class R>
		static R sqrt_func(const R& x, const R& y, const R& z)
		{
			return Sqrt(2*x * x + 1) + Sqrt(9 * y * y - 6 * y * z + z * z + 1);
		}

		template <class R>
		void PerformTest(const std::function<R(R, R, R)>& funcReal,
			const std::function<dual <R, 3>(dual <R, 3>, dual <R, 3>, dual <R, 3>)>& funcDual, R tolerance)
		{
			const R x = get_random_double();
			const R y = get_random_double();
			const R z = get_random_double();

			const double h = 1e-7;

			const auto dx_reference = (funcReal(x + h, y, z) - funcReal(x - h, y, z)) / (2 * h);
			const auto dy_reference = (funcReal(x, y + h, z) - funcReal(x, y - h, z)) / (2 * h);
			const auto dz_reference = (funcReal(x, y, z + h) - funcReal(x, y, z - h)) / (2 * h);

			Logger::WriteMessage((std::string("dx_reference = ") + std::to_string(dx_reference) + "\n").c_str());
			Logger::WriteMessage((std::string("dy_reference = ") + std::to_string(dy_reference) + "\n").c_str());
			Logger::WriteMessage((std::string("dz_reference = ") + std::to_string(dz_reference) + "\n").c_str());

			const auto result_dual = funcDual(
				dual <R, 3>{x, { 1, 0, 0 }},
				dual <R, 3>{y, { 0, 1, 0 }},
				dual <R, 3>{z, { 0, 0, 1 }});

			Logger::WriteMessage((std::string("dx_calc = ") + std::to_string(result_dual.Dual()[0]) + "\n").c_str());
			Logger::WriteMessage((std::string("dy_calc = ") + std::to_string(result_dual.Dual()[1]) + "\n").c_str());
			Logger::WriteMessage((std::string("dz_calc = ") + std::to_string(result_dual.Dual()[2]) + "\n").c_str());

			const auto dx_diff = std::abs(dx_reference - result_dual.Dual()[0]);
			const auto dy_diff = std::abs(dy_reference - result_dual.Dual()[1]);
			const auto dz_diff = std::abs(dz_reference - result_dual.Dual()[2]);

			Logger::WriteMessage((std::string("dx_diff = ") + auxutils::to_string_with_precision(dx_diff) + "\n").c_str());
			Logger::WriteMessage((std::string("dy_diff = ") + auxutils::to_string_with_precision(dy_diff) + "\n").c_str());
			Logger::WriteMessage((std::string("dz_diff = ") + auxutils::to_string_with_precision(dz_diff) + "\n").c_str());

			Assert::IsTrue(dx_diff < tolerance, L"Too big deviation for the partial derivative with respect to x");
			Assert::IsTrue(dy_diff < tolerance, L"Too big deviation for the partial derivative with respect to y");
			Assert::IsTrue(dz_diff < tolerance, L"Too big deviation for the partial derivative with respect to z");
		}
	public:

		TEST_CLASS_INITIALIZE(Initialize)
		{
			std::srand(std::time(nullptr));
		}
		TEST_METHOD(SimpleArithmeticsTest)
		{
			PerformTest<double>(arithm_func<double>, arithm_func<dual <double, 3>>, 1e-7);
		}

		TEST_METHOD(TrigonometryTest)
		{
			PerformTest<double>(trigonometric_func<double>, trigonometric_func<dual <double, 3>>, 1e-7);
		}

		TEST_METHOD(ExponentTest)
		{
			PerformTest<double>(exponent_func<double>, exponent_func<dual <double, 3>>, 6e-5);
		}

		TEST_METHOD(SquareRootTest)
		{
			PerformTest<double>(sqrt_func<double>, sqrt_func<dual <double, 3>>, 1e-7);
		}
	};
}