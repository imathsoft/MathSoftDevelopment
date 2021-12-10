
#include "CppUnitTest.h"
#include "UnitTestAux.h"
#include "../BVP/Utils/AuxUtils.h"
#include <functional>
#include "../BVP/FunctionApproximation/DerivativeEvaluator/DifferentiableFunction.h"
#include <algorithm>
#include <numeric>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{
	struct testFunction
	{
		template <class T>
		T operator () (const T& x)
		{
			return (4 + 2*x*x*4 + 6)*(7.0 + 3.0*x*x*5.0 + 10.0)*(9 - (6/(x*x*x))/7 - 11)*(2.0 - (3.0/(x*x))/4.0 - 1.0)/(1.2*x*x*x*x - 3.5*x*x*x + 7.7*x*x + 6.8*x + 9.5);
		}
	};

	struct testFunctionDerivative
	{
		template <class T>
		T operator () (const T& x)
		{
			T t1 = x * x;
			T t3 = 0.170e2 + 0.1500e2 * t1;
			T t5 = t1 * x;
			T t6 = 0.1e1 / t5;
			T t8 = -0.2e1 - 0.85714285714285714286e0 * t6;
			T t11 = 0.10e1 - 0.75000000000000000000e0 / t1;
			T t13 = t1 * t1;
			T t18 = 0.12e1 * t13 - 0.35e1 * t5 + 0.77e1 * t1 + 0.68e1 * x + 0.95e1;
			T t19 = 0.1e1 / t18;
			T t20 = t8 * t11 * t19;
			T t24 = 0.8e1 * t1 + 0.10e2;
			T t28 = t24 * t3;
			T t39 = t18 * t18;
			return 0.16e2 * x * t3 * t20 + 0.3000e2 * t24 * x * t20 + 0.25714285714285714286e1 * t28 / t13 * t11 * t19 + 0.15000000000000000000e1 * t28 * t8 * t6 * t19 - 0.1e1 * t28 * t8 * t11 / t39 * (0.48e1 * t5 - 0.105e2 * t1 + 0.154e2 * x + 0.68e1);
		}
	};

	struct testFunctionSecondDerivative
	{
		template <class T>
		T operator () (const T& x)
		{
			T t1 = x * x;
			T t3 = 0.170e2 + 0.1500e2 * t1;
			T t5 = t1 * x;
			T t6 = 0.1e1 / t5;
			T t8 = -0.2e1 - 0.85714285714285714286e0 * t6;
			T t10 = 0.1e1 / t1;
			T t12 = 0.10e1 - 0.75000000000000000000e0 * t10;
			T t13 = t1 * t1;
			T t18 = 0.12e1 * t13 - 0.35e1 * t5 + 0.77e1 * t1 + 0.68e1 * x + 0.95e1;
			T t19 = t18 * t18;
			T t20 = 0.1e1 / t19;
			T t21 = t12 * t20;
			T t25 = 0.48e1 * t5 - 0.105e2 * t1 + 0.154e2 * x + 0.68e1;
			T t26 = t21 * t25;
			T t30 = 0.8e1 * t1 + 0.10e2;
			T t35 = t30 * t3;
			T t39 = 0.1e1 / t18;
			T t43 = 0.1e1 / t13;
			T t51 = t35 * t8;
			T t65 = t25 * t25;
			T t70 = t12 * t39;
			T t83 = t8 * t39;
			return -0.32e2 * x * t3 * t8 * t26 - 0.6000e2 * t30 * x * t8 * t26 - 0.10285714285714285714e2 * t35 / t13 / x * t12 * t39 - 0.51428571428571428571e1 * t35 * t43 * t26 - 0.45000000000000000000e1 * t35 * t8 * t43 * t39 - 0.30000000000000000000e1 * t51 * t6 * t20 * t25 - 0.1e1 * t51 * t21 * (0.144e2 * t1 - 0.210e2 * x + 0.154e2) + 0.2e1 * t51 * t12 / t19 / t18 * t65 + 0.16e2 * t3 * t8 * t70 + 0.3000e2 * t30 * t8 * t70 + 0.96000e3 * t1 * t8 * t70 + 0.82285714285714285714e2 * t6 * t3 * t70 + 0.48000000000000000000e2 * t10 * t3 * t83 + 0.15428571428571428571e3 * t30 * t6 * t70 + 0.90000000000000000000e2 * t30 * t10 * t83 + 0.77142857142857142858e1 * t35 / t13 / t5 * t39;
		}
	};

	struct testFunctionThirdDerivative
	{
		template <class T>
		T operator () (const T& x)
		{
			T t1 = x * x;
			T t2 = t1 * x;
			T t3 = 0.1e1 / t2;
			T t5 = -0.2e1 - 0.85714285714285714286e0 * t3;
			T t7 = 0.1e1 / t1;
			T t9 = 0.10e1 - 0.75000000000000000000e0 * t7;
			T t10 = t1 * t1;
			T t15 = 0.12e1 * t10 - 0.35e1 * t2 + 0.77e1 * t1 + 0.68e1 * x + 0.95e1;
			T t16 = 0.1e1 / t15;
			T t17 = t9 * t16;
			T t21 = 0.8e1 * t1 + 0.10e2;
			T t22 = 0.1e1 / t10;
			T t27 = 0.170e2 + 0.1500e2 * t1;
			T t28 = t21 * t27;
			T t29 = t10 * t10;
			T t34 = t3 * t27;
			T t35 = t5 * t16;
			T t38 = t21 * t3;
			T t48 = 0.1e1 / t10 / t1;
			T t59 = t28 * t22;
			T t60 = t15 * t15;
			T t61 = 0.1e1 / t60;
			T t62 = t9 * t61;
			T t65 = 0.144e2 * t1 - 0.210e2 * x + 0.154e2;
			T t66 = t62 * t65;
			T t69 = t28 * t5;
			T t74 = 0.48e1 * t2 - 0.105e2 * t1 + 0.154e2 * x + 0.68e1;
			T t88 = 0.1e1 / t60 / t15;
			T t89 = t9 * t88;
			T t90 = t74 * t74;
			T t91 = t89 * t90;
			T t98 = 0.288000e4 * x * t5 * t17 - 0.69428571428571428570e3 * t21 * t22 * t17 - 0.81000000000000000001e2 * t28 / t29 * t16 - 0.14400000000000000000e3 * t34 * t35 - 0.27000000000000000000e3 * t38 * t35 - 0.37028571428571428571e3 * t22 * t27 * t17 + 0.74057142857142857142e4 * t7 * t9 * t16 + 0.37028571428571428571e3 * t48 * t27 * t16 + 0.43200000000000000000e4 / x * t5 * t16 + 0.69428571428571428570e3 * t21 * t48 * t16 - 0.77142857142857142857e1 * t59 * t66 + 0.13500000000000000000e2 * t69 * t22 * t61 * t74 - 0.45000000000000000000e1 * t69 * t3 * t61 * t65 - 0.1e1 * t69 * t62 * (0.288e2 * x - 0.210e2) + 0.15428571428571428571e2 * t59 * t91 + 0.90000000000000000000e1 * t69 * t3 * t88 * t90;
			T t100 = t21 * x * t5;
			T t104 = x * t27 * t5;
			T t107 = t60 * t60;
			T t119 = 0.1e1 / t10 / x;
			T t121 = t62 * t74;
			T t142 = t5 * t61 * t74;
			T t163 = 0.18000e3 * t100 * t91 + 0.96e2 * t104 * t91 - 0.6e1 * t69 * t9 / t107 * t90 * t74 - 0.48e2 * t104 * t66 - 0.9000e2 * t100 * t66 + 0.30857142857142857143e2 * t28 * t119 * t121 - 0.24685714285714285714e3 * t34 * t121 - 0.9000e2 * t21 * t5 * t121 + 0.51428571428571428571e2 * t28 * t48 * t9 * t16 + 0.18000000000000000000e2 * t28 * t5 * t119 * t16 - 0.288000e4 * t1 * t5 * t121 - 0.14400000000000000000e3 * t7 * t27 * t142 - 0.46285714285714285713e3 * t38 * t121 - 0.27000000000000000000e3 * t21 * t7 * t142 - 0.23142857142857142857e2 * t28 / t10 / t2 * t61 * t74 - 0.48e2 * t27 * t5 * t121 + 0.6e1 * t69 * t89 * t65 * t74;
			return t98 + t163;
		}
	};


	TEST_CLASS(DifferentiableFunctionTest)
	{
	private:
		template <class T, int Size>
		T GetMaxAbsDeviationfromValue(std::array<T, Size> arr, const T value)
		{
			T maxDiff = 0.0;
			std::for_each(arr.begin(), arr.end(), [&](T& d) { maxDiff = std::max( maxDiff, std::abs(d - value));});

			return maxDiff;
		}

	public:
		TEST_METHOD(DifferentiableFunctionFullTest)
		{

			auto func = DifferentiableFunctionFactory<double>(testFunction());
			std::array<double, 4> maxRelDiff; maxRelDiff.fill(0.0);

			std::array<double, 10> arguments = { -5.6, 3, 11.32, 1, 0.3, 3.5, -3.78, -20.1, 15.4, -7.777};

			for (int iteration = 0; iteration < arguments.size(); iteration ++)
			{
				std::array<double, 4> result;
				double argument = arguments[iteration];
				result[0] = func.Evaluate(argument, result[1], result[2], result[3]);

				std::array<double, 4> resultRef = {testFunction()(argument), 
					testFunctionDerivative()(argument), 
					testFunctionSecondDerivative()(argument), testFunctionThirdDerivative()(argument)};
				for (int i = 0; i < 4; i++)
				{
					maxRelDiff[i] = std::max(std::abs(result[i] - resultRef[i])/std::max<double>(1, std::abs(resultRef[i])),maxRelDiff[i]);
				}
			}

			Assert::IsTrue(maxRelDiff[0] <= 5e-16,  UnitTestAux::Message("Too big deviation for the function value : " + auxutils::ToString(maxRelDiff[0])));
			Assert::IsTrue(maxRelDiff[1] <= 6e-15,  UnitTestAux::Message("Too big deviation for the first derivative value : " + auxutils::ToString(maxRelDiff[1])));
			Assert::IsTrue(maxRelDiff[2] <= 2e-14,  UnitTestAux::Message("Too big deviation for the second derivative value : " + auxutils::ToString(maxRelDiff[2])));
			Assert::IsTrue(maxRelDiff[3] <= 3e-14,  UnitTestAux::Message("Too big deviation for the third derivative value : " + auxutils::ToString(maxRelDiff[3])));
		}

		TEST_METHOD(DifferentiableFunctionSyncTest)
		{
			auto func = DifferentiableFunctionFactory<double>(testFunction());
			std::array<double, 4> values;
			std::array<double, 3> derivatives;
			std::array<double, 2> secondDerivatives;
			std::array<double, 1> thirdDerivatives;

			std::array<double, 10> arguments = { -5.6, 3, 11.32, 1, 0.3, 3.5, -3.78, -20.1, 15.4, -7.777};

			for (int iteration = 0; iteration < arguments.size(); iteration ++)
			{
				double argument = arguments[iteration];
				values[0] = func.Evaluate(argument, derivatives[0], secondDerivatives[0], thirdDerivatives[0]);
				values[1] = func.Evaluate(argument, derivatives[1], secondDerivatives[1]);
				values[2] = func.Evaluate(argument, derivatives[2]);
				values[3] = func.Evaluate(argument);

				double maxDiffValues = GetMaxAbsDeviationfromValue(values, values[0])/std::abs(values[0]);
				double maxDiffDerivatives = GetMaxAbsDeviationfromValue(derivatives, derivatives[0])/std::abs(derivatives[0]);
				double maxDiffSecDerivatives = GetMaxAbsDeviationfromValue(secondDerivatives, secondDerivatives[0])/std::abs(secondDerivatives[0]);

				Assert::IsTrue(maxDiffValues < 5e-16, UnitTestAux::Message("Too big deviation for the function value : " + auxutils::ToString(maxDiffValues)));
				Assert::IsTrue(maxDiffDerivatives <= 0.0, UnitTestAux::Message("Too big deviation for the derivative value : " + auxutils::ToString(maxDiffDerivatives)));
				Assert::IsTrue(maxDiffSecDerivatives <= 0.0, UnitTestAux::Message("Too big deviation for the second derivative value : " + auxutils::ToString(maxDiffSecDerivatives)));
			}
		}
	};
}