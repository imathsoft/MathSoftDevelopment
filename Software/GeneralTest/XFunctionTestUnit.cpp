
#include "CppUnitTest.h"
#include "..\BVP\FunctionApproximation\X_Function.h"
#include "UnitTestAux.h"
#include <numeric>
#include "../BVP/Utils/AuxUtils.h"

using namespace UnitTestAux;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

typedef number<cpp_dec_float<40>, et_off> multiprec_float;

namespace GeneralTest
{		
	TEST_CLASS(XFunctionUnitTest)
	{
	public:
		
		TEST_METHOD(InvertXFunctionTestMethod)
		{
		  auto invertFunc = XI_Func<multiprec_float>(10,-2,-30,14, multiprec_float("0.1"), 1e-36);

		  auto deriv_mismatch = abs(multiprec_float("-25.821239292751734216871012936299241006") - invertFunc.Derivative);
		  Assert::IsTrue(deriv_mismatch < 1e-36, Message("Derivative mismatch : " + auxutils::ToString(deriv_mismatch)));

		  auto value_mismatch = abs(multiprec_float("11.237250563288499476414732181933848865") - invertFunc.Value);
		  Assert::IsTrue(value_mismatch < 1e-36, Message("Value mismatch : " + auxutils::ToString(value_mismatch)));
		}

		TEST_METHOD(StraightX3FunctionTestMethod)
		{
		  auto streightFunc = X3_Func<multiprec_float>(10,-2,-30,14,multiprec_float("0.1"),1e-36);

		  auto deriv_mismatch = abs(multiprec_float("-31.895112005308260774136632962176948035") - streightFunc.Derivative);
		  Assert::IsTrue(deriv_mismatch < 1e-36, Message("Derivative mismatch : " + auxutils::ToString(deriv_mismatch)));

		  auto value_mismatch = abs(multiprec_float("10.890975443956699164876695239161704776") - streightFunc.Value);
		  Assert::IsTrue(value_mismatch < 1e-36, Message("Value mismatch : " + auxutils::ToString(value_mismatch)));
		}

		TEST_METHOD(StraightX4FunctionTestMethod)
		{
		  auto streightFunc = X4_Func<multiprec_float>(10,-2,-30, 14, 0, 0, multiprec_float("0.1"),1e-36);

		  auto deriv_mismatch = abs(multiprec_float("-31.895112005308260774136632962176948035") - streightFunc.Derivative);
		  Assert::IsTrue(deriv_mismatch < 1e-36, Message("Derivative mismatch : " + auxutils::ToString(deriv_mismatch)));

		  auto value_mismatch = abs(multiprec_float("10.890975443956699164876695239161704776") - streightFunc.Value);
		  Assert::IsTrue(value_mismatch < 1e-36, Message("Value mismatch : " + auxutils::ToString(value_mismatch)));
		}

		TEST_METHOD(StraightX4FunctionNonuniformTestMethod)
		{
		  auto streightFunc = X4_Func<multiprec_float>(10,-2,-30, 14, 13, 17, multiprec_float("0.1"),1e-36);

		  auto deriv_mismatch = abs(multiprec_float("-30.133716147148905398140439228392955594") - streightFunc.Derivative);
		  Assert::IsTrue(deriv_mismatch < 1e-36, Message("Derivative mismatch : " + auxutils::ToString(deriv_mismatch)));
		  auto value_mismatch = abs(multiprec_float("10.978041548052233737621628620689023378") - streightFunc.Value);
		  Assert::IsTrue(value_mismatch < 1e-36, Message("Value mismatch : " + auxutils::ToString(value_mismatch)));
		  auto second_deriv_mismatch = abs(multiprec_float("7.32195845194776626237837137931097662187") - streightFunc.SecDerivative);
		  Assert::IsTrue(second_deriv_mismatch < 1e-36, Message("Second derivative mismatch : " + auxutils::ToString(second_deriv_mismatch)));
		}

		TEST_METHOD(StraightXFunctionDerivativesTestMethod)
		{
			GVTypes<multiprec_float>::FullGradientVector X; X << multiprec_float("0.1") << 1 << 0 << 0 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector A; A << 10 << 0 << 1 << 0 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector B; B << -2 << 0 << 0 << 1 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector C; C << -30 << 0 << 0 << 0 << 1 << 0;
			GVTypes<multiprec_float>::FullGradientVector D; D << 14 << 0 << 0 << 0 << 0 << 1;

			auto straightFunc = X3_Func<GVTypes<multiprec_float>::FullGradientVector>(A, B, C, D, X, 1e-38);

			auto f_mismatch = abs(multiprec_float("10.8909754439566991648766952391617047763") - straightFunc.Value[0, 0]);
			Assert::IsTrue(f_mismatch < 2e-38, Message("Value Mismatch : " + auxutils::ToString(f_mismatch)));

			auto dfdx_mismatch = abs(multiprec_float("-31.8951120053082607741366329621769480352") - straightFunc.Value[0, 1]);
			Assert::IsTrue(dfdx_mismatch < 6e-38, Message("Derivative Mismatch : " + auxutils::ToString(dfdx_mismatch)));

			auto dfdA_mismatch = abs(multiprec_float("0.00207594273574432631847956266999583199") - straightFunc.Value[0, 2]);
			Assert::IsTrue(dfdA_mismatch < 1e-38, Message("dA Mismatch : " + auxutils::ToString(dfdA_mismatch)));

			auto dfdB_mismatch = abs(multiprec_float("0.06482097233719439617574302934815502327") - straightFunc.Value[0, 3]);
			Assert::IsTrue(dfdB_mismatch < 1e-38, Message("dB Mismatch : " + auxutils::ToString(dfdB_mismatch)));

			auto dfdC_mismatch = abs(multiprec_float("0.0997501864418126094088359352974983223384") - straightFunc.Value[0, 4]);
			Assert::IsTrue(dfdC_mismatch < 1e-38, Message("dC Mismatch : " + auxutils::ToString(dfdC_mismatch)));

			auto dfdD_mismatch = abs(multiprec_float("0.991677216943648389081555235577618174747") - straightFunc.Value[0, 5]);
			Assert::IsTrue(dfdD_mismatch < 1e-38, Message("dD Mismatch : " + auxutils::ToString(dfdD_mismatch)));

			auto dfdx_mismatch1 = abs(straightFunc.Value[0, 1] - straightFunc.Derivative[0, 0]);
			Assert::IsTrue(dfdx_mismatch1 < 1e-38, Message("dX Mismatch : " + auxutils::ToString(dfdx_mismatch1)));

			auto dfddx_mismatch = abs(multiprec_float("-10.8909754439566991648766952391617047763") - straightFunc.Derivative[0, 1]);
			Assert::IsTrue(dfddx_mismatch < 1e-37, Message("ddX Mismatch : " + auxutils::ToString(dfddx_mismatch)));

			auto dfdXdA_mismatch = abs(multiprec_float("0.059648736118014891805556428787273815261") - straightFunc.Derivative[0, 2]);
			Assert::IsTrue(dfdXdA_mismatch < 1e-38, Message("dXdA Mismatch : " + auxutils::ToString(dfdXdA_mismatch)));

			auto dfdXdB_mismatch = abs(multiprec_float("1.24335752123318609230345690920365851415") - straightFunc.Derivative[0, 3]);
			Assert::IsTrue(dfdXdB_mismatch < 1e-37, Message("dXdB Mismatch : " + auxutils::ToString(dfdXdB_mismatch)));

			auto dfdXdC_mismatch = abs(multiprec_float("0.993341384943721102519238923763223297139") - straightFunc.Derivative[0, 4]);
			Assert::IsTrue(dfdXdC_mismatch < 1e-38, Message("dXdC Mismatch : " + auxutils::ToString(dfdXdC_mismatch)));

			auto dfdXdD_mismatch = abs(multiprec_float("-0.149633604071187692754247517805732080072") - straightFunc.Derivative[0, 5]);
			Assert::IsTrue(dfdXdD_mismatch < 1e-38, Message("dXdD Mismatch : " + auxutils::ToString(dfdXdD_mismatch)));

			auto dfddX_mismatch1 = abs(straightFunc.SecDerivative[0, 0] - straightFunc.Derivative[0, 1]);
			Assert::IsTrue(dfddX_mismatch1 < 1e-37, Message("dXdX Mismatch : " + auxutils::ToString(dfddX_mismatch1)));
		}

		TEST_METHOD(InverseXFunctionDerivativesTestMethod)
		{
			GVTypes<multiprec_float>::FullGradientVector X; X << multiprec_float("0.1") << 1 << 0 << 0 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector A; A << 10 << 0 << 1 << 0 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector B; B << -2 << 0 << 0 << 1 << 0 << 0;
			GVTypes<multiprec_float>::FullGradientVector C; C << -30 << 0 << 0 << 0 << 1 << 0;
			GVTypes<multiprec_float>::FullGradientVector D; D << 14 << 0 << 0 << 0 << 0 << 1;

			auto invertFunc = XI_Func<GVTypes<multiprec_float>::FullGradientVector>(A, B, C, D, X, 1e-38);

			auto f_mismatch = abs(multiprec_float("11.2372505632884994764147321819338488645") - invertFunc.Value[0, 0]);
			Assert::IsTrue(f_mismatch < 1e-38, Message("Value Mismatch : " + auxutils::ToString(f_mismatch)));

			auto dfdX_mismatch = abs(multiprec_float("-25.8212392927517342168710129362992410056") - invertFunc.Value[0, 1]);
			Assert::IsTrue(dfdX_mismatch < 2e-37, Message("Derivative Mismatch : " + auxutils::ToString(dfdX_mismatch)));

			auto dfdA_mismatch = abs(multiprec_float("-0.0044361062899309975455071595025040810224") - invertFunc.Value[0, 2]);
			Assert::IsTrue(dfdA_mismatch < 1e-38, Message("dA Mismatch : " + auxutils::ToString(dfdA_mismatch)));

			auto dfdB_mismatch = abs(multiprec_float("-0.134673816617473526404154857243154327654") - invertFunc.Value[0, 3]);
			Assert::IsTrue(dfdB_mismatch < 1e-37, Message("dB Mismatch : " + auxutils::ToString(dfdB_mismatch)));

			auto dfdC_mismatch = abs(multiprec_float("0.0920916478903833507861755939355383711840") - invertFunc.Value[0, 4]);
			Assert::IsTrue(dfdC_mismatch < 1e-38, Message("dC Mismatch : " + auxutils::ToString(dfdC_mismatch)));

			auto dfdD_mismatch = abs(multiprec_float("1.0") - invertFunc.Value[0, 5]);
			Assert::IsTrue(dfdD_mismatch < 1e-38, Message("dD Mismatch : " + auxutils::ToString(dfdD_mismatch)));

			auto dfdX_mismatch1 = abs(invertFunc.Value[0, 1] - invertFunc.Derivative[0, 0]);
			Assert::IsTrue(dfdX_mismatch1 < 1e-38, Message("dX Mismatch : " + auxutils::ToString(dfdX_mismatch1)));

			auto dfddX_mismatch = abs(multiprec_float("25.8212392927517342168710129362992410056") - invertFunc.Derivative[0, 1]);
			Assert::IsTrue(dfddX_mismatch < 2e-37, Message("ddX Mismatch : " + auxutils::ToString(dfddX_mismatch)));

			auto dfdXdA_mismatch = abs(multiprec_float("-0.129106196463758671084355064681496205028") - invertFunc.Derivative[0, 2]);
			Assert::IsTrue(dfdXdA_mismatch < 1e-38, Message("dXdA Mismatch : " + auxutils::ToString(dfdXdA_mismatch)));

			auto dfdXdB_mismatch = abs(multiprec_float("-2.58212392927517342168710129362992410056") - invertFunc.Derivative[0, 3]);
			Assert::IsTrue(dfdXdB_mismatch < 1e-37, Message("dXdB Mismatch : " + auxutils::ToString(dfdXdB_mismatch)));

			auto dfdXdC_mismatch = abs(multiprec_float("0.860707976425057807229033764543308033519") - invertFunc.Derivative[0, 4]);
			Assert::IsTrue(dfdXdC_mismatch < 1e-38, Message("dXdC Mismatch : " + auxutils::ToString(dfdXdC_mismatch)));

			auto dfdXdD_mismatch = abs(0 - invertFunc.Derivative[0, 5]);
			Assert::IsTrue(dfdXdD_mismatch < 1e-38, Message("dXdD Mismatch : " + auxutils::ToString(dfdXdD_mismatch)));

			auto dfddX_mismatch1 = abs(invertFunc.SecDerivative[0, 0] - invertFunc.Derivative[0, 1]);
			Assert::IsTrue(dfddX_mismatch1 < 1e-38, Message("dXdX Mismatch : " + auxutils::ToString(dfddX_mismatch1)));
		}

		template <class T>
		T XI_Bernoulli_check_function(const T A, const T B, const T C, const T D, const T h)
		{
			return log((h*A*C+B*C+auxutils::Sqrt(A*A*C*C*h*h+2*A*B*C*C*h+A))/(B*C+auxutils::Sqrt(A)))/auxutils::Sqrt(A) + D;
		}

		template <class T>
		T XI_Bernoulli_derivative_check_function(const T A, const T B, const T C, const T h)
		{
			return 1/(auxutils::Sqrt(1/(C*C) - A*h*h - 2*B*h));
		}

		TEST_METHOD(XI_Bernoulli_test)
		{
			//Arrange
			double A = 10;
			double B = 100;
			double C = 0.5;
			double D = 3;
			double h = 0.015;
			auto reference_func = XI_Bernoulli_check_function(A, B, C, D, h);
			auto reference_deriv = XI_Bernoulli_derivative_check_function(-A, -B, C, h);
			auto reference_deriv_sec = -(A*h + B)*reference_deriv*reference_deriv*reference_deriv;

			//Act
			auto result = XI_Bernoulli<double>(-A, -B, C, D, h, std::numeric_limits<double>::epsilon());

			//Assert
			Assert::IsTrue(abs((result.Value - reference_func)/reference_func) <= std::numeric_limits<double>::epsilon() && 
						   abs((result.Derivative - reference_deriv)/reference_deriv) <= std::numeric_limits<double>::epsilon() && 
						   abs((result.SecDerivative - reference_deriv_sec)/reference_deriv_sec) <= 3*std::numeric_limits<double>::epsilon(),
						   L"Mismatch between calculated and reference values");
		}
	};
}