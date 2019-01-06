
#include <mpreal.h>
#include "CppUnitTest.h"
#include "..\BVP\FunctionApproximation\X_Function.h"
#include "UnitTestAux.h"
#include <numeric>

using namespace mpfr;
using namespace UnitTestAux;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{		
	TEST_CLASS(XFunctionUnitTest)
	{
	public:
		
		TEST_METHOD(InvertXFunctionTestMethod)
		{
		  mpreal::set_default_prec(128);
		  auto invertFunc = XI_Func<mpreal>(10,-2,-30,14,"0.1",1e-36);

		  Assert::IsTrue(abs("-25.821239292751734216871012936299241006" - invertFunc.Derivative) < 
			  1e-36);
		  Assert::IsTrue(abs("11.237250563288499476414732181933848865" - invertFunc.Value) < 1e-36);
		}

		TEST_METHOD(StraightX3FunctionTestMethod)
		{
		  mpreal::set_default_prec(128);

		  auto streightFunc = X3_Func<mpreal>(10,-2,-30,14,"0.1",1e-36);

		  Assert::IsTrue(abs("-31.895112005308260774136632962176948035" - streightFunc.Derivative) < 
			  1e-36);
		  Assert::IsTrue(abs("10.890975443956699164876695239161704776" - streightFunc.Value) < 1e-36);
		}

		TEST_METHOD(StraightX4FunctionTestMethod)
		{
		  mpreal::set_default_prec(128);

		  auto streightFunc = X4_Func<mpreal>(10,-2,-30, 14, 0, 0, "0.1",1e-36);

		  Assert::IsTrue(abs("-31.895112005308260774136632962176948035" - streightFunc.Derivative) < 
			  1e-36);
		  Assert::IsTrue(abs("10.890975443956699164876695239161704776" - streightFunc.Value) < 1e-36);
		}

		TEST_METHOD(StraightX4FunctionNonuniformTestMethod)
		{
		  mpreal::set_default_prec(128);

		  auto streightFunc = X4_Func<mpreal>(10,-2,-30, 14, 13, 17, "0.1",1e-36);

		  Assert::IsTrue(abs("-30.133716147148905398140439228392955594" - streightFunc.Derivative) < 1e-36);
		  Assert::IsTrue(abs( "10.978041548052233737621628620689023378" - streightFunc.Value) < 1e-36);
		  Assert::IsTrue(abs( "7.32195845194776626237837137931097662187" - streightFunc.SecDerivative) < 1e-36);
		}

		TEST_METHOD(StraightXFunctionDerivativesTestMethod)
		{
			mpreal::set_default_prec(128);
			GVTypes<mpreal>::FullGradientVector X; X << "0.1" << 1 << 0 << 0 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector A; A << 10 << 0 << 1 << 0 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector B; B << -2 << 0 << 0 << 1 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector C; C << -30 << 0 << 0 << 0 << 1 << 0;
			GVTypes<mpreal>::FullGradientVector D; D << 14 << 0 << 0 << 0 << 0 << 1;

			auto straightFunc = X3_Func<GVTypes<mpreal>::FullGradientVector>(A, B, C, D, X, 1e-38);

			Assert::IsTrue(abs("10.8909754439566991648766952391617047763" - straightFunc.Value[0, 0]) < 
				1e-38, Message("Value Mismatch : " + 
				("10.8909754439566991648766952391617047763" - straightFunc.Value[0, 0]).toString()));

			Assert::IsTrue(abs("-31.8951120053082607741366329621769480352" - straightFunc.Value[0, 1]) < 
				1e-38, Message("Derivative Mismatch : " + 
				("-31.8951120053082607741366329621769480352" - straightFunc.Value[0, 1]).toString()));

			Assert::IsTrue(abs("0.00207594273574432631847956266999583199" - straightFunc.Value[0, 2]) < 
				1e-38, Message("dA Mismatch : " + 
				("0.00207594273574432631847956266999583199" - straightFunc.Value[0, 2]).toString()));

			Assert::IsTrue(abs("0.06482097233719439617574302934815502327" - straightFunc.Value[0, 3]) < 
				1e-38, Message("dB Mismatch : " + 
				("0.06482097233719439617574302934815502327" - straightFunc.Value[0, 3]).toString()));

			Assert::IsTrue(abs("0.0997501864418126094088359352974983223384" - straightFunc.Value[0, 4]) < 
				1e-38, Message("dC Mismatch : " + 
				("0.0997501864418126094088359352974983223384" - straightFunc.Value[0, 4]).toString()));

			Assert::IsTrue(abs("0.991677216943648389081555235577618174747" - straightFunc.Value[0, 5]) < 
				1e-38, Message("dD Mismatch : " + 
				("0.991677216943648389081555235577618174747" - straightFunc.Value[0, 5]).toString()));

			Assert::IsTrue(abs(straightFunc.Value[0, 1] - straightFunc.Derivative[0, 0]) < 
				1e-38, Message("dX Mismatch : " + 
				(straightFunc.Value[0, 1] - straightFunc.Derivative[0, 0]).toString()));

			Assert::IsTrue(abs("-10.8909754439566991648766952391617047763" - straightFunc.Derivative[0, 1]) < 
				1e-37, Message("ddX Mismatch : " + 
				("-10.8909754439566991648766952391617047763" - straightFunc.Derivative[0, 1]).toString()));

			Assert::IsTrue(abs("0.059648736118014891805556428787273815261" - straightFunc.Derivative[0, 2]) < 
				1e-38, Message("dXdA Mismatch : " + 
				("0.059648736118014891805556428787273815261" - straightFunc.Derivative[0, 2]).toString()));

			Assert::IsTrue(abs("1.24335752123318609230345690920365851415" - straightFunc.Derivative[0, 3]) < 
				1e-37, Message("dXdB Mismatch : " + 
				("1.24335752123318609230345690920365851415" - straightFunc.Derivative[0, 3]).toString()));

			Assert::IsTrue(abs("0.993341384943721102519238923763223297139" - straightFunc.Derivative[0, 4]) < 
				1e-38, Message("dXdC Mismatch : " + 
				("0.993341384943721102519238923763223297139" - straightFunc.Derivative[0, 4]).toString()));

			Assert::IsTrue(abs("-0.149633604071187692754247517805732080072" - straightFunc.Derivative[0, 5]) < 
				1e-38, Message("dXdD Mismatch : " + 
				("-0.149633604071187692754247517805732080072" - straightFunc.Derivative[0, 5]).toString()));

			Assert::IsTrue(abs(straightFunc.SecDerivative[0, 0] - straightFunc.Derivative[0, 1]) < 
				1e-37, Message("dXdX Mismatch : " + 
				(straightFunc.SecDerivative[0, 0] - straightFunc.Derivative[0, 1]).toString()));
		}

		TEST_METHOD(InverseXFunctionDerivativesTestMethod)
		{
			mpreal::set_default_prec(128);
			GVTypes<mpreal>::FullGradientVector X; X << "0.1" << 1 << 0 << 0 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector A; A << 10 << 0 << 1 << 0 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector B; B << -2 << 0 << 0 << 1 << 0 << 0;
			GVTypes<mpreal>::FullGradientVector C; C << -30 << 0 << 0 << 0 << 1 << 0;
			GVTypes<mpreal>::FullGradientVector D; D << 14 << 0 << 0 << 0 << 0 << 1;

			auto invertFunc = XI_Func<GVTypes<mpreal>::FullGradientVector>(A, B, C, D, X, 1e-38);

			Assert::IsTrue(abs("11.2372505632884994764147321819338488645" - invertFunc.Value[0, 0]) < 
				1e-38, Message("Value Mismatch : " + 
				("11.2372505632884994764147321819338488645" - invertFunc.Value[0, 0]).toString()));

			Assert::IsTrue(abs("-25.8212392927517342168710129362992410056" - invertFunc.Value[0, 1]) < 
				1e-37, Message("Derivative Mismatch : " + 
				("-25.8212392927517342168710129362992410056" - invertFunc.Value[0, 1]).toString()));

			Assert::IsTrue(abs("-0.0044361062899309975455071595025040810224" - invertFunc.Value[0, 2]) < 
				1e-38, Message("dA Mismatch : " + 
				("-0.0044361062899309975455071595025040810224" - invertFunc.Value[0, 2]).toString()));

			Assert::IsTrue(abs("-0.134673816617473526404154857243154327654" - invertFunc.Value[0, 3]) < 
				1e-37, Message("dB Mismatch : " + 
				("-0.134673816617473526404154857243154327654" - invertFunc.Value[0, 3]).toString()));

			Assert::IsTrue(abs("0.0920916478903833507861755939355383711840" - invertFunc.Value[0, 4]) < 
				1e-38, Message("dC Mismatch : " + 
				("0.0920916478903833507861755939355383711840" - invertFunc.Value[0, 4]).toString()));

			Assert::IsTrue(abs("1.0" - invertFunc.Value[0, 5]) < 
				1e-38, Message("dD Mismatch : " + 
				("1.0" - invertFunc.Value[0, 5]).toString()));

			Assert::IsTrue(abs(invertFunc.Value[0, 1] - invertFunc.Derivative[0, 0]) < 
				1e-38, Message("dX Mismatch : " + 
				(invertFunc.Value[0, 1] - invertFunc.Derivative[0, 0]).toString()));

			Assert::IsTrue(abs("25.8212392927517342168710129362992410056" - invertFunc.Derivative[0, 1]) < 
				1e-37, Message("ddX Mismatch : " + 
				("25.8212392927517342168710129362992410056" - invertFunc.Derivative[0, 1]).toString()));

			Assert::IsTrue(abs("-0.129106196463758671084355064681496205028" - invertFunc.Derivative[0, 2]) < 
				1e-38, Message("dXdA Mismatch : " + 
				("-0.129106196463758671084355064681496205028" - invertFunc.Derivative[0, 2]).toString()));

			Assert::IsTrue(abs("-2.58212392927517342168710129362992410056" - invertFunc.Derivative[0, 3]) < 
				1e-37, Message("dXdB Mismatch : " + 
				("-2.58212392927517342168710129362992410056" - invertFunc.Derivative[0, 3]).toString()));

			Assert::IsTrue(abs("0.860707976425057807229033764543308033519" - invertFunc.Derivative[0, 4]) < 
				1e-38, Message("dXdC Mismatch : " + 
				("0.860707976425057807229033764543308033519" - invertFunc.Derivative[0, 4]).toString()));

			Assert::IsTrue(abs("0.0" - invertFunc.Derivative[0, 5]) < 
				1e-38, Message("dXdD Mismatch : " + 
				("0.0" - invertFunc.Derivative[0, 5]).toString()));

			Assert::IsTrue(abs(invertFunc.SecDerivative[0, 0] - invertFunc.Derivative[0, 1]) < 
				1e-38, Message("dXdX Mismatch : " + 
				(invertFunc.SecDerivative[0, 0] - invertFunc.Derivative[0, 1]).toString()));
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