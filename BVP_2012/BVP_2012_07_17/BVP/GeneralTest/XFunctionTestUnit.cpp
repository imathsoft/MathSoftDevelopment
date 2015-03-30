
#include <mpreal.h>
#include "CppUnitTest.h"
#include "..\BVP\FunctionApproximation\X_Function.h"
#include "UnitTestAux.h"

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

		TEST_METHOD(StraightXFunctionTestMethod)
		{
		  mpreal::set_default_prec(128);
		  auto invertFunc = X3_Func<mpreal>(10,-2,-30,14,"0.1",1e-36);

		  Assert::IsTrue(abs("-31.895112005308260774136632962176948035" - invertFunc.Derivative) < 
			  1e-36);
		  Assert::IsTrue(abs("10.890975443956699164876695239161704776" - invertFunc.Value) < 1e-36);
		}

		TEST_METHOD(StraightXFunctionDerivativesTestMethod)
		{
			const int c = 6;
			mpreal::set_default_prec(128);
			GradientVector<mpreal, c> X; X << "0.1", 1, 0, 0, 0, 0;
			GradientVector<mpreal, c> A; A << 10, 0, 1, 0, 0, 0;
			GradientVector<mpreal, c> B; B << -2, 0, 0, 1, 0, 0;
			GradientVector<mpreal, c> C; C << -30, 0, 0, 0, 1, 0;
			GradientVector<mpreal, c> D; D << 14, 0, 0, 0, 0, 1;

			auto invertFunc = X3_Func<GradientVector<mpreal, c>>(A, B, C, D, X, 1e-38);

			Assert::IsTrue(abs("10.8909754439566991648766952391617047763" - invertFunc.Value[0, 0]) < 
				1e-38, Message("Value Mismatch : " + 
				("10.8909754439566991648766952391617047763" - invertFunc.Value[0, 0]).toString()));

			Assert::IsTrue(abs("-31.8951120053082607741366329621769480352" - invertFunc.Value[0, 1]) < 
				1e-38, Message("Derivative Mismatch : " + 
				("-31.8951120053082607741366329621769480352" - invertFunc.Value[0, 1]).toString()));

			Assert::IsTrue(abs("0.00207594273574432631847956266999583199" - invertFunc.Value[0, 2]) < 
				1e-38, Message("dA Mismatch : " + 
				("0.00207594273574432631847956266999583199" - invertFunc.Value[0, 2]).toString()));

			Assert::IsTrue(abs("0.06482097233719439617574302934815502327" - invertFunc.Value[0, 3]) < 
				1e-38, Message("dB Mismatch : " + 
				("0.06482097233719439617574302934815502327" - invertFunc.Value[0, 3]).toString()));

			Assert::IsTrue(abs("0.0997501864418126094088359352974983223384" - invertFunc.Value[0, 4]) < 
				1e-38, Message("dC Mismatch : " + 
				("0.0997501864418126094088359352974983223384" - invertFunc.Value[0, 4]).toString()));

			Assert::IsTrue(abs("0.991677216943648389081555235577618174747" - invertFunc.Value[0, 5]) < 
				1e-38, Message("dD Mismatch : " + 
				("0.991677216943648389081555235577618174747" - invertFunc.Value[0, 5]).toString()));

			Assert::IsTrue(abs("-10.8909754439566991648766952391617047763" - invertFunc.Derivative[0, 1]) < 
				1e-37, Message("ddX Mismatch : " + 
				("-10.8909754439566991648766952391617047763" - invertFunc.Derivative[0, 1]).toString()));

			Assert::IsTrue(abs("0.059648736118014891805556428787273815261" - invertFunc.Derivative[0, 2]) < 
				1e-38, Message("dXdA Mismatch : " + 
				("0.059648736118014891805556428787273815261" - invertFunc.Derivative[0, 2]).toString()));

			Assert::IsTrue(abs("1.24335752123318609230345690920365851415" - invertFunc.Derivative[0, 3]) < 
				1e-37, Message("dXdB Mismatch : " + 
				("1.24335752123318609230345690920365851415" - invertFunc.Derivative[0, 3]).toString()));

			Assert::IsTrue(abs("0.993341384943721102519238923763223297139" - invertFunc.Derivative[0, 4]) < 
				1e-38, Message("dXdC Mismatch : " + 
				("0.993341384943721102519238923763223297139" - invertFunc.Derivative[0, 4]).toString()));

			Assert::IsTrue(abs("-0.149633604071187692754247517805732080072" - invertFunc.Derivative[0, 5]) < 
				1e-38, Message("dXdD Mismatch : " + 
				("-0.149633604071187692754247517805732080072" - invertFunc.Derivative[0, 5]).toString()));
		}
	};
}