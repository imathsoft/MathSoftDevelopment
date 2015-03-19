
#include <mpreal.h>
#include "CppUnitTest.h"
#include "..\BVP\FunctionApproximation\X_Function.h"
#include "UnitTestAux.h"

using namespace mpfr;

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
	};
}