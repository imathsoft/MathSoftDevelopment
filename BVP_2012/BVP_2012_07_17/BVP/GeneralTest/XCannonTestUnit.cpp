#include "CppUnitTest.h"

#include <mpreal.h>
#include "..\BVP\Cannon\XCannon.h"
#include "..\BVP\Problems\Troesch.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\Cannon\TroeschHybridCannon.h"
#include "UnitTestAux.h"

using namespace mpfr;
using namespace UnitTestAux;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{
	TEST_CLASS(XCannonTestUnit)
	{
	public:
		
		TEST_METHOD(XCanonStraightTestMethod)
		{
			 mpreal::set_default_prec(128);
			 const int l = 1;
             XCannon<mpreal> xc(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, "0.01", &Troesch<mpreal, l>::StepFunc, &Troesch<mpreal, l>::CheckFunc, 1e-40);

			 auto shootingResult = xc.Shoot("0", "1", "0", "0.1");		
		
			 Assert::IsTrue(abs("0.15436783738642197865293282790786071397" - shootingResult.Derivative) < 
			 1e-38);
		     Assert::IsTrue(abs("0.11753094321729447380730293828117240339" - shootingResult.Value) < 1e-38);
			 Assert::IsTrue(abs("1.0" - shootingResult.Argument) < 1e-38);
		}

		static bool CheckFunc(InitCondition<mpfr::mpreal>& ic)
		{
			return true;
		}
	
		static bool CheckFuncDouble(InitCondition<double>& ic)
		{
			return true;
		}

		TEST_METHOD(XCanonStraightDerivativeBisectionTestMethod)
		{
			 mpreal::set_default_prec(128);
			 const int l = 1;
             XCannon<mpreal> xc(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, "0.1", &Troesch<mpreal, l>::StepFunc, &CheckFunc, 1e-38);
			 BisectionComponent<mpreal> bc(xc);
			 Assert::IsTrue(bc.DerivativeBisection("0.0", "1.0", "0.0", "1.0", "0.0", "1.0"));
			 auto knots = xc.GetKnotVector();
			 Assert::IsTrue(abs(knots[knots.size() - 1].Value - "1.0") < 1e-37);
			 Assert::IsTrue(abs(knots[knots.size() - 1].Derivative - "1.3414695048710118466622871534936817101") < 1e-37);
		}

		TEST_METHOD(XCanonInverseDerivativeBisectionTestMethod)
		{
			 mpreal::set_default_prec(128);
			 const int l = 3;
			 mpreal h = "0.001";
			 XCannon<mpreal> xcs(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, h, &Troesch<mpreal, l>::StepFunc, &CheckFunc, 1e-40);
			 auto shootingResults = xcs.Shoot("0", "1", "0", "0.2");

			 //xcs.SaveToFile("d:\\StraightFunction.txt");

             XCannonInverse<mpreal> xc(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, h, &Troesch<mpreal, l>::StepFunc, &CheckFunc, 1e-38);
			 auto shootingResultInverse = xc.Shoot("0.0", shootingResults.Value, "0.0", "5.0");	

			 Assert::IsTrue(abs(shootingResultInverse.Value - "1.0") < 10*h*h);
			 Assert::IsTrue(abs(shootingResultInverse.Derivative - 1/shootingResults.Derivative) < 10*h*h);

			 //xc.SaveToFile("d:\\InverseFunction.txt");

             XCannonInverse<mpreal> xcb(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, -h, &Troesch<mpreal, l>::StepFunc, &CheckFunc, 1e-38);
			 auto shootingResultInverseBack = xcb.Shoot(shootingResults.Value, "0.0", shootingResults.Argument, 1/shootingResults.Derivative);	

			 Assert::IsTrue(abs(shootingResultInverseBack.Value) < 100*h*h);
			 Assert::IsTrue(abs(1/shootingResultInverseBack.Derivative - "0.2") < 100*h*h, Message("Derivatives are too different"));

			 //xcb.SaveToFile("d:\\InverseFunctionBack.txt");
		}

		TEST_METHOD(XCanonHybridBisectionTestMethod)
		{
			 mpreal::set_default_prec(64);
			 const int l = 10;
			 mpreal h = "0.01";
			 TroeschHybridCannon<mpreal> thc(&Troesch<mpreal, l>::Nonlin, &Troesch<mpreal, l>::dNonlin, h, 
				 &Troesch<mpreal, l>::StepFunc, &CheckFunc, 1e-15);

			 std::function<int(const InitCondition<mpreal>&)> evalFunc = 
				 [](const InitCondition<mpreal>& ic) { return sgn(ic.Value - ic.Argument); };
			 BisectionComponent<mpreal> bc(thc);
			 Assert::IsTrue(bc.DerivativeBisectionGen("0.0", "1.0", "0.0", "1.0", "0.0", "1.0", evalFunc));
			 auto knots = thc.GetKnotVector();
			 Assert::IsTrue(abs(( knots[knots.size() - 1].Derivative - "148.4062956832580325")/knots[knots.size() - 1].Derivative) < 10*h*h);
			 Assert::IsTrue(abs( knots[knots.size() - 1].Argument - "1.0") < thc.GetPrecision());
			 Assert::IsTrue(abs( knots[knots.size() - 1].Value - "1.0") < thc.GetPrecision());
			 thc.SaveToFile("f:\\HybridBisection.txt");
		}

		TEST_METHOD(XCanonHybridBisectionDoubleTestMethod)
		{
			 const int l = 10;
			 double h = 0.01;
			 TroeschHybridCannon<double> thc(&Troesch<double, l>::Nonlin, &Troesch<double, l>::dNonlin, h, 
			 	 &Troesch<double, l>::StepFunc, &CheckFuncDouble, 1e-15);

			 std::function<int(const InitCondition<double>&)> evalFunc = 
				 [](const InitCondition<double>& ic) { return sgn(ic.Value - ic.Argument); };
			 BisectionComponent<double> bc(thc);
			 Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, evalFunc));
			 auto knots = thc.GetKnotVector();
			 Assert::IsTrue(abs(( knots[knots.size() - 1].Derivative - 148.4062956832580325)/knots[knots.size() - 1].Derivative) < 10*h*h);
			 Assert::IsTrue(abs( knots[knots.size() - 1].Argument - 1.0) < thc.GetPrecision());
			 Assert::IsTrue(abs( knots[knots.size() - 1].Value - 1.0) < thc.GetPrecision());
			 thc.SaveToFile("f:\\HybridBisectionDouble.txt");			
		}
	};
}