#include "CppUnitTest.h"

#include "..\BVP\Cannon\XCannon.h"
#include "..\BVP\Problems\TroeschProblem.h"
#include "..\BVP\ShootingSimple\BisectionComponent.h"
#include "..\BVP\Cannon\HybridCannon.h"
#include "UnitTestAux.h"
#include "..\BVP\Utils\AuxUtils.h"

using namespace UnitTestAux;

typedef float_50_noet numTypeMp;

typedef number<cpp_dec_float<40>, et_off> multiprec_float40;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GeneralTest
{
	TEST_CLASS(XCannonTestUnit)
	{
	public:
		
		TEST_METHOD(XCanonStraightTestMethod)
		{
			 const int l = 1;
			 TroeschProblem<multiprec_float40> tpf(l);
			 XCannon<multiprec_float40> xc(tpf, multiprec_float40("0.01"), 1e-40);

			 auto shootingResult = xc.Shoot(multiprec_float40("0"), multiprec_float40("1"), multiprec_float40("0"), multiprec_float40("0.1"));	

			 auto derivativeDifference = abs(multiprec_float40("0.15436783738642197865293282790786071397") - shootingResult.Derivative);
		
			 Assert::IsTrue(derivativeDifference < 
				 1e-38, Message("Derivative is different" + auxutils::ToString(derivativeDifference)));
		     Assert::IsTrue(abs(multiprec_float40("0.11753094321729447380730293828117240339") - shootingResult.Value) < 1e-38);
			 Assert::IsTrue(abs(multiprec_float40("1.0") - shootingResult.Argument) < 1e-38);
		}

		TEST_METHOD(XCanonStraightDerivativeBisectionTestMethod)
		{
			 const int l = 1;
			 TroeschProblem<multiprec_float40> tpf(l);
			 XCannon<multiprec_float40> xc(tpf, multiprec_float40("0.1"), 1e-38);
			 BisectionComponent<multiprec_float40> bc(xc);
			 Assert::IsTrue(bc.DerivativeBisection(multiprec_float40("0.0"), multiprec_float40("1.0"), multiprec_float40("0.0"), 
				 multiprec_float40("1.0"), multiprec_float40("0.0"), multiprec_float40("1.0")));
			 auto knots = xc.GetKnotVector();
			 Assert::IsTrue(abs(knots[knots.size() - 1].Value - multiprec_float40("1.0")) < 1e-37);
			 Assert::IsTrue(abs(knots[knots.size() - 1].Derivative - multiprec_float40("1.3414695048710118466622871534936817101")) < 1e-37);
		}

		TEST_METHOD(XCanonInverseDerivativeBisectionTestMethod)
		{
			 const int l = 3;
			 multiprec_float40 h(0.001);
			 TroeschProblem<multiprec_float40> tpf(l);
			 XCannon<multiprec_float40> xcs(tpf, h, 1e-40);
			 auto shootingResults = xcs.Shoot(multiprec_float40("0"), multiprec_float40("1"), multiprec_float40("0"), multiprec_float40("0.2"));

			 XCannonInverse<multiprec_float40> xc(tpf, h, 1e-38);
			 auto shootingResultInverse = xc.Shoot(multiprec_float40("0.0"), shootingResults.Value, multiprec_float40("0.0"), multiprec_float40("5.0"));	

			 Assert::IsTrue(abs(shootingResultInverse.Value - multiprec_float40("1.0")) < 10*h*h);
			 Assert::IsTrue(abs(shootingResultInverse.Derivative - 1/shootingResults.Derivative) < 10*h*h);

			 XCannonInverse<multiprec_float40> xcb(tpf, -h, 1e-38);
			 auto shootingResultInverseBack = xcb.Shoot(shootingResults.Value, multiprec_float40("0.0"), 
				 shootingResults.Argument, 1/shootingResults.Derivative);	

			 Assert::IsTrue(abs(shootingResultInverseBack.Value) < 100*h*h);
			 Assert::IsTrue(abs(1/shootingResultInverseBack.Derivative - multiprec_float40("0.2")) < 100*h*h, Message("Derivatives are too different"));
		}

		TEST_METHOD(XCanonHybridBisectionTestMethod)
		{
			 const int l = 10;
			 TroeschProblem<numTypeMp> tpf(l);
			 numTypeMp h = (numTypeMp)"0.01";
			 HybridCannon<numTypeMp> thc(tpf, h, 1e-29);

			 std::function<int(const InitCondition<numTypeMp>&)> evalFunc = 
				 [](const InitCondition<numTypeMp>& ic) { return sgn(ic.Value - ic.Argument); };
			 BisectionComponent<numTypeMp> bc(thc);
			 Assert::IsTrue(bc.DerivativeBisectionGen((numTypeMp)"0.0", (numTypeMp)"1.0", (numTypeMp)"0.0", (numTypeMp)"1.0", (numTypeMp)"0.0", (numTypeMp)"1.0", evalFunc));
			 auto knots = thc.GetKnotVector();

			 Assert::IsTrue(abs( knots[knots.size() - 1].Derivative - (numTypeMp)"148.406295683258037233643960143") < knots[knots.size() - 1].Derivative*thc.GetPrecision());
			 Assert::IsTrue(abs( knots[knots.size() - 1].Argument - (numTypeMp)"1.0") < thc.GetPrecision());
			 Assert::IsTrue(abs( knots[knots.size() - 1].Value - (numTypeMp)"1.0") < thc.GetPrecision());
		}

		TEST_METHOD(XCanonHybridBisectionDoubleTestMethod)
		{
			 const int l = 10;
			 double h = 0.01;
			 TroeschProblem<double> tpf(l);
			 HybridCannon<double> thc(tpf, h, 1e-15);

			 std::function<int(const InitCondition<double>&)> evalFunc = 
				 [](const InitCondition<double>& ic) { return sgn(ic.Value - ic.Argument); };
			 BisectionComponent<double> bc(thc);
			 Assert::IsTrue(bc.DerivativeBisectionGen(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, evalFunc));
			 auto knots = thc.GetKnotVector();
			 Assert::IsTrue(abs(( knots[knots.size() - 1].Derivative - 148.4062956832580325)/knots[knots.size() - 1].Derivative) < 10*h*h);
			 Assert::IsTrue(abs( knots[knots.size() - 1].Argument - 1.0) < thc.GetPrecision());
			 Assert::IsTrue(abs( knots[knots.size() - 1].Value - 1.0) < thc.GetPrecision());
		}
	};
}