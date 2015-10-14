#ifndef GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT
#define GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT

#include "..\Problems\ProblemAbstract.h"
#include "..\FunctionApproximation\GradientVector.h"
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\PointSimple.h"
#include "..\Cannon\TroeschHybridCannon.h"
#include "..\ShootingSimple\BisectionComponent.h"

#include <Eigen/Sparse>
#include <Eigen/MPRealSupport>

#include "..\Utils\Exceptions.h"
#include "..\FunctionApproximation\X_Function.h"
#include <vector>

/// A component that implements hybrid multiple shooting algorithm
template <class T>
class HybridMultipleShootingComponent
{
private: 
	ProblemAbstract<T>* _problem;
	vector<InitCondition<T>> _meshData;
	PointSimple<T> _ptLeft;
	PointSimple<T> _ptRight;
	double _precision;

	///Generates Jacobi matrix from the given meth data vector
	Eigen::SparseMatrix<T> GenerateJacobiMatrix(vector<InitCondition<T>> meshData)
	{
		auto dim = 2 * (meshData.size()-1);
		Eigen::SparseMatrix<T> JM((int)dim, (int)dim);


		auto A = _problem->GetACoeff();
		auto B = _problem->GetBCoeff();
		auto A_grad = _problem->GetACoeffGradient();
		auto B_grad = _problem->GetBCoeffGradient();

		auto AI = _problem->GetACoeffInverse();
		auto BI = _problem->GetBCoeffInverse();
		auto AI_grad = _problem->GetACoeffInverseGradient();
		auto BI_grad = _problem->GetBCoeffInverseGradient();

		InitCondition<T> firstKnot = meshData[0];

		typename X_Func_Gradient<T> dX;
		if (abs(firstKnot.Derivative) <= 1.0)
		{
			dX = X_Func_Gradient<T>::X3_Func_Gradient(
				A(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument),
				B(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument), 
				firstKnot.Derivative, 
				firstKnot.Value, 
				meshData[1].Argument - firstKnot.Argument, _precision);
			
			auto A_grad_value = A_grad(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument);
			auto B_gard_value = B_grad(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument);

			JM.coeffRef(0,0) = dX.dA * A_grad_value[0] +
				               dX.dB * B_gard_value[0] +
				               dX.dC;

			JM.coeffRef(1,0) = dX.dhdA * A_grad_value[0] +
							   dX.dhdB * B_gard_value[0] +
				               dX.dhdC;
		} else
		{
			dX = X_Func_Gradient<T>::XI_Func_Gradient(
				AI(1/firstKnot.Derivative, firstKnot.Argument, firstKnot.Value),
				BI(1/firstKnot.Derivative, firstKnot.Argument, firstKnot.Value), 
				1/firstKnot.Derivative, 
				firstKnot.Argument, 
				meshData[1].Value - firstKnot.Value, _precision);

			auto AI_grad_value = AI_grad(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument);
			auto BI_grad_value = BI_grad(firstKnot.Derivative, firstKnot.Value, firstKnot.Argument);

			JM.coeffRef(0,0) = dX.dA * AI_grad_value[0] +
				               dX.dB * BI_grad_value[0] +
				               dX.dC;

			JM.coeffRef(1,0) = dX.dhdA * AI_grad_value[0] +
				               dX.dhdB * BI_grad_value[0] +
				               dX.dhdC;
		}

		//auxutils::SaveToFile(JM, "f:\\matr.txt");

		for(std::vector<InitCondition<T>>::size_type i = 1; i != meshData.size() - 1; i++)
		{
			InitCondition<T> curKnot = meshData[i];

			if (abs(curKnot.Derivative) <= 1.0)
			{
				dX = X_Func_Gradient<T>::X3_Func_Gradient(
					A(curKnot.Derivative, curKnot.Value, curKnot.Argument),
					B(curKnot.Derivative, curKnot.Value, curKnot.Argument), 
					curKnot.Derivative, 
					curKnot.Value, 
					meshData[i + 1].Argument - curKnot.Argument, _precision);
			} else
			{
				dX = X_Func_Gradient<T>::XI_Func_Gradient(
					AI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value),
					BI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value), 
					1/curKnot.Derivative, 
					curKnot.Argument, 
					meshData[i + 1].Value - curKnot.Value, _precision);
			}
		}
		
		return JM;
	}
public:
	///Constructor
	HybridMultipleShootingComponent(ProblemAbstract<T>& problem, 
		PointSimple<T> ptLeft, PointSimple<T> ptRight, double precision)
	{
		//throw NotImplementedException();
		_problem = &problem;
		_ptLeft = ptLeft;
		_ptRight = ptRight;
		_precision = precision;

		T h = 0.1;
		TroeschHybridCannon<T> thc((*_problem), h, precision);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[](const InitCondition<T>& ic) { return sgn(ic.Value - ic.Argument); };
		BisectionComponent<T> bc(thc);
		bc.DerivativeBisectionGen(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, evalFunc);
		_meshData = thc.GetKnotVector();

		GenerateJacobiMatrix(_meshData);
	}
};

#endif