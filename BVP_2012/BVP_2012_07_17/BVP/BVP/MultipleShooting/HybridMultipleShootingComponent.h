#ifndef GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT
#define GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT

#include "..\Problems\ProblemAbstract.h"
#include "..\FunctionApproximation\GradientVector.h"
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\PointSimple.h"
#include "..\Cannon\TroeschHybridCannon.h"
#include "..\ShootingSimple\BisectionComponent.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
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
		Eigen::Matrix<T, Eigen::Dynamic,1> F((int)dim);

		auto A = _problem->GetACoeff();
		auto B = _problem->GetBCoeff();
		auto A_grad = _problem->GetACoeffGradient();
		auto B_grad = _problem->GetBCoeffGradient();

		auto AI = _problem->GetACoeffInverse();
		auto BI = _problem->GetBCoeffInverse();
		auto AI_grad = _problem->GetACoeffInverseGradient();
		auto BI_grad = _problem->GetBCoeffInverseGradient();

		InitCondition<T> curKnot;
		InitCondition<T> prevKnot;
		InitCondition<T> nextKnot;

		array<T, 3> A_grad_value;
		array<T, 3> B_grad_value;

		curKnot = meshData[0];
		nextKnot = meshData[1];
		typename X_Func_Gradient<T> dX;
		if (abs(curKnot.Derivative) <= 1.0)
		{
			dX = X_Func_Gradient<T>::X3_Func_Gradient(
				A(curKnot.Derivative, curKnot.Value, curKnot.Argument),
				B(curKnot.Derivative, curKnot.Value, curKnot.Argument), 
				curKnot.Derivative, 
				curKnot.Value, 
				nextKnot.Argument - curKnot.Argument, _precision);
			
			A_grad_value = A_grad(curKnot.Derivative, curKnot.Value, curKnot.Argument);
			B_grad_value = B_grad(curKnot.Derivative, curKnot.Value, curKnot.Argument);

			F(0) = dX.X - nextKnot.Value;

		} else
		{
			dX = X_Func_Gradient<T>::XI_Func_Gradient(
				AI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value),
				BI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value), 
				1/curKnot.Derivative, 
				curKnot.Argument, 
				nextKnot.Value - curKnot.Value, _precision);

			A_grad_value = AI_grad(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);
			B_grad_value = BI_grad(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);

			F(0) = dX.X - nextKnot.Argument;
		}

		JM.coeffRef(0,0) = dX.dA * A_grad_value[0] +
				            dX.dB * B_grad_value[0] +
				            dX.dC;

		JM.coeffRef(1,0) = dX.dhdA * A_grad_value[0] +
				            dX.dhdB * B_grad_value[0] +
				            dX.dhdC;

		if ((abs(nextKnot.Derivative) <= 1.0) ^ (abs(curKnot.Derivative) <= 1.0))
			F(1) = dX.dh - 1/nextKnot.Derivative;
		else
			F(1) = dX.dh - nextKnot.Derivative;

		for(std::vector<InitCondition<T>>::size_type i = 1; i < meshData.size() - 1; i++)
		{
			curKnot = meshData[i];
			prevKnot = meshData[i - 1];
			nextKnot = meshData[i + 1];

			int currRow = 2 * (int)i;
			int currCol = 2 * (int)i - 1;

			if (abs(curKnot.Derivative) <= 1.0)
			{
				dX = X_Func_Gradient<T>::X3_Func_Gradient(
					A(curKnot.Derivative, curKnot.Value, curKnot.Argument),
					B(curKnot.Derivative, curKnot.Value, curKnot.Argument), 
					curKnot.Derivative, 
					curKnot.Value, 
					nextKnot.Argument - curKnot.Argument, _precision);

				A_grad_value = A_grad(curKnot.Derivative, curKnot.Value, curKnot.Argument);
				B_grad_value = B_grad(curKnot.Derivative, curKnot.Value, curKnot.Argument);

			} else
			{
				dX = X_Func_Gradient<T>::XI_Func_Gradient(
					AI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value),
					BI(1/curKnot.Derivative, curKnot.Argument, curKnot.Value), 
					1/curKnot.Derivative, 
					curKnot.Argument, 
					nextKnot.Value - curKnot.Value, _precision);

				A_grad_value = AI_grad(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);
				B_grad_value = BI_grad(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);

			}

			JM.coeffRef(currRow, currCol)         = dX.dA * A_grad_value[0] +
						                            dX.dB * B_grad_value[0] +
						                            dX.dC;

			JM.coeffRef(currRow + 1, currCol)     = dX.dhdA * A_grad_value[0] +
						                            dX.dhdB * B_grad_value[0] + 
													dX.dhdC;

			JM.coeffRef(currRow - 2, currCol + 1) = 1.0;

			if ((abs(prevKnot.Derivative) <= 1.0) ^ (abs(curKnot.Derivative) <= 1.0))
			{
				JM.coeffRef(currRow - 1, currCol)     = - 1.0 / auxutils::sqr(curKnot.Derivative);

				JM.coeffRef(currRow, currCol + 1)     = dX.dA * A_grad_value[2] +
						                                dX.dB * B_grad_value[2] - 
														dX.dh;

				JM.coeffRef(currRow + 1, currCol + 1) = dX.dhdA * A_grad_value[2] +
						                                dX.dhdB * B_grad_value[2] - 
														dX.dhdh;
			}
			else
			{
				JM.coeffRef(currRow - 1, currCol)     = 1.0;

				JM.coeffRef(currRow, currCol + 1)     = dX.dA * A_grad_value[1] +
						                                dX.dB * B_grad_value[1] + 
						                                dX.dD;

				JM.coeffRef(currRow + 1, currCol + 1) = dX.dhdA * A_grad_value[1] +
						                                dX.dhdB * B_grad_value[1] + 
														dX.dhdD;
			}
		}
		
		InitCondition<T> penultKnot = meshData[meshData.size() - 2];
		InitCondition<T> lastKnot = meshData[meshData.size() - 1];

		if ((abs(penultKnot.Derivative) <= 1) ^ (abs(lastKnot.Derivative) <= 1))
			JM.coeffRef((int)dim - 1, (int)dim - 1) = - 1.0 / auxutils::sqr(lastKnot.Derivative);
		else
			JM.coeffRef((int)dim - 1, (int)dim - 1) = 1.0;

		//auxutils::SaveToFile(JM, "f:\\matr.txt");
		//auxutils::SaveToFile(F, "f:\\F.txt");

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