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
	typedef Eigen::Matrix<T, Eigen::Dynamic,1> Vector;
	typedef Eigen::SparseMatrix<T> Matrix;

	//Method that converts mesh data to vector
	Vector MeshDataToVector(vector<InitCondition<T>> meshData)
	{
		auto dim = 2 * (meshData.size() - 1);
		Vector result(dim);
		InitCondition<T> curKnot = meshData[0];
		InitCondition<T> prevKnot;

		result(0) = abs(curKnot.Derivative) <= 1.0 ? meshData[0].Derivative : 
			1/meshData[0].Derivative;

		for (std::vector<InitCondition<T>>::size_type i = 1; i < meshData.size() - 1; i++)
		{
			prevKnot = curKnot;
			curKnot = meshData[i];

			int curRow = 2 * (int)i - 1;

			result(curRow) = abs(curKnot.Derivative) <= 1.0 ? 
				curKnot.Derivative : 1/curKnot.Derivative;

			result(curRow + 1) = abs(prevKnot.Derivative) <= 1.0 ? 
				curKnot.Value : curKnot.Argument;
		}

		curKnot = meshData[meshData.size() - 1];

		result(dim - 1) = abs(curKnot.Derivative) <= 1.0 ? 
				curKnot.Derivative : 1/curKnot.Derivative;

		return result;
	}

	///Returns mesh data obtained from the previous mesh data and the given vector
	vector<InitCondition<T>> VectorToMeshData(vector<InitCondition<T>> baseMeshData, Vector vect)
	{
		vector<InitCondition<T>> result(std::begin(baseMeshData), std::end(baseMeshData));

		auto dim = baseMeshData.size();
		auto vectDim = 2 * (baseMeshData.size() - 1);

		if (vect.rows() != vectDim)
			throw exception("Nuexpected vector size");

        result[0].Derivative = abs(baseMeshData[0].Derivative) <= 1.0 ? vect(0) : 1/vect(0);

		for (std::vector<InitCondition<T>>::size_type i = 1; i < result.size() - 1; i++)
		{
			int curRow = 2 * (int)i - 1;

            result[i].Derivative = abs(baseMeshData[i].Derivative) <= 1.0 ? 
				vect(curRow) : 1/vect(curRow);

			if (abs(baseMeshData[i - 1].Derivative) <= 1.0)
				result[i].Value = vect(curRow + 1);
			else
				result[i].Argument = vect(curRow + 1);
		}

            result[dim - 1].Derivative = abs(baseMeshData[dim - 1].Derivative) <= 1.0 ? 
				vect(vectDim - 1) : 1/vect(vectDim - 1);

		return result;
	}

	///Method to compare two mesh data
	bool MeshDatasAreEqual(vector<InitCondition<T>> meshData1, vector<InitCondition<T>> meshData2)
	{
		if (meshData1.size() != meshData2.size())
			return false;

		for (std::vector<InitCondition<T>>::size_type i = 0; i < meshData1.size(); i++)
		{
			InitCondition<T> ic1 = meshData1[i];
			InitCondition<T> ic2 = meshData2[i];

			if (ic1.Derivative != ic2.Derivative || 
				ic1.Value != ic2.Value || 
				ic1.Argument != ic2.Argument)
				return false;
		}

		return true;
	}

	///A struct that contains data for a Newton method'd iteration 
	struct NewtonData
	{
	public:
		Matrix matrix;
		Vector F;
		Vector U;
	};

	///Generates Jacobi matrix from the given meth data vector
	void GenerateNewtonData(vector<InitCondition<T>> meshData, NewtonData& result)
	{
		auto dim = 2 * (meshData.size()-1);
		Matrix JM((int)dim, (int)dim);
		Vector F((int)dim);

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
			F(1) = dX.dh - nextKnot.Derivative;
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
			F(1) = dX.dh - 1/nextKnot.Derivative;
		}

		JM.coeffRef(0,0) = dX.dA * A_grad_value[0] +
				            dX.dB * B_grad_value[0] +
				            dX.dC;

		JM.coeffRef(1,0) = dX.dhdA * A_grad_value[0] +
				            dX.dhdB * B_grad_value[0] +
				            dX.dhdC;

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

				F(currRow) = dX.X - nextKnot.Value;
				F(currRow + 1) = dX.dh - nextKnot.Derivative;
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

				F(currRow) = dX.X - nextKnot.Argument;
				F(currRow + 1) = dX.dh - 1/nextKnot.Derivative;
			}

			JM.coeffRef(currRow, currCol)         = dX.dA * A_grad_value[0] +
						                            dX.dB * B_grad_value[0] +
						                            dX.dC;

			JM.coeffRef(currRow + 1, currCol)     = dX.dhdA * A_grad_value[0] +
						                            dX.dhdB * B_grad_value[0] + 
													dX.dhdC;

			JM.coeffRef(currRow - 2, currCol + 1) = - 1;

			if ((abs(prevKnot.Derivative) <= 1.0) ^ (abs(curKnot.Derivative) <= 1.0))
			{
				auto derivative = abs(curKnot.Derivative) <= 1.0 ? curKnot.Derivative : 1/curKnot.Derivative;
				JM.coeffRef(currRow - 1, currCol)     = 1.0 / auxutils::sqr(derivative);

				JM.coeffRef(currRow, currCol + 1)     = dX.dA * A_grad_value[2] +
						                                dX.dB * B_grad_value[2] - 
														dX.dh;

				JM.coeffRef(currRow + 1, currCol + 1) = dX.dhdA * A_grad_value[2] +
						                                dX.dhdB * B_grad_value[2] - 
														dX.dhdh;
			}
			else
			{
				JM.coeffRef(currRow - 1, currCol)     = - 1;

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

		if ((abs(penultKnot.Derivative) <= 1.0) ^ (abs(lastKnot.Derivative) <= 1.0))
		{
			auto derivative = abs(lastKnot.Derivative) <= 1.0 ? lastKnot.Derivative : 1/lastKnot.Derivative;
			JM.coeffRef((int)dim - 1, (int)dim - 1) = 1 / auxutils::sqr(derivative); 
		}
		else
			JM.coeffRef((int)dim - 1, (int)dim - 1) = - 1;


		result.matrix = JM;
		result.F = F;
		result.U = MeshDataToVector(meshData);

		//auxutils::SaveToFile(result.matrix, "F:\\matr.txt");
		//auxutils::SaveToFile(result.F, "F:\\F.txt");
		//auxutils::SaveToFile(result.U, "F:\\V.txt");
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

		T h = 0.01;
		TroeschHybridCannon<T> thc((*_problem), h, 0.001);

		std::function<int(const InitCondition<T>&)> evalFunc = 
			[](const InitCondition<T>& ic) { return sgn(ic.Value - ic.Argument); };
		BisectionComponent<T> bc(thc);
		bc.DerivativeBisectionGen(0, 1, 0, 1, 0, 1, evalFunc);
		_meshData = thc.GetKnotVector();

		//auxutils::SaveToFile(_meshData, "f:\\ExactSolution.txt");

		_meshData[_meshData.size() - 1].Value = 1;;
		_meshData[_meshData.size() - 1].Argument = 1;;

		NewtonData ND;
		vector<InitCondition<T>> MD(std::begin(_meshData), std::end(_meshData));

		Vector newVect;
		T norm = 1;
		for (int i = 0 ; i < 10 && norm > _precision; i++)
		{
			GenerateNewtonData(MD, ND);

			Eigen::BiCGSTAB<Matrix> solver; // Our BiCGStab solver;
			solver.setMaxIterations(100000);

			solver.compute(ND.matrix);// Initialization of the solver
			Vector correction = solver.solve(ND.F);

			T norm1;
			if (newVect.rows() == ND.U.rows())
			    norm1 = (newVect - ND.U).norm();

			newVect = ND.U - correction;

			//auxutils::SaveToFile(newVect, "f:\\Sol.txt");
			//auxutils::SaveToFile(correction, "f:\\correction.txt");

		    MD = VectorToMeshData(MD, newVect);

			norm = correction.norm();
		}

		//auxutils::SaveToFile(MD, "f:\\NewtonExactSolution.txt");

		//bool sanityCheck = MeshDatasAreEqual(_meshData, newMeshData);
	}
};

#endif