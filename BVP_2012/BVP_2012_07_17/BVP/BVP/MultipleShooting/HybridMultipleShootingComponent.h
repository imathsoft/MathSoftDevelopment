#ifndef GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT
#define GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT

#include "..\Problems\ProblemAbstract.h"
#include "..\FunctionApproximation\GradientVector.h"
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\PointSimple.h"
#include "..\Cannon\HybridCannon.h"
#include "..\ShootingSimple\BisectionComponent.h"

#include "..\Utils\Exceptions.h"
#include "..\FunctionApproximation\X_Function.h"
#include <vector>
#include "FourDiagonalSweepMethodStruct.h"

/// A component that implements hybrid multiple shooting algorithm
template <class T>
class HybridMultipleShootingComponent
{
private: 
	const std::unique_ptr<ProblemAbstract<T>> _problem;
	T _precision;

	const T _critical_d_value;

	///Method to refine mesh according to the given step
	vector<InitCondition<T>> RefineMesh(const vector<InitCondition<T>>& meshData, T step)
	{
		vector<InitCondition<T>> result; 

		if (meshData.size() == 0 || step <= 0)
			return result;

		InitCondition<T> nextKnot = meshData[0];
		result.push_back(nextKnot);
		InitCondition<T> prevKnot;

		T maxSquaredDist = 2 * step * step;

		for (std::vector<InitCondition<T>>::size_type i = 1; i < meshData.size(); i++)
		{
			prevKnot = nextKnot;
			nextKnot = meshData[i];

			InitCondition<T> vect = nextKnot - prevKnot;

			if (vect.NormSquaredPartial<false>() >= maxSquaredDist)
			{
				T naiveNorm = auxutils::Sqrt(vect.NormSquaredPartial<false>()); //TODO : avoid using the square root function 
				vect = vect / naiveNorm;
				T tau = step;

				while (tau < naiveNorm)
				{
					InitCondition<T> newKnot = prevKnot + vect * tau;
					result.push_back(newKnot);
					tau += step;
				}
			}
        	result.push_back(nextKnot);
		}

		return result;
	}

	//Method that converts mesh data to vector
	vector<T> MeshDataToStdVector(const vector<InitCondition<T>>& meshData)
	{
		auto dim = 2 * (meshData.size() - 1);
		vector<T> result;
		result.resize(dim);
		InitCondition<T> curKnot = meshData[0];
		InitCondition<T> prevKnot;

		result[0] = abs(curKnot.Derivative) <= _critical_d_value ? meshData[0].Derivative : 
			1/meshData[0].Derivative;

		for (std::vector<InitCondition<T>>::size_type i = 1; i < meshData.size() - 1; i++)
		{
			prevKnot = curKnot;
			curKnot = meshData[i];

			int curRow = 2 * (int)i - 1;

			result[curRow + 1] = abs(curKnot.Derivative) <= _critical_d_value ? 
				curKnot.Derivative : 1/curKnot.Derivative;

			result[curRow] = abs(prevKnot.Derivative) <= _critical_d_value ? 
				curKnot.Value : curKnot.Argument;
		}

		curKnot = meshData[meshData.size() - 1];

		result[dim - 1] = abs(curKnot.Derivative) <= _critical_d_value ? 
				curKnot.Derivative : 1/curKnot.Derivative;

		return result;
	}

	/*
	//Method that converts mesh data to vector
	Vector MeshDataToVector(vector<InitCondition<T>> meshData)
	{
		auto dim = 2 * (meshData.size() - 1);
		Vector result(dim);

		auto stdVector = MeshDataToStdVector(meshData);

		for (int i = 0; i < dim; i++)
			result(i) = stdVector[i];

		return result;
	}
	*/

	///Returns mesh data obtained from the previous mesh data and the given vector
	vector<InitCondition<T>> VectorToMeshData(const vector<InitCondition<T>>& baseMeshData, const vector<T>& vect)
	{
		vector<InitCondition<T>> result(std::begin(baseMeshData), std::end(baseMeshData));

		auto dim = baseMeshData.size();
		auto vectDim = 2 * (baseMeshData.size() - 1);

		if (vect.size() != vectDim)
			throw exception("Nuexpected vector size");

        result[0].Derivative = abs(baseMeshData[0].Derivative) <= _critical_d_value ? vect[0] : 1/vect[0];

		for (std::vector<InitCondition<T>>::size_type i = 1; i < result.size() - 1; i++)
		{
			int curRow = 2 * (int)i - 1;

            result[i].Derivative = abs(baseMeshData[i].Derivative) <= _critical_d_value ? 
				vect[curRow + 1] : 1/vect[curRow + 1];

			if (abs(baseMeshData[i - 1].Derivative) <= _critical_d_value)
				result[i].Value = vect[curRow];
			else
				result[i].Argument = vect[curRow];
		}

            result[dim - 1].Derivative = abs(baseMeshData[dim - 1].Derivative) <= _critical_d_value ? 
				vect[vectDim - 1] : 1/vect[vectDim - 1];

		return result;
	}

	/*
	///Returns mesh data obtained from the previous mesh data and the given vector
	vector<InitCondition<T>> VectorToMeshData(vector<InitCondition<T>> baseMeshData, Vector vect)
	{
		vector<T> vec;
		vec.resize(vect.rows());
		for (int i = 0; i < vect.rows(); i++)
			vec[i] = vect[i];

		return VectorToMeshData(baseMeshData, vec);
	}
	*/

	///Method to compare two mesh data
	T DiffMeshDatas(const vector<InitCondition<T>>& meshData1, const vector<InitCondition<T>>& meshData2)
	{
		if (meshData1.size() != meshData2.size())
			return 100;

		T diff = 0;

		for (std::vector<InitCondition<T>>::size_type i = 0; i < meshData1.size(); i++)
		{
			InitCondition<T> ic1 = meshData1[i];
			InitCondition<T> ic2 = meshData2[i];

			diff = max(diff, abs(ic1.Derivative - ic2.Derivative));
			diff = max(diff, abs(ic1.Value - ic2.Value));
			diff = max(diff, abs(ic1.Argument - ic2.Argument));
			diff = max(diff, abs(ic1.SecDerivative - ic2.SecDerivative));
		}

		return diff;
	}

	///Method to run sweep method
	vector<T> RunSweepMethod(FourDiagonalSweepMethodStruct<T>& fDSMS, T& absCorrection)
	{
		vector<T> result;
		result.resize(fDSMS.RowCount());

		for (int i = 0; i < fDSMS.RowCount(); i++)
		{
			T denominator = fDSMS[i][2];
			if (denominator != 0)
			{
				fDSMS[i][0] /= denominator;
				fDSMS[i][1] /= denominator;
				fDSMS[i][2] = 1;
				fDSMS[i][3] /= denominator;
			}
		}

		int dim = fDSMS.RowCount() / 2;
		for (int  i = 1; i < dim; i++)
		{
			int curRow = 2*i;
			for (int j = curRow; j < curRow + 2; j++)
			{
				fDSMS[j][3] -= fDSMS[j][0]*fDSMS[curRow - 2][3] + fDSMS[j][1]*fDSMS[curRow - 1][3];
				fDSMS[j][0] = -(fDSMS[j][0]*fDSMS[curRow - 2][0] + fDSMS[j][1]*fDSMS[curRow - 1][0]);
			}
		}

		int preLastRow = fDSMS.RowCount() - 2;

		T det = fDSMS[preLastRow][0]*fDSMS[preLastRow + 1][2] - fDSMS[preLastRow][2]*fDSMS[preLastRow + 1][0];

		if (det == 0)
			throw exception("Invalid determinant");

		T C0 = (fDSMS[preLastRow][3]*fDSMS[preLastRow + 1][2] - fDSMS[preLastRow][2]*fDSMS[preLastRow + 1][3])/det;
		result[result.size() - 1] = 
			(fDSMS[preLastRow][0]*fDSMS[preLastRow + 1][3] - fDSMS[preLastRow][3]*fDSMS[preLastRow + 1][0])/det;

		result[0] = C0;

		absCorrection = max(abs(C0), abs(result[result.size() - 1]));

		for (int i = 0; i < fDSMS.RowCount() - 2; i++)
		{
			result[i + 1] =  - fDSMS[i][0]*C0 + fDSMS[i][3];
			absCorrection = max(absCorrection, abs(result[i + 1]));
		}

		for (int i = 0; i < fDSMS.RowCount(); i++)
			result[i] = fDSMS[i][4] - result[i];

		return result;
	}

	///Calculates and returns the "acheivable precision"
	T GetAcheivablePrecision(size_t numberOfKnots)
	{
		return 10 * max(auxutils::RoughSqrt((T)numberOfKnots) * std::numeric_limits<T>::epsilon(), _precision);
	}

	///Method to run Newton iterations
	vector<InitCondition<T>> RunNewtonIterations(const vector<InitCondition<T>>& meshData, T desiredMaxStepSize, bool& succeeded )
	{
		succeeded = false;
		_debugData.clear();
		vector<InitCondition<T>> MD(std::begin(meshData), std::end(meshData));

		T absCorrection = 1;
		T oldAbsCorrection;

		///DD-160228: Achievable precision is limited by round-off errors
		///which, in turn, depend on the number of knots
		auto achievablePrec = GetAcheivablePrecision(MD.size());

		for (int i = 0; i < MaxNumberOfNewtonIterations && 
			(absCorrection > achievablePrec || oldAbsCorrection >= 1); i++)
		{
			oldAbsCorrection = absCorrection;

			if (absCorrection*auxutils::RoughSqrt((T)MD.size()) < 1) // the process is starting to converge
			{
				MD = RefineMesh(MD, desiredMaxStepSize);

				achievablePrec = GetAcheivablePrecision(MD.size());
			}

			FourDiagonalSweepMethodStruct<T> fDSMS = GenerateSweepMethodStruct(MD);

			vector<T> newVector = RunSweepMethod(fDSMS, absCorrection);

			_debugData.push_back(std::make_pair((int)MD.size(), absCorrection));

			_correctionMagnitudes.push_back(absCorrection);

			if (absCorrection / oldAbsCorrection > 100 && absCorrection > 100)
				break; // the process is divergent

			MD = VectorToMeshData(MD, newVector);
		}

		succeeded = absCorrection <= achievablePrec;

		return MD;
	}

	///Method to generate sweep method struct
	FourDiagonalSweepMethodStruct<T> GenerateSweepMethodStruct(const vector<InitCondition<T>>& meshData)
	{
		auto dim = 2 * (meshData.size()-1);
		FourDiagonalSweepMethodStruct<T> result((int)dim);

		InitCondition<T> curKnot;
		InitCondition<T> prevKnot;
		InitCondition<T> nextKnot;

		array<T, 3> A_grad_value;
		array<T, 3> B_grad_value;

		curKnot = meshData[0];
		nextKnot = meshData[1];
		typename X_Func_Gradient<T> dX;
		if (abs(curKnot.Derivative) <= _critical_d_value)
		{
			dX = _problem->step_streight_gradient(
				_problem->GetACoeff(curKnot.Derivative, curKnot.Value, curKnot.Argument),
				_problem->GetBCoeff(curKnot.Derivative, curKnot.Value, curKnot.Argument), 
				curKnot.Derivative, 
				curKnot.Value, 
				_problem->GetECoeff(curKnot.Argument),
				_problem->GetFCoeff(curKnot.Argument),
				nextKnot.Argument - curKnot.Argument, _precision);
			
			A_grad_value = _problem->GetACoeffGradient(curKnot.Derivative, curKnot.Value, curKnot.Argument);
			B_grad_value = _problem->GetBCoeffGradient(curKnot.Derivative, curKnot.Value, curKnot.Argument);

			result[0][3] = dX.X - nextKnot.Value;
			result[1][3] = dX.dh - nextKnot.Derivative;
		} else
		{
			dX = _problem->step_inverse_gradient(
				_problem->GetACoeffInverse(1/curKnot.Derivative, curKnot.Argument, curKnot.Value),
				_problem->GetBCoeffInverse(1/curKnot.Derivative, curKnot.Argument, curKnot.Value), 
				1/curKnot.Derivative, 
				curKnot.Argument, 
				nextKnot.Value - curKnot.Value, _precision);

			A_grad_value = _problem->GetACoeffInverseGradient(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);
			B_grad_value = _problem->GetBCoeffInverseGradient(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);

			result[0][3] = dX.X - nextKnot.Argument;
			result[1][3] = dX.dh - 1/nextKnot.Derivative;
		}

		result[0][0] = dX.dA * A_grad_value[0] +
				            dX.dB * B_grad_value[0] +
				            dX.dC;

		result[1][0] = dX.dhdA * A_grad_value[0] +
				            dX.dhdB * B_grad_value[0] +
				            dX.dhdC;

		for(std::vector<InitCondition<T>>::size_type i = 1; i < meshData.size() - 1; i++)
		{
			curKnot = meshData[i];
			prevKnot = meshData[i - 1];
			nextKnot = meshData[i + 1];

			int currRow = 2 * (int)i;

			if (abs(curKnot.Derivative) <= _critical_d_value)
			{
				dX = _problem->step_streight_gradient(
					_problem->GetACoeff(curKnot.Derivative, curKnot.Value, curKnot.Argument),
					_problem->GetBCoeff(curKnot.Derivative, curKnot.Value, curKnot.Argument), 
					curKnot.Derivative, 
					curKnot.Value, 
					_problem->GetECoeff(curKnot.Argument),
					_problem->GetFCoeff(curKnot.Argument),
					nextKnot.Argument - curKnot.Argument, _precision);

				A_grad_value = _problem->GetACoeffGradient(curKnot.Derivative, curKnot.Value, curKnot.Argument);
				B_grad_value = _problem->GetBCoeffGradient(curKnot.Derivative, curKnot.Value, curKnot.Argument);

				result[currRow][3] = dX.X - nextKnot.Value;
				result[currRow + 1][3] = dX.dh - nextKnot.Derivative;
			} else
			{
				dX = _problem->step_inverse_gradient(
					_problem->GetACoeffInverse(1/curKnot.Derivative, curKnot.Argument, curKnot.Value),
					_problem->GetBCoeffInverse(1/curKnot.Derivative, curKnot.Argument, curKnot.Value), 
					1/curKnot.Derivative, 
					curKnot.Argument, 
					nextKnot.Value - curKnot.Value, _precision);

				A_grad_value = _problem->GetACoeffInverseGradient(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);
				B_grad_value = _problem->GetBCoeffInverseGradient(1/curKnot.Derivative, curKnot.Argument, curKnot.Value);

				result[currRow][3] = dX.X - nextKnot.Argument;
				result[currRow + 1][3] = dX.dh - 1/nextKnot.Derivative;
			}

			result[currRow][1]         = dX.dA * A_grad_value[0] +
						                            dX.dB * B_grad_value[0] +
						                            dX.dC;

			result[currRow + 1][1]     = dX.dhdA * A_grad_value[0] +
						                            dX.dhdB * B_grad_value[0] + 
													dX.dhdC;

			result[currRow - 2][2] = - 1;

			if ((abs(prevKnot.Derivative) <= _critical_d_value) ^ (abs(curKnot.Derivative) <= _critical_d_value))
			{
				auto derivative = abs(curKnot.Derivative) <= _critical_d_value ? curKnot.Derivative : 1/curKnot.Derivative;
				result[currRow - 1][2]     = 1.0 / auxutils::sqr(derivative);

				result[currRow][0]     = dX.dA * A_grad_value[2] +
						                                dX.dB * B_grad_value[2] - 
														dX.dh + 
														dX.dE*_problem->GetdEdX(curKnot.Argument) + 
														dX.dF * _problem->GetdFdX(curKnot.Argument);

				result[currRow + 1][0] = dX.dhdA * A_grad_value[2] +
						                                dX.dhdB * B_grad_value[2] - 
														dX.dhdh +
														dX.dhdE*_problem->GetdEdX(curKnot.Argument) + 
														dX.dhdF * _problem->GetdFdX(curKnot.Argument);

			}
			else
			{
				result[currRow - 1][2]     = - 1;

				result[currRow][0]     = dX.dA * A_grad_value[1] +
						                                dX.dB * B_grad_value[1] + 
						                                dX.dD;

				result[currRow + 1][0] = dX.dhdA * A_grad_value[1] +
						                                dX.dhdB * B_grad_value[1] + 
														dX.dhdD;
			}
		}
		
		InitCondition<T> penultKnot = meshData[meshData.size() - 2];
		InitCondition<T> lastKnot = meshData[meshData.size() - 1];

		if ((abs(penultKnot.Derivative) <= _critical_d_value) ^ (abs(lastKnot.Derivative) <= _critical_d_value))
		{
			auto derivative = abs(lastKnot.Derivative) <= _critical_d_value ? lastKnot.Derivative : 1/lastKnot.Derivative;
			result[(int)dim - 1][2] = 1 / auxutils::sqr(derivative); 
		}
		else
			result[(int)dim - 1][2] = - 1;

		auto uVector = MeshDataToStdVector(meshData);

		for (std::vector<T>::size_type i = 0; i < uVector.size(); i++)
			result[(int)i][4] = uVector[i];

		return result;
	};

	/*
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

			JM.coeffRef(currRow, currCol + 1)         = dX.dA * A_grad_value[0] +
						                            dX.dB * B_grad_value[0] +
						                            dX.dC;

			JM.coeffRef(currRow + 1, currCol + 1)     = dX.dhdA * A_grad_value[0] +
						                            dX.dhdB * B_grad_value[0] + 
													dX.dhdC;

			JM.coeffRef(currRow - 2, currCol) = - 1;

			if ((abs(prevKnot.Derivative) <= 1.0) ^ (abs(curKnot.Derivative) <= 1.0))
			{
				auto derivative = abs(curKnot.Derivative) <= 1.0 ? curKnot.Derivative : 1/curKnot.Derivative;
				JM.coeffRef(currRow - 1, currCol + 1)     = 1.0 / auxutils::sqr(derivative);

				JM.coeffRef(currRow, currCol)     = dX.dA * A_grad_value[2] +
						                                dX.dB * B_grad_value[2] - 
														dX.dh;

				JM.coeffRef(currRow + 1, currCol) = dX.dhdA * A_grad_value[2] +
						                                dX.dhdB * B_grad_value[2] - 
														dX.dhdh;
			}
			else
			{
				JM.coeffRef(currRow - 1, currCol + 1)     = - 1;

				JM.coeffRef(currRow, currCol)     = dX.dA * A_grad_value[1] +
						                                dX.dB * B_grad_value[1] + 
						                                dX.dD;

				JM.coeffRef(currRow + 1, currCol) = dX.dhdA * A_grad_value[1] +
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
	*/
	std::vector<std::pair<int, T>> _debugData;

	std::vector<T> _correctionMagnitudes;
public:
	int MaxNumberOfNewtonIterations;

	///Returns container with absolute values of approximation corrections added at each iteration of the Newton's method
	///It is supposed to be used in testing purposes to assert quadratic convergence rate
	std::vector<T> GetCorrectionMgnitudes()
	{
		return _correctionMagnitudes;
	}

	//Method to save debug data in the specified file (as text)
	void SaveDebugData(const char* filename)
	{
		std::ofstream file;
		file.precision(std::numeric_limits<T>::digits10);
		file.open (filename);
		for (std::vector<std::pair<int, T>>::const_iterator m = _debugData.begin(); m != _debugData.end(); ++m)
		{
			auto ic = (*m);
			file << ic.first << " " << ic.second << endl;
		}
        file.close();
	}

	///Constructor
	HybridMultipleShootingComponent(const ProblemAbstract<T>& problem, const T critical_d_value = T(1)) : 
		_problem(problem.copy()), 
		_critical_d_value(critical_d_value)
	{
		_precision = 10 * std::numeric_limits<T>::epsilon();
		MaxNumberOfNewtonIterations = 20;
	}

	vector<InitCondition<T>> Run(
		const PointSimple<T>& ptLeft, const PointSimple<T>& ptRight, 
		const T desiredStepSize,
		bool& succeeded,
		const T customPreStep = (T)0)
	{
		T h = desiredStepSize;
		_correctionMagnitudes.clear();

		T preStep = customPreStep > 0 ? customPreStep : min((T)1/100, 10*h);
		HybridCannon<T> thc((*HybridMultipleShootingComponent::_problem), preStep, min((T)1/100, h/10));
		std::function<int(const InitCondition<T>&)> evalFunc = 
			[](const InitCondition<T>& ic) { return sgn(ic.Value - ic.Argument); };
		BisectionComponent<T> bc(thc);
		bc.DerivativeBisectionGen(ptLeft, ptRight, 0, 1, evalFunc);
		auto meshData = thc.GetKnotVector();

		meshData[meshData.size() - 1].Value = ptRight.Value;
		meshData[meshData.size() - 1].Argument = ptRight.Argument;
		
		auto sol = RunNewtonIterations(meshData, h, succeeded);

		return sol;
	}

	vector<InitCondition<T>> Run(const vector<InitCondition<T>>& initialGuess, 
		const T desiredStepSize,
		bool& succeeded)
	{
		T h = desiredStepSize;
		_correctionMagnitudes.clear();

		auto sol = RunNewtonIterations(initialGuess, h, succeeded);

		return sol;
	}


};

#endif