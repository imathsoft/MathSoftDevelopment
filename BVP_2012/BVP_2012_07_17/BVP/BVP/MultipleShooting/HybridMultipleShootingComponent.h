#ifndef GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT
#define GUARD_HYBRID_MULTIPLE_SHOOTING_COMPONENT

#include "..\Problems\ProblemAbstract.h"
#include "..\FunctionApproximation\GradientVector.h"
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\FunctionApproximation\PointSimple.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/MPRealSupport>

#include "..\Utils\Exceptions.h"
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

	///Generates Jacobi matrix from the given meth data vector
	Eigen::SparseMatrix<T> GenerateJacobiMatrix(vector<InitCondition<T>> meshData)
	{
		auto dim = 2 * (meshData.size()-1);
		Eigen::SparseMatrix<T> JM((int)dim, (int)dim);

		return JM;
	}
public:
	///Constructor
	HybridMultipleShootingComponent(ProblemAbstract<T>& problem, 
		PointSimple<T> ptLeft, PointSimple<T> ptRight)
	{
		throw NotImplementedException();
		_problem = &problem;
		_ptLeft = ptLeft;
		_ptRight = ptRight;

		GenerateJacobiMatrix(_meshData);
	}
};

#endif