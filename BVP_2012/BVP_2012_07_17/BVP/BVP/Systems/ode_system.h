#pragma once
#include <array>
#include <functional>
#include <iostream>
#include <fstream>
#include "../FunctionApproximation/DerivativeEvaluator/Dual.h"


/// <summary>
/// Representation of a mesh point, which is a point on the 
/// solution surface in the (number of unknown functions + 1)-dimensional space
/// </summary>
template <class R, int varCnt>
struct mesh_point
{
	std::array<R, varCnt> pt{};

	/// <summary>
	/// Sub-script operator
	/// </summary>
	R& operator [](const int i)
	{
		return pt[i];
	}

	/// <summary>
	/// sub-script operator (const)
	/// </summary>
	const R& operator [](const int i) const
	{
		return pt[i];
	}

	/// <summary>
	/// Returns maximal absolute value of the mesh point coordinates
	/// </summary>
	R max_abs() const
	{
		return std::abs<R>(*std::max_element(pt.begin(), pt.end(), [](const auto& a, const auto& b) { return std::abs<R>(a) < std::abs<R>(b); }));
	}

	/// <summary>
	/// += operator
	/// </summary>
	mesh_point<R, varCnt>& operator += (const mesh_point<R, varCnt>& arg)
	{
		for (int i = 0; i < varCnt; i++)
			pt[i] += arg.pt[i];

		return *this;
	}

	/// <summary>
	/// -= operator
	/// </summary>
	mesh_point<R, varCnt>& operator -= (const mesh_point<R, varCnt>& arg)
	{
		for (int i = 0; i < varCnt; i++)
			pt[i] -= arg.pt[i];

		return *this;
	}

	/// <summary>
	/// *= operator
	/// </summary>
	mesh_point<R, varCnt>& operator *= (const R& scalar)
	{
		for (int i = 0; i < varCnt; i++)
			pt[i] *= scalar;

		return *this;
	}

	friend std::ostream& operator<<(std::ostream& os, const mesh_point<R, varCnt>& pt)
	{
		for (int varId = 0; varId < varCnt; varId++)
			os << pt[varId] << " ";

		return os;
	}
};

template <class R, int Dim >
void SaveMeshPoints(const char* filename, const std::vector<mesh_point<R, Dim>>& pts)
{
	std::ofstream file;
	file.precision(std::numeric_limits<R>::digits10);
	file.open(filename);
	for (int ptId = 0; ptId < pts.size(); ptId++)
		file << ptId << " " << pts[ptId] << std::endl;

	file.close();
}

template <class R, int varCnt>
mesh_point<R, varCnt> operator + (mesh_point<R, varCnt> lhs, const mesh_point<R, varCnt>& rhs)
{
	return lhs += rhs;
}

template <class R, int varCnt>
mesh_point<R, varCnt> operator - (mesh_point<R, varCnt> lhs, const mesh_point<R, varCnt>& rhs)
{
	return lhs -= rhs;
}

template <class R, int varCnt>
mesh_point<R, varCnt> operator * (mesh_point<R, varCnt> lhs, const R& rhs)
{
	return lhs *= rhs;
}

template <class R, int varCnt>
mesh_point<R, varCnt> operator * (const R& lhs, mesh_point<R, varCnt> rhs)
{
	return rhs *= lhs;
}


/// <summary>
/// Dara structure to contain value of a scalar function together with its gradient at some point
/// </summary>
template <class R, int varCnt>
struct func_value_with_gradient
{
	/// <summary>
	/// Value of the function at the current point
	/// </summary>
	R v;

	/// <summary>
	/// Gradient of the function at the current point
	/// </summary>
	mesh_point<R, varCnt> grad;

	/// <summary>
	/// += operator
	/// </summary>
	func_value_with_gradient<R, varCnt>& operator += (const func_value_with_gradient<R, varCnt>& arg)
	{
		v += arg.v;
		grad += arg.grad;

		return this*;
	}

	/// <summary>
	/// -= operator
	/// </summary>
	func_value_with_gradient<R, varCnt>& operator -= (const func_value_with_gradient<R, varCnt>& arg)
	{
		v -= arg.v;
		grad -= arg.grad;

		return this*;
	}

	/// <summary>
	/// *= operator
	/// </summary>
	func_value_with_gradient<R, varCnt>& operator *= (const R& scalar)
	{
		v *= scalar;
		grad *= scalar;

		return this*;
	}
};


template <class R, int varCnt>
func_value_with_gradient<R, varCnt> operator +(func_value_with_gradient<R, varCnt> lhs, const func_value_with_gradient<R, varCnt>& rhs)
{
	return lhs += rhs;
}

template <class R, int varCnt>
func_value_with_gradient<R, varCnt> operator -(func_value_with_gradient<R, varCnt> lhs, const func_value_with_gradient<R, varCnt>& rhs)
{
	return lhs -= rhs;
}

template <class R, int varCnt>
func_value_with_gradient<R, varCnt> operator *(func_value_with_gradient<R, varCnt> lhs, const R& rhs)
{
	return lhs *= rhs;
}

template <class R, int varCnt>
func_value_with_gradient<R, varCnt> operator *(const R& lhs, func_value_with_gradient<R, varCnt> rhs)
{
	return rhs *= lhs;
}

/// <summary>
/// Data structure representing a system of ordinary differential equations of the first order
/// </summary>
template <class R, int eqCnt>
class ode_system
{
private:
	/// <summary>
	/// Number of variables is equal to the number of equations plus one (for the independent variable)
	/// </summary>
	const static int varCnt = eqCnt + 1;

	/// <summary>
	/// The type for right hand side functions
	/// </summary>
	typedef std::function < dual<R, varCnt>(std::array<dual<R, varCnt>, varCnt >)> TFunc;

	/// <summary>
	/// Right hand side functions
	/// </summary>
	std::array<TFunc, eqCnt> rhs_functions;

public:

	/// <summary>
	/// Right hand side of the system ant its gradients evaluated at some point 
	/// </summary>
	struct eval_result
	{
	private:
		std::array<func_value_with_gradient<R, varCnt>, eqCnt> _values_and_gradients;

	public:

		/// <summary>
		/// The point (set of arguments) at which the right hand side of the system was evaluated
		/// </summary>
		mesh_point<R, varCnt> pt;

		/// <summary>
		/// Subscript operator
		/// </summary>
		func_value_with_gradient<R, varCnt>& operator[](const int i)
		{
			return _values_and_gradients[i];
		}

		/// <summary>
		/// Subscript operator (const version)
		/// </summary>
		const func_value_with_gradient<R, varCnt>& operator[](const int i) const
		{
			return _values_and_gradients[i];
		}
	};

	/// <summary>
	/// Returns number of variables involved (number of equations + 1)
	/// </summary>
	/// <returns></returns>
	int get_var_count() const
	{
		return varCnt;
	}

	/// <summary>
	/// Constructor
	/// </summary>
	ode_system(const std::array<TFunc, eqCnt>& functions) : rhs_functions{ functions }
	{}

	/// <summary>
	/// Evaluates right hand side functions together with their gradient at the given set of arguments
	/// </summary>
	eval_result Evaluate(const mesh_point<R, varCnt>& pt) const
	{
		std::array<dual<R, varCnt>, varCnt> arguments_dual;
		std::array<R, varCnt> temp{};

		for (int i = 0; i < varCnt; i++)
		{
			temp[i] = 1;
			arguments_dual[i] = dual<R, varCnt>(pt.pt[i], temp);
			temp[i] = 0;
		}

		eval_result result;
		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
		{
			const auto result_dual = rhs_functions[eq_id](arguments_dual);
			result[eq_id] = { result_dual.Real(), result_dual.Dual() };
		}

		result.pt = pt;

		return result;
	}
};
