#pragma once
#include <array>
#include <functional>
#include <iostream>
#include <fstream>
#include "poly_function.h"


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
		return auxutils::Abs(*std::max_element(pt.begin(), pt.end(), [](const auto& a, const auto& b) { return auxutils::Abs(a) < auxutils::Abs(b); }));
	}

	/// <summary>
	/// Conditional version of the max abs function, where maximum is taken between the variables that correspond to "true" items in the given map
	/// </summary>
	R max_abs(const std::array<bool, varCnt>& map)
	{
		R result = -std::numeric_limits<R>::max();

		for (int var_id = 0; var_id < varCnt; var_id++)
		{
			if (map[var_id])
				result = std::max<R>(result, auxutils::Abs(pt[var_id]));
		}

		return result;
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

	/// <summary>
	/// Write to stream (in binary format)
	/// </summary>
	friend std::ostream& operator<<(std::ostream& os, const mesh_point<R, varCnt>& pt)
	{
		for (int varId = 0; varId < varCnt; varId++)
			os.write((char*)&pt[varId], sizeof(R));

		return os;
	}

	//Read from stream (in binary format)
	friend std::istream& operator>>(std::istream& is, mesh_point<R, varCnt>& pt)
	{
		for (int varId = 0; varId < varCnt; varId++)
			is.read((char*)&pt[varId], sizeof(R));

		return is;
	}

	/// <summary>
	/// Returns string representation of the current instance
	/// </summary>
	std::string to_string() const
	{
		std::string result;
		for (int varId = 0; varId < varCnt; varId++)
			result += auxutils::to_string_with_precision(pt[varId]) + " ";

		return result;
	}
};

/// <summary>
/// Saves mesh points in text format
/// </summary>
template <class R, int Dim >
void SaveMeshPoints(const std::vector<mesh_point<R, Dim>>& pts, const char* filename)
{
	std::ofstream file;
	file.precision(std::numeric_limits<R>::digits10);
	file.open(filename);
	for (int ptId = 0; ptId < pts.size(); ptId++)
		file << ptId << " " << pts[ptId].to_string() << std::endl;

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
/// Value of a scalar function of multiple scalarar arguments
/// </summary>
template <class R>
struct func_value
{
	/// <summary>
	/// Value of the function at the current point
	/// </summary>
	R v;
};


/// <summary>
/// Dara structure to contain value of a scalar function together with its gradient at some point
/// </summary>
template <class R, int varCnt>
struct func_value_with_gradient : public func_value<R>
{
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
	/// Right hand side functions
	/// </summary>
	std::array<poly_func<R, varCnt>, eqCnt> rhs_functions;

public:

	/// <summary>
	/// A factory to produce the "poly-functions" 
	/// </summary>
	template <class F>
	poly_func<R, varCnt> static create_func(const F& func)
	{
		return poly_func<R, varCnt>::create(func);
	}

	/// <summary>
	/// Right hand side of the system ant its gradients evaluated at some point 
	/// </summary>
	template <class V>
	struct eval_result_base
	{
	private:
		std::array<V, eqCnt> _values;

	public:

		/// <summary>
		/// The point (set of arguments) at which the right hand side of the system was evaluated
		/// </summary>
		mesh_point<R, varCnt> pt;

		/// <summary>
		/// Subscript operator
		/// </summary>
		V& operator[](const int i)
		{
			return _values[i];
		}

		/// <summary>
		/// Subscript operator (const version)
		/// </summary>
		const V& operator[](const int i) const
		{
			return _values[i];
		}
	};

	using eval_result = eval_result_base<func_value_with_gradient<R, varCnt>>;

	using eval_result_minimal = eval_result_base<func_value<R>>;

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
	ode_system(const std::array<poly_func<R, varCnt>, eqCnt>& functions) : rhs_functions{ functions }
	{}

	/// <summary>
	/// Evaluates right hand side functions together with their gradient at the given set of arguments
	/// </summary>
	eval_result evaluate(const mesh_point<R, varCnt>& pt) const
	{
		std::array<dual<R, varCnt>, varCnt> arguments_dual;
		std::array<R, varCnt> temp{};

		for (int i = 0; i < varCnt; i++)
		{
			temp[i] = R(1);
			arguments_dual[i] = dual<R, varCnt>(pt.pt[i], temp);
			temp[i] = R(0);
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

	/// <summary>
	/// Evaluates right hand side functions together with their gradient at the given set of arguments
	/// </summary>
	eval_result_minimal evaluate_minimal(const mesh_point<R, varCnt>& pt) const
	{
		eval_result_minimal result;
		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
			result[eq_id] = { rhs_functions[eq_id].get_func()(pt.pt) };

		result.pt = pt;

		return result;
	}
};
