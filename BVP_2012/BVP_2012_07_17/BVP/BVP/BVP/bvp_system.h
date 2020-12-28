#pragma once
#include <array>
#include <functional>
#include "../FunctionApproximation/DerivativeEvaluator/Dual.h"

/// <summary>
/// A structure to represet an argument-value pair of a scalar function of scalar argument
/// </summary>
template <class R>
struct func_value
{
	/// <summary>
	/// Function argument
	/// </summary>
	R x;

	/// <summary>
	/// Function value
	/// </summary>
	R v;
};

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
	std::array<R, varCnt> grad;
};

/// <summary>
/// A structure to holde the "simplest" boundary conditions
/// </summary>
template <class R>
struct bnd_cond_simple
{
	/// <summary>
	/// Boundary condition at the "left" argument point
	/// </summary>
	func_value<R> left;

	/// <summary>
	/// Boundary condition ar the "right" argument point
	/// </summary>
	func_value<R> right;
};

/// <summary>
/// Data structure representing a system of ordinary differential equations of the first order
/// </summary>
template <class R, int eqCnt>
class bvp_system
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

	/// <summary>
	/// Boundary conditions
	/// </summary>
	bnd_cond_simple<R> boundary_conditions;

public:

	/// <summary>
	/// Constructor
	/// </summary>
	bvp_system(const std::array<TFunc, eqCnt>& functions, const bnd_cond_simple<R> bc)
	{
		rhs_functions = functions;
		boundary_conditions = bc;
	}

	/// <summary>
	/// Evaluates right hand side function represented with its id together with the gradient at the given set of arguments
	/// </summary>
	func_value_with_gradient<R, varCnt> Evaluate(const int funcId, const std::array<R, varCnt>& args) const
	{
		if (funcId < 0 || funcId >= eqCnt)
			throw std::exception("Invalid function id");

		std::array<dual<R, varCnt>, varCnt> arguments_dual;
		std::array<R, varCnt> temp{};

		for (int i = 0; i < varCnt; i++)
		{
			temp[i] = 1;
			arguments_dual[i] = dual<R, varCnt>(args[i], temp);
			temp[i] = 0;
		}

		const auto dual_result = rhs_functions[funcId](arguments_dual);

		return { dual_result.Real(), dual_result.Dual() };
	}

	/// <summary>
	/// Getter for the boundary conditions field
	/// </summary>
	bnd_cond_simple<R> get_bc() const
	{
		return boundary_conditions;
	}
};
