#pragma once
#include <functional>
#include "../FunctionApproximation/DerivativeEvaluator/Dual.h"

/// <summary>
/// A "poly-function" implementation
/// i.e. a function that can take different arguments types and return results of different types respectively 
/// </summary>
template <class R, int varCnt>
struct poly_func
{
private:
	/// <summary>
	/// Private default constructor, so that an instance of this class can be created
	/// exclusively through the factory method below
	/// </summary>
	poly_func() = default;

	using func_dual_t = std::function < dual<R, varCnt>(std::array<dual<R, varCnt>, varCnt >)>;
	using func_t = std::function <R(std::array<R, varCnt >)>;

	func_dual_t func_dual{};
	func_t func{};

	/// <summary>
	/// Two parameter constructor
	/// </summary>
	poly_func(const func_dual_t& f_dual, const func_t& f)
	{
		func_dual = f_dual;
		func = f;
	}

public:

	/// <summary>
	/// A facorty : the only way how an instance of the object can be created
	/// </summary>
	template<class F>
	static poly_func<R, varCnt> create(const F& func)
	{
		return poly_func<R, varCnt>(func, func);
	}

	/// <summary>
	/// Getter for the dual function
	/// </summary>
	const func_dual_t& get_gual() const
	{
		return func_dual;
	}

	/// <summary>
	/// Getter for the regular function
	/// </summary>
	const func_t& get_func() const
	{
		return func;
	}

	/// <summary>
	/// Operator calling the "dual" function
	/// </summary>
	dual<R, varCnt> operator()(const std::array<dual<R, varCnt>, varCnt>& args) const
	{
		return func_dual(args);
	}
};
