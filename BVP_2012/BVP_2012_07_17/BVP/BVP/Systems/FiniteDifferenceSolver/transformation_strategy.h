#pragma once
#include "../ode_system.h"

/// <summary>
/// Data describing a transformation that was applied to the system in question on one of the discretization intervals
/// </summary>
template <int eqCnt>
struct transformation_maker
{
	/// <summary>
	/// Index of the "pivot" unknown
	/// </summary>
	int pivot_id{};

	/// <summary>
	/// The map defining which unknowns have been inverted during the transformation
	/// </summary>
	std::array<bool, eqCnt> inversion_map{};

	/// <summary>
	/// Equality operator
	/// </summary>
	bool operator ==(const transformation_maker& another_marker) const
	{
		return another_marker.pivot_id == pivot_id && another_marker.inversion_map == inversion_map;
	}

	/// <summary>
	/// Returns string representation of the current instance
	/// </summary>
	std::string to_string() const
	{
		std::string result;
		result += std::to_string(pivot_id) + " ";

		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
			result += std::string(inversion_map[eq_id] ? "True" : "False") + " ";

		return result;
	}
};

/// <summary>
/// Restrictions which are to be teken into account when deciding about transformation to be applied
/// </summary>
template <class R, int eqCnt>
struct transform_restrictions
{
	/// <summary>
	/// A map defining what variables can be considered as independent when doing a transformation (i.e., "swap" transformations)
	/// </summary>
	std::array<bool, eqCnt> swap_map{};

	/// <summary>
	/// A map defining what variables can be inverted ("flipped")
	/// </summary>
	std::array<bool, eqCnt> flip_map{};

	/// <summary>
	/// Derivative threshold using to make a decision about what variable should be treated as "inrependent"
	/// </summary>
	R derivative_threshold{};
};


/// <summary>
/// A "standard" transformation strategy
/// </summary>
class ts_standard
{
	
public:
	/// <summary>
	/// Returns the transformation marker for the given result of evaluatioin of the right hand side of the given system
	/// </summary>
	template <class R, class V, int eqCnt>
	static transformation_maker<eqCnt> get_transform_marker(const typename ode_system<R, eqCnt>::template eval_result_base<V>& res, const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		int independent_var_id = eqCnt;
		std::array<bool, eqCnt> inversion_map{};
		const auto& pt = res.pt;

		R max_abs_val = trans_restrict.derivative_threshold;

		for (int var_id = 0; var_id < eqCnt; var_id++)
		{
			const auto trial_abs_val = auxutils::Abs(res[var_id].v);

			if (trial_abs_val > max_abs_val && trans_restrict.swap_map[var_id])//only variable that is "marked" in the "flip" map can serve as "independent"
			{
				independent_var_id = var_id;
				max_abs_val = trial_abs_val;
			}
			else if (trans_restrict.flip_map[var_id] &&
				(auxutils::Abs(pt[var_id]) > R(1)))//there is no point in applying inversion ("flip" transformation) if the absolute
												   //value of the corresponding variavle is less or equal to 1
			{
				inversion_map[var_id] = true;
			}
		}

		return { independent_var_id, inversion_map };
	}
};

/// <summary>
/// A transformation strategy designed specifically for the systems of 2 ODEs that represents a second order ordinary differential equaion 
/// (under some additional assumptions)
/// </summary>
class ts_experimental
{
public:
	/// <summary>
	/// Returns the transformation marker for the given result of evaluatioin of the right hand side of the given system
	/// </summary>
	template <class R, class V, int eqCnt>
	static transformation_maker<eqCnt> get_transform_marker(const typename ode_system<R, eqCnt>::template eval_result_base<V>& res, const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		static_assert(eqCnt == 2, "Unexpected number of equations");

		int independent_var_id = eqCnt;
		std::array<bool, eqCnt> inversion_map{};
		const auto& pt = res.pt;

		if (auxutils::Abs(pt[1]) > R(1) && auxutils::Abs(pt[1] * pt[1] * pt[1]) > auxutils::Abs(res[1].v))
			return { 0, {false, true} };

		if (auxutils::Abs(res[1].v) > R(1) && pt[2] > R(0.2) && pt[2] < R(0.7))
			return { 1, {false, false} };

		return { 2, {false, false} };
	}
};

