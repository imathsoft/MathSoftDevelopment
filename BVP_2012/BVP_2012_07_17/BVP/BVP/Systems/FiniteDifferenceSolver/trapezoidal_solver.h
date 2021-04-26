#pragma once
#include <vector>
#include <array>
#include <exception>
#include <algorithm>
#include "../ode_system.h"
#include "../../LinearAlgebra/Matrix.h"
#include "../../Utils/AuxUtils.h"

/// <summary>
/// A "boundary conditions marker" is a simple way to define some sub-set of two-point boundary conditions
/// when values of (some) unknown functions are explicitely given at the endpoint(-s)
/// The marker consists of two sub-parts: one for the "left" endpoint and another one for the "right" endpoint
/// each containing bool values defining which function value is given (as a boundary condition) and which is not
/// </summary>
template <int eqCnt>
struct bc_marker
{
	std::array<bool, eqCnt> left_marker{};
	std::array<bool, eqCnt> right_marker{};

	/// <summary>
	/// It is supposed that the number of defined values for both endpoints should be equal to number of equations
	/// </summary>
	bool is_valid() const
	{
		return eqCnt == std::count(left_marker.begin(), left_marker.end(), true) + std::count(right_marker.begin(), right_marker.end(), true);
	}
};

/// <summary>
/// Data structure to contain information needed to diagnose problems in the convergence rate of the Newton's method
/// </summary>
template <class R>
struct ConvergenceInfo
{
	/// <summary>
	/// Magnitude of Newtn's correction on the given iteration
	/// </summary>
	R CorrectionMagnitude;

	/// <summary>
	/// Indicates whether mesh refinement was applied on the given iteration
	/// </summary>
	bool RefinementApplied;

	/// <summary>
	/// Write to stream
	/// </summary>
	friend std::ostream& operator<<(std::ostream& os, const ConvergenceInfo& conv_info)
	{
		os << conv_info.CorrectionMagnitude << " " << conv_info.RefinementApplied;

		return os;
	}
};

/// <summary>
/// Finite difference solver implementing the trapezoidal scheme for systems of nonlinear 
/// first order ordinary differential euations
/// </summary>
template <class R, int eqCnt = 2>
class trapezoidal_solver
{
	static_assert(eqCnt == 2, "Does not support other number of equaions");

	typedef ode_system<R, eqCnt> sys;

	/// <summary>
	/// Data describing a transformation that was applied to the system in question on one of the discretization intervals
	/// </summary>
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
	};

	/// <summary>
	/// A struct representing extended atrix of the form [m|b], where "m"
	/// is a square matrix and "b" is a column-vector of the corresponding size
	/// </summary>
	struct stripe
	{
		LinAlg::Matrix<R, eqCnt, eqCnt> m;
		LinAlg::Matrix<R, eqCnt, 1> b;
		transformation_maker trans_marker{};
	};

	/// <summary>
	/// Type to hold extended block-diagonal matrix
	/// </summary>
	typedef std::vector<stripe> block_matrix;

	inline static const R OneHalf = R(1) / R(2);

	/// <summary>
	/// A collectiion that gets filled during the Newton's iteration process 
	///  and contains magnitudes of each Newton's correction
	/// Setves diagnostic purposes (if the iteration process is properly implemented then the corrctions should exhibit quadratic decay)
	/// </summary>
	std::vector<ConvergenceInfo<R>> correction_magnitudes{};

	int max_iterations_count{ 15 };

	/// <summary>
	/// Field to store "status" of the previous call of the "solve" method
	/// </summary>
	bool succeeded{true};

	/// <summary>
	/// A data struct to hold evaluation result transformed via a "pivot" transformation
	/// </summary>
	struct eval_result_transformed
	{
		/// <summary>
		/// Transformed result
		/// </summary>
		typename ode_system<R, eqCnt>::eval_result res{};

		/// <summary>
		/// Marker defining a transformation
		/// </summary>
		transformation_maker trans_marker{};
	};

	/// <summary>
	/// Performs transformation by proper change of independent variable
	/// </summary>
	static typename ode_system<R, eqCnt>::eval_result transform_independent_var(const typename ode_system<R, eqCnt>::eval_result& res_to_transform, const int independent_var_id)
	{
		if (independent_var_id < 0 || independent_var_id >= eqCnt)
			return res_to_transform; // identity transformation

		const auto& pivot_eq_data = res_to_transform[independent_var_id];

		const auto one_over_denominator = R(1) / (pivot_eq_data.v * pivot_eq_data.v);

		auto res_transformed = res_to_transform;

		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
		{
			auto& data = res_transformed[eq_id];

			if (eq_id == independent_var_id)
			{
				data.v = R(1) / pivot_eq_data.v;

				for (int var_id = 0; var_id <= eqCnt; var_id++)
					data.grad[var_id] = -pivot_eq_data.grad[var_id] * one_over_denominator;

				continue;
			}

			const auto& curr_eq_data = res_to_transform[eq_id];

			data.v = curr_eq_data.v / pivot_eq_data.v;

			for (int var_id = 0; var_id <= eqCnt; var_id++)
				data.grad[var_id] = (curr_eq_data.grad[var_id] * pivot_eq_data.v - curr_eq_data.v * pivot_eq_data.grad[var_id]) * one_over_denominator;
		}

		return res_transformed;
	}

	/// <summary>
	/// Performs transformation by inversing unknown variables according to the given inversion map
	/// </summary>
	static typename ode_system<R, eqCnt>::eval_result invert_unknowns(const typename ode_system<R, eqCnt>::eval_result& res_to_transform, std::array<bool, eqCnt> inversion_map)
	{
		auto res_transformed = res_to_transform;
		const auto& pt = res_to_transform.pt;

		//Iterate through the derivatives of unknown functions (i.e., equations) 
		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
		{
			if (!inversion_map[eq_id])
				continue; //skip those equations that do not correspond to inverted unknown functions

			auto& eq_data = res_transformed[eq_id];
			const auto factor = R(1) / pt[eq_id];
			const auto minus_factor_squared = -factor * factor;

			eq_data.v *= minus_factor_squared;
			for (int var_id = 0; var_id <= eqCnt; var_id++)
			{
				if (eq_id != var_id)
					eq_data.grad[var_id] *= minus_factor_squared;
				else
					eq_data.grad[var_id] -= 2 * res_to_transform[eq_id].v * factor;
			}
		}

		//Iterate though the unknown functions
		for (int var_id = 0; var_id < eqCnt; var_id++)
		{
			if (!inversion_map[var_id])
				continue; //skip those that are not inverted

			const auto factor = - pt[var_id] * pt[var_id];

			for (int eq_id = 0; eq_id < eqCnt; eq_id++)
			{
				if (eq_id == var_id)
					continue; //skip the equation that corresponds to the current inverted unknown function, since all of them have been already processed above

				res_transformed[eq_id].grad[var_id] *= factor;
			}
		}

		return res_transformed;
	}

	/// <summary>
	/// Performs transformation of the system evalueation result according to the given "pivot" variable
	/// </summary>
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res_to_transform, const transformation_maker& trans_marker)
	{
		const auto res_independent_var_transformed = transform_independent_var(res_to_transform, trans_marker.pivot_id);

		const auto res_unknowns_inverted = invert_unknowns(res_independent_var_transformed, trans_marker.inversion_map);

		return { res_unknowns_inverted, trans_marker };
	}

	/// <summary>
	/// Performs transformation of the system evaluation result based on the transformation map and the value of the corresponding derivatives
	/// Transformation map defines what unknownc can be chosen as "pivot" ones (to be swaped with the independent variable via inverting)
	/// Returns intex of the chosen pivot variable
	/// </summary>
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res,
		const std::array<bool, eqCnt>& transform_map, const bool use_inversion, const R& derivative_threshold)
	{
		int independent_var_id = eqCnt;
		std::array<bool, eqCnt> inversion_map{};
		const auto& pt = res.pt;

		R max_abs_val = derivative_threshold;

		for (int var_id = 0; var_id < eqCnt; var_id++)
		{
			const auto trial_abs_val = auxutils::Abs(res[var_id].v);

			if (trial_abs_val > max_abs_val && transform_map[var_id])//only variable that is "marked" in the transformation map can serve as "independent"
			{
				independent_var_id = var_id;
				max_abs_val = trial_abs_val;
			} else if (use_inversion && !transform_map[var_id] &&
				trial_abs_val > derivative_threshold && (auxutils::Abs(pt[var_id]) > R(1)))//there is no point in applying inversion if the absolute
																						 //value of the corresponding variavle is less or equal to 1
			{
				inversion_map[var_id] = true;
			}
		}

		return transform(res, {independent_var_id, inversion_map});
	}

	/// <summary>
	/// Calculates difference of the two mesh points according to the given inversion map
	/// </summary>
	static mesh_point<R, eqCnt + 1> get_mesh_point_diff(const mesh_point<R, eqCnt + 1>& pt_prev,
		const mesh_point<R, eqCnt + 1>& pt_next, const std::array<bool, eqCnt>& inversion_map)
	{
		mesh_point<R, eqCnt + 1> result{};

		for (int unknown_id = 0; unknown_id < eqCnt; unknown_id++)
		{
			result[unknown_id] = inversion_map[unknown_id] ?
				(R(1) / pt_next[unknown_id] - R(1) / pt_prev[unknown_id]) :
				(pt_next[unknown_id] - pt_prev[unknown_id]);
		}

		result[eqCnt] = pt_next[eqCnt] - pt_prev[eqCnt];

		return result;
	}

	/// <summary>
	/// Calculates "agreement factor" that is aimed to handle the "transition" situatuion, when inversion transformation is applied on the current interval
	/// but was not applied on the previous interval (or vise versa). This transitions state requires us to still work with not transformed variable 
	/// despite of the fact that equations that we built correspond to the transformed one.
	/// In the other words, we have calculated partial derivative that corresponds to the transformed variable `w = 1/v`, however we need to have 
	/// the partial derivative with respect to `v`. This can be handled by multiplying the partial derivative with respect to `w` by `-1/v^2`.
	/// The latter is actually the "agreement" factor that the function calculates with respect to each unknown
	/// </summary>
	static std::array<R, eqCnt> calc_inversion_agreement_factor(const std::array<bool, eqCnt>& inversion_map_prev,
		const std::array<bool, eqCnt>& inversion_map_next,
		const mesh_point<R, eqCnt + 1>& current_pt)
	{
		std::array<R, eqCnt> result{};

		for (int unknown_id = 0; unknown_id < eqCnt; unknown_id++)
		{
			if (inversion_map_prev[unknown_id] == inversion_map_next[unknown_id])
			{
				result[unknown_id] = R(1);
				continue;
			}

			if (inversion_map_prev[unknown_id])
				result[unknown_id] = -current_pt[unknown_id] * current_pt[unknown_id];
			else {
				const auto temp = R(1) / current_pt[unknown_id];
				result[unknown_id] = -temp * temp;
			}
		}

		return result;
	}

	/// <summary>
	/// Constructs the extended gradient matrix, which is essentially
	/// [F'(s_{i}) | F(s_{i})],
	/// where F(s) = 0 is the system of nonliner equations that we get applying the trapezoidal scheme to the given system of ODEs 
	/// s_{i} is the "initial guess"
	/// </summary>
	static block_matrix construct_extended_gradient_matrix(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess,
		const std::array<bool, eqCnt>& transform_map, const bool use_inversion, const R& derivative_threshold)
	{
		if (init_guess.size() <= 1)
			throw std::exception("Invalid input");

		block_matrix result(init_guess.size() - 1);

		auto res_prev = system.evaluate(init_guess[0]);
		transformation_maker trans_marker_prev{ -1 };
		for (auto pt_id = 0; pt_id < init_guess.size() - 1; pt_id++)
		{
			const auto res_prev_transformed = transform(res_prev, transform_map, use_inversion, derivative_threshold);
			const auto trans_marker_next = res_prev_transformed.trans_marker;
			if (trans_marker_prev.pivot_id < 0)
				trans_marker_prev = trans_marker_next;

			const auto res_next = system.evaluate(init_guess[pt_id + 1]);
			const auto res_next_transformed = transform(res_next, trans_marker_next);

			const auto independent_var_id_next = trans_marker_next.pivot_id;
			const auto& inversion_map_next = trans_marker_next.inversion_map;

			const auto independent_var_id_prev = trans_marker_prev.pivot_id;
			const auto& inversion_map_prev = trans_marker_prev.inversion_map;

			const auto pt_diff = get_mesh_point_diff(init_guess[pt_id], init_guess[pt_id + 1], inversion_map_next);
			const auto one_over_h = R(1) / pt_diff.pt[independent_var_id_next];
			const auto divided_diff = one_over_h * pt_diff;

			const auto agreement_factor = calc_inversion_agreement_factor(inversion_map_prev, inversion_map_next, init_guess[pt_id]);

			auto& current_stripe = result[pt_id];
			current_stripe.trans_marker = trans_marker_next;
			LinAlg::Matrix<R, eqCnt, eqCnt> temp;

			for (int row_id = 0; row_id < eqCnt; row_id++)
			{
				for (int col_id = 0; col_id < eqCnt; col_id++)
				{
					current_stripe.m[row_id][col_id] = -OneHalf * agreement_factor[col_id] * res_prev_transformed.res[row_id].grad[col_id != independent_var_id_prev ? col_id : eqCnt];
					temp[row_id][col_id] = -OneHalf * res_next_transformed.res[row_id].grad[col_id != independent_var_id_next ? col_id : eqCnt];
				}

				const int actual_var_id = row_id != independent_var_id_next ? row_id : eqCnt;

				current_stripe.m[row_id][row_id] -= actual_var_id != independent_var_id_prev ? agreement_factor[row_id] * one_over_h : R(0);
				temp[row_id][row_id] += one_over_h;

				current_stripe.b[row_id][0] = divided_diff[actual_var_id] - OneHalf * (res_prev_transformed.res[row_id].v + res_next_transformed.res[row_id].v);
			}

			if (independent_var_id_prev != independent_var_id_next)//"transition" step
			{
				const int swapped_var_col_id = independent_var_id_next != eqCnt ? independent_var_id_next : independent_var_id_prev;
				for (int row_id = 0; row_id < eqCnt; row_id++)
				{
					const int actual_var_id = row_id != independent_var_id_next ? row_id : eqCnt;
					current_stripe.m[row_id][swapped_var_col_id] += divided_diff[actual_var_id] * one_over_h;
				}
			}

			const auto det = temp.Determinant();

			if (auxutils::Abs(det) < 1e3 * std::numeric_limits<R>::epsilon())
				throw std::exception("Singular block");

			const auto temp_inverted = temp.Inverse();

			current_stripe.m = -temp_inverted * current_stripe.m;//inverse sign so that now matrix "m" and vector "b" are on the "same side"
			current_stripe.b = temp_inverted * current_stripe.b;

			res_prev = res_next;
			trans_marker_prev = trans_marker_next;
		}

		return result;
	}

	/// <summary>
	/// Representation of a mesh point correction
	/// </summary>
	struct correction
	{
		/// <summary>
		/// The correction itself
		/// </summary>
		mesh_point<R, eqCnt + 1> point;

		/// <summary>
		/// Transformation marker
		/// </summary>
		transformation_maker trans_marker{};

		/// <summary>
		/// Correction magnitude
		/// </summary>
		R magnitude() const
		{
			return point.max_abs();
		}
	};

	/// <summary>
	/// A helper method to copy result from vector-column to the mesh point
	/// </summary>
	static void copy(correction& dest, const LinAlg::Matrix<R, eqCnt, 1>& source, const transformation_maker& trans_marker)
	{
		dest.trans_marker = trans_marker;
		for (int i = 0; i < eqCnt; i++)
			dest.point[trans_marker.pivot_id != i ? i : eqCnt] = source[i][0];
	}

	/// <summary>
	/// Returns correction for the left endpoint.
	/// The calculations are done based on the input matrix "m" and column-vector "b" as well
	/// as on the input "boundary conditions marker" (see summary of the corresponding data struct).
	/// Without boundary conditions, the system that we are working with is underdetermined 
	/// (eqCnt equations with 2*eqCnt unknowns) and its extended matrix looks like this:
	/// [-m,-I|b], where "I" is the identity matrix of the corresponding dimension.
	/// Boundary marker allows us to "discard" those columns of -[m.I] that correspond to unknowns 
	/// specified via the boundary conditions, so that eventually we end up with a
	/// system of eqCnt equations with eqCnt unknowns
	/// </summary>
	static LinAlg::Matrix<R, eqCnt, 1> resolve_endpoint_corrections(
		const stripe& final_block, const bc_marker<eqCnt>& bcm)
	{
		LinAlg::Matrix<R, eqCnt, eqCnt> temp{};
		LinAlg::Matrix<R, eqCnt, 1> result{};
		std::vector<R*> map{};

		int col_id = 0;
		const auto& m = final_block.m;
		for (int m_id = 0; m_id < bcm.left_marker.size(); m_id++)
		{
			if (bcm.left_marker[m_id])
				continue;

			map.push_back(&result[m_id][0]);

			for (int row_id = 0; row_id < eqCnt; row_id++)
				temp[row_id][col_id] = m[row_id][m_id];

			col_id++;
		}

		for (int m_id = 0; m_id < bcm.right_marker.size(); m_id++)
		{
			if (bcm.right_marker[m_id])
				continue;

			temp[m_id][col_id] = R(1);

			col_id++;
		}

		if (col_id != eqCnt)//sanity check
			throw std::exception("Something went wrong");

		const auto det = temp.Determinant();

		if (auxutils::Abs(det) < 100 * std::numeric_limits<R>::epsilon())
			throw std::exception("Singular system");

		const auto temp_inverted = temp.Inverse();

		const auto solution = - temp_inverted * final_block.b;

		for (int i = 0; i < map.size(); i++)
			*map[i] = solution[i][0];

		return result;
	}

	/// <summary>
	/// Solves system of linear equation wiith respect to the correction of the Newton's method
	/// The correction is returned
	/// Modifies the input gradient matrix 
	/// </summary>
	static std::vector<correction> get_newton_correction(block_matrix& gradient_matrix, const bc_marker<eqCnt>& bcm)
	{
		for (auto block_id = 1; block_id < gradient_matrix.size(); block_id++)
		{
			gradient_matrix[block_id].b += gradient_matrix[block_id].m * gradient_matrix[block_id - 1].b;
			gradient_matrix[block_id].m = gradient_matrix[block_id].m * gradient_matrix[block_id - 1].m;
		}

		static_assert(eqCnt == 2, "Current implementation supports only systems of ODEs with 2 equations and 2 unknown functions.");

		const auto& final_block = *gradient_matrix.rbegin();
		const auto u_0 = resolve_endpoint_corrections(final_block, bcm);

		std::vector<correction> result(gradient_matrix.size() + 1);

		copy(result[0], u_0, gradient_matrix[0].trans_marker);

		for (auto block_id = 0; block_id < gradient_matrix.size(); block_id++)
			copy(result[block_id + 1], gradient_matrix[block_id].m * u_0 + gradient_matrix[block_id].b, gradient_matrix[block_id].trans_marker);

		return result;
	}

	/// <summary>
	/// Applies correction to the given solution
	/// Returns true if actual refinement was applied to the corrected solution
	/// </summary>
	static bool apply_correction_and_refine_solution(std::vector<mesh_point<R, eqCnt + 1>>& solution, const std::vector<correction>& correction, const R desired_step_size)
	{
		if (solution.size() != correction.size())
			throw std::exception("Invalid input");

		for (auto pt_id = 0; pt_id < correction.size(); pt_id++)
		{
			const auto& marker = correction[pt_id].trans_marker;
			for (int var_id = 0; var_id < eqCnt; var_id++)
			{
				if (marker.inversion_map[var_id])
					solution[pt_id][var_id] = solution[pt_id][var_id] /(R(1) - solution[pt_id][var_id] * correction[pt_id].point[var_id]);
				else
					solution[pt_id][var_id] -= correction[pt_id].point[var_id];
			}
			solution[pt_id][eqCnt] -= correction[pt_id].point[eqCnt];
		}

		std::vector<mesh_point<R, eqCnt + 1>> solution_refined;
		solution_refined.reserve(solution.size());

		const auto& argument_min = solution[0][2];
		const auto& argument_max = solution[solution.size() - 1][2];

		bool refinement_applied = false;

		solution_refined.push_back(solution[0]);

		for (auto pt_id = 1; pt_id < solution.size(); pt_id++)
		{
			const auto prev_pt = *solution_refined.rbegin();

			if (solution[pt_id][2] < argument_min || solution[pt_id][2] > argument_max || prev_pt[2] >= solution[pt_id][2])
			{
				refinement_applied = true;
				continue;
			}

			const auto ind_var_id = correction[pt_id].trans_marker.pivot_id;
			const auto actual_step = auxutils::Abs(prev_pt[ind_var_id] - solution[pt_id][ind_var_id]);
			if (actual_step < R(2) * desired_step_size)
			{
				if (actual_step > R(0.1) * desired_step_size)
					solution_refined.push_back(solution[pt_id]);
				else
					refinement_applied = true;

				continue;
			}

			refinement_applied = true;

			const int points_to_add = static_cast<int>(actual_step / desired_step_size);
			const auto h = R(1) / (points_to_add + 1);
			const auto position_increment = (solution[pt_id] - prev_pt)*h;

			auto temp = prev_pt;
			for (int extra_pt_id = 0; extra_pt_id < points_to_add; extra_pt_id++)
			{
				temp += position_increment;
				solution_refined.push_back(temp);
			}

			solution_refined.push_back(solution[pt_id]);
		}

		solution = solution_refined;

		return refinement_applied;
	}

public:

	/// <summary>
	/// Getter for the "status" field
	/// </summary>
	bool success() const
	{
		return succeeded;
	}

	/// <summary>
	/// Getter for the maximal iterations count field
	/// </summary>
	R& max_iter_count()
	{
		return max_iterations_count;
	}

	/// <summary>
	/// Getter for the collection of correction magnitudes
	/// </summary>
	std::vector<ConvergenceInfo<R>> get_correcion_magnitudes() const
	{
		return correction_magnitudes;
	}

	/// <summary>
	/// The main solving method
	/// </summary>
	std::vector<mesh_point<R, eqCnt + 1>> solve(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const bc_marker<eqCnt>& bcm,
		const R& precision, const R& desired_step, const bool use_reparametrization, const bool use_inversion, const R& derivative_threshold = R(1))
	{
		if (!bcm.is_valid())
			throw std::exception("Invalid boundary condition marker");

		std::array<bool, eqCnt> transform_map{};

		if (use_reparametrization)
			for (int var_id = 0; var_id < eqCnt; var_id++)
				transform_map[var_id] = bcm.left_marker[var_id] && bcm.right_marker[var_id];

		auto solution = init_guess;
		correction_magnitudes.clear();

		R correction_magnitude = std::numeric_limits<R>::max();

		int iter_count = 0;
		while (correction_magnitude > precision && iter_count < max_iterations_count)
		{
			auto g_matrix = construct_extended_gradient_matrix(system, solution, transform_map, use_inversion, derivative_threshold);
			const auto correction = get_newton_correction(g_matrix, bcm);
			correction_magnitude = std::max_element(correction.begin(), correction.end(),
				[](const auto& a, const auto& b) { return a.magnitude() < b.magnitude(); })->magnitude();

			const auto refinement_applied = apply_correction_and_refine_solution(solution, correction, desired_step);

			correction_magnitudes.emplace_back( ConvergenceInfo<R>{ correction_magnitude, refinement_applied });

			iter_count++;
		}

		succeeded = (correction_magnitude <= precision);

		return solution;
	}
};