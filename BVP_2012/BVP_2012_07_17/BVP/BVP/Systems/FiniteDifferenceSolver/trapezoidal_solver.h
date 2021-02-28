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
/// Finite difference solver implementing the trapezoidal scheme for systems of nonlinear 
/// first order ordinary differential euations
/// </summary>
template <class R, int eqCnt = 2>
class trapezoidal_solver
{
	static_assert(eqCnt == 2, "Does not support other number of equaions");

	typedef ode_system<R, eqCnt> sys;

	/// <summary>
	/// A struct representing extended atrix of the form [m|b], where "m"
	/// is a square matrix and "b" is a column-vector of the corresponding size
	/// </summary>
	struct stripe
	{
		LinAlg::Matrix<R, eqCnt, eqCnt> m;
		LinAlg::Matrix<R, eqCnt, 1> b;
		int independent_var_id;
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
	std::vector<R> correction_magnitudes{};

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
		typename ode_system<R, eqCnt>::eval_result res;

		/// <summary>
		/// Index of the "pivot" unknown
		/// </summary>
		int pivot_id;
	};

	/// <summary>
	/// Performs transformation of the system evalueation result according to the given "pivot" variable
	/// </summary>
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res_to_transform, const int independent_var_id)
	{
		if (independent_var_id < 0 || independent_var_id >= eqCnt)
			return { res_to_transform,  eqCnt }; // identity transformation

		const auto& pivot_eq_data = res_to_transform[independent_var_id];

		const auto one_over_denominator = R(1) / (pivot_eq_data.v * pivot_eq_data.v);

		ode_system<R, eqCnt>::eval_result res_transformed{};

		for (int eq_id = 0; eq_id < eqCnt; eq_id++)
		{
			auto& data = res_transformed[eq_id];

			if (eq_id == independent_var_id)
			{
				data.v = R(1) / pivot_eq_data.v;

				for (int var_id = 0; var_id <= eqCnt; var_id++)
					data.grad[var_id] = - pivot_eq_data.grad[var_id] * one_over_denominator;

				continue;
			}

			const auto& curr_eq_data = res_to_transform[eq_id];

			data.v = curr_eq_data.v / pivot_eq_data.v;

			for (int var_id = 0; var_id <= eqCnt; var_id++)
				data.grad[var_id] = (curr_eq_data.grad[var_id] * pivot_eq_data.v - curr_eq_data.v * pivot_eq_data.grad[var_id]) * one_over_denominator;
		}

		return { res_transformed, independent_var_id };
	}

	/// <summary>
	/// Performs transformation of the system evaluation result based on the transformation map and the value of the corresponding derivatives
	/// Transformation map defines what unknownc can be chosen as "pivot" ones (to be swaped with the independent variable via inverting)
	/// Returns intex of the chosen pivot variable
	/// </summary>
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res, const std::array<bool, eqCnt>& transform_map, const R& derivative_threshold = R(1))
	{
		int independent_var_id = eqCnt;

		R max_abs_val = derivative_threshold;

		for (int var_id = 0; var_id < eqCnt; var_id++)
		{
			if (!transform_map[var_id])
				continue;

			const auto trial_abs_val = std::abs<R>(res[var_id].v);

			if (trial_abs_val > max_abs_val)
			{
				independent_var_id = var_id;
				max_abs_val = trial_abs_val;
			}
		}

		return transform(res, independent_var_id);
	}

	/// <summary>
	/// Constructs the extended gradient matrix, which is essentially
	/// [F'(s_{i}) | F(s_{i})],
	/// where F(s) = 0 is the system of nonliner equations that we get applying the trapezoidal scheme to the given system of ODEs 
	/// s_{i} is the "initial guess"
	/// </summary>
	static block_matrix construct_extended_gradient_matrix(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const std::array<bool, eqCnt>& transform_map)
	{
		if (init_guess.size() <= 1)
			throw std::exception("Invalid input");

		block_matrix result(init_guess.size() - 1);

		auto res_prev = system.Evaluate(init_guess[0]);
		int independent_var_id_prev = -1;
		for (int pt_id = 0; pt_id < init_guess.size() - 1; pt_id++)
		{
			const auto res_prev_transformed = transform(res_prev, transform_map);
			const int independent_var_id_next = res_prev_transformed.pivot_id;
			if (independent_var_id_prev < 0)
				independent_var_id_prev = independent_var_id_next;

			const auto res_next = system.Evaluate(init_guess[pt_id + 1]);
			const auto res_next_transformed = transform(res_next, independent_var_id_next);

			const auto pt_diff = init_guess[pt_id + 1] - init_guess[pt_id];
			const auto one_over_h = R(1) / pt_diff.pt[independent_var_id_next];
			const auto divided_diff = one_over_h * pt_diff;

			auto& current_stripe = result[pt_id];
			current_stripe.independent_var_id = independent_var_id_next;
			LinAlg::Matrix<R, eqCnt, eqCnt> temp;

			for (int row_id = 0; row_id < eqCnt; row_id++)
			{
				for (int col_id = 0; col_id < eqCnt; col_id++)
				{
					current_stripe.m[row_id][col_id] = -OneHalf * res_prev_transformed.res[row_id].grad[col_id != independent_var_id_prev ? col_id : eqCnt];
					temp[row_id][col_id] = -OneHalf * res_next_transformed.res[row_id].grad[col_id != independent_var_id_next ? col_id : eqCnt];
				}

				const int actual_var_id = row_id != independent_var_id_next ? row_id : eqCnt;

				current_stripe.m[row_id][row_id] -= actual_var_id != independent_var_id_prev ? one_over_h : R(0);
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

			if (std::abs<R>(det) < 1e3 * std::numeric_limits<R>::epsilon())
				throw std::exception("Singular block");

			const auto temp_inverted = temp.Inverse();

			current_stripe.m = -temp_inverted * current_stripe.m;//inverse sign so that now matrix "m" and vector "b" are on the "same side"
			current_stripe.b = temp_inverted * current_stripe.b;

			res_prev = res_next;
			independent_var_id_prev = independent_var_id_next;
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
		/// Index of independent variable used on the interval to the right
		/// </summary>
		int independent_var_id;

		/// <summary>
		/// Correction magnitude
		/// </summary>
		/// <returns></returns>
		R magnitude() const
		{
			return point.max_abs();
		}
	};

	/// <summary>
	/// A helper method to copy result from vector-column to the mesh point
	/// </summary>
	static void copy(correction& dest, const LinAlg::Matrix<R, eqCnt, 1>& source, const int independent_var_id)
	{
		dest.independent_var_id = independent_var_id;
		for (int i = 0; i < eqCnt; i++)
		{
			dest.point[independent_var_id != i ? i : eqCnt] = source[i][0];
		}
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

		if (std::abs<R>(det) < 100 * std::numeric_limits<R>::epsilon())
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
		for (int block_id = 1; block_id < gradient_matrix.size(); block_id++)
		{
			gradient_matrix[block_id].b += gradient_matrix[block_id].m * gradient_matrix[block_id - 1].b;
			gradient_matrix[block_id].m = gradient_matrix[block_id].m * gradient_matrix[block_id - 1].m;
		}

		static_assert(eqCnt == 2, "Current implementation supports only systems of ODEs with 2 equations and 2 unknown functions.");

		const auto& final_block = *gradient_matrix.rbegin();
		const auto u_0 = resolve_endpoint_corrections(final_block, bcm);

		std::vector<correction> result(gradient_matrix.size() + 1);

		copy(result[0], u_0, gradient_matrix[0].independent_var_id);

		for (int block_id = 0; block_id < gradient_matrix.size(); block_id++)
			copy(result[block_id + 1], gradient_matrix[block_id].m * u_0 + gradient_matrix[block_id].b, gradient_matrix[block_id].independent_var_id);

		return result;
	}

	/// <summary>
	/// Applies correction to the given solution
	/// </summary>
	static void apply_correction(std::vector<mesh_point<R, eqCnt + 1>>& solution, const std::vector<correction>& correction, const R desired_step_size)
	{
		if (solution.size() != correction.size())
			throw std::exception("Invalid input");

		for (int pt_id = 0; pt_id < correction.size(); pt_id++)
			solution[pt_id] -= correction[pt_id].point;

		std::vector<mesh_point<R, eqCnt + 1>> solution_refined;
		solution_refined.reserve(solution.size());

		for (int pt_id = 0; pt_id < solution.size() - 1; pt_id++)
		{
			solution_refined.push_back(solution[pt_id]);
			const auto ind_var_id = correction[pt_id + 1].independent_var_id;
			const auto actual_step = std::abs(solution[pt_id][ind_var_id] - solution[pt_id + 1][ind_var_id]);
			if (actual_step < R(1.5) * desired_step_size)
				continue;

			const int points_to_add = static_cast<int>(actual_step / desired_step_size);
			const auto h = R(1) / (points_to_add + 1);
			const auto position_increment = (solution[pt_id + 1] - solution[pt_id])*h;

			auto temp = solution[pt_id];
			for (int extra_pt_id = 0; extra_pt_id < points_to_add; extra_pt_id++)
			{
				temp += position_increment;
				solution_refined.push_back(temp);
			}
		}

		solution_refined.push_back(*(solution.rbegin()));//put the last point to the refined collection

		solution = solution_refined;
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
	std::vector<R> get_correcion_magnitudes() const
	{
		return correction_magnitudes;
	}

	/// <summary>
	/// The main solving method
	/// </summary>
	std::vector<mesh_point<R, eqCnt + 1>> solve(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const bc_marker<eqCnt>& bcm,
		const R& precision, const R& desired_step, const bool use_reparametrization)
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
			auto g_matrix = construct_extended_gradient_matrix(system, solution, transform_map);
			const auto correction = get_newton_correction(g_matrix, bcm);
			correction_magnitude = std::max_element(correction.begin(), correction.end(),
				[](const auto& a, const auto& b) { return a.magnitude() < b.magnitude(); })->magnitude();

			correction_magnitudes.push_back(correction_magnitude);

			apply_correction(solution, correction, desired_step);

			iter_count++;
		}

		succeeded = (correction_magnitude <= precision);

		return solution;

	}

};