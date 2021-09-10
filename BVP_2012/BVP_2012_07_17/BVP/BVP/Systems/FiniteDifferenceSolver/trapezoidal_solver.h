#pragma once
#include <vector>
#include <array>
#include <exception>
#include <algorithm>
#include "../ode_system.h"
#include "../../LinearAlgebra/Matrix.h"
#include "../../Utils/AuxUtils.h"
#include "transformation_strategy.h"

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
template <class R, class S = ts_standard, int eqCnt = 2>
class trapezoidal_solver
{
	//TODO: it seems that now when `Matrix` class got general implementation of the matrix inversion functionality,
	//this limitation on the number of equations (`eqCnt == 2`) can be finally removed.
	//Of course, this requires proper coverage of automated tests
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
		transformation_maker<eqCnt> trans_marker{};
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
	/// Logs the transformation-related data to the fiven file on disk
	/// </summary>
	static void LogTransformationStrategy(const block_matrix& matrix, const char* filename)
	{
		std::ofstream file(filename);
		for (const auto& block : matrix)
		{
			file << block.trans_marker.pivot_id << " ";

			for (int eq_id = 0; eq_id < eqCnt; eq_id++)
				file << std::to_string(block.trans_marker.inversion_map[eq_id]) << " ";

			file << std::endl;
		}

		file.close();
	}

	/// <summary>
	/// A data struct to hold evaluation result transformed via a "pivot" transformation
	/// </summary>
	template <class V>
	struct eval_result_transformed_base
	{
		/// <summary>
		/// Transformed result
		/// </summary>
		V res{};

		/// <summary>
		/// Marker defining a transformation
		/// </summary>
		transformation_maker<eqCnt> trans_marker{};

		/// <summary>
		/// Returns values of the righ hand side vector together with the value of the independent variable used in the transformation
		/// </summary>
		mesh_point<R, eqCnt + 1> values_to_mesh_point() const
		{
			mesh_point<R, eqCnt + 1> result;

			for (auto eq_id = 0; eq_id < eqCnt; eq_id++)
				result[eq_id] = res[eq_id].v;

			result[eqCnt] = res.pt[trans_marker.pivot_id];

			return result;
		}
	};

	using eval_result_transformed = eval_result_transformed_base<typename ode_system<R, eqCnt>::eval_result>;
	using eval_result_transformed_minimal = eval_result_transformed_base<typename ode_system<R, eqCnt>::eval_result_minimal>;

	/// <summary>
	/// Applies intependent variable transformation to the given result of the right hand side avaluation according to the given index of independent variable
	/// Affects only values of the right hand side fucntions and not their gradients 
	/// </summary>
	template <class V>
	static typename ode_system<R, eqCnt>::template eval_result_base<V> transform_independent_var_values_only(
		const typename ode_system<R, eqCnt>::template eval_result_base<V>& res, const int& independent_var_id)
	{
		auto result_transformed = res;

		if (independent_var_id != eqCnt)
		{
			const auto denominator_inverted = R(1) / res[independent_var_id].v;
			for (auto eq_id = 0; eq_id < eqCnt; eq_id++)
			{
				if (independent_var_id == eq_id)
					result_transformed[eq_id].v = denominator_inverted;
				else
					result_transformed[eq_id].v *= denominator_inverted;
			}
		}

		return result_transformed;
	}

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
	/// Applies unknowns inversion transformation to the result of the right hand side evaluation according to the given inversion map
	/// Affects only values of the right hand side fucntions and not their gradients
	/// </summary>
	template<class V>
	static typename ode_system<R, eqCnt>::template eval_result_base<V> invert_unknowns_values_only(const typename ode_system<R, eqCnt>::template eval_result_base<V>& res,
		std::array<bool, eqCnt> inversion_map)
	{
		auto result_transformed = res;

		for (auto eq_id = 0; eq_id < eqCnt; eq_id++)
		{
			if (!inversion_map[eq_id])
				continue;

			const auto factor = R(-1) / (res.pt[eq_id] * res.pt[eq_id]);
			result_transformed[eq_id].v *= factor;
		}

		return result_transformed;
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
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res_to_transform, const transformation_maker<eqCnt>& trans_marker)
	{
		const auto res_independent_var_transformed = transform_independent_var(res_to_transform, trans_marker.pivot_id);

		const auto res_unknowns_inverted = invert_unknowns(res_independent_var_transformed, trans_marker.inversion_map);

		return { res_unknowns_inverted, trans_marker };
	}

	/// <summary>
	/// Returns transformations parameters, i.e. independent variable id and map of unknowns to be inverted
	/// </summary>
	template <class V>
	static transformation_maker<eqCnt> get_transform_marker(const typename ode_system<R, eqCnt>::template eval_result_base<V>& res, const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		return S::get_transform_marker<R, V, eqCnt>(res, trans_restrict);
	}

	/// <summary>
	/// Performs transformation of the system evaluation result based on the transformation map and the value of the corresponding derivatives
	/// Transformation map defines what unknownc can be chosen as "pivot" ones (to be swaped with the independent variable via inverting)
	/// </summary>
	static eval_result_transformed transform(const typename ode_system<R, eqCnt>::eval_result& res,	const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		const auto trans_marker = get_transform_marker(res, trans_restrict);

		return transform(res, trans_marker);
	}

	/// <summary>
	/// Subroutine to perform transformation of a minimal result of evaluation (only values of right hand side functions)
	/// </summary>
	static eval_result_transformed_minimal transform_values_only(const typename ode_system<R, eqCnt>::eval_result_minimal& res, const transformation_maker<eqCnt>& trans_marker)
	{
		const auto independent_var_transform_result = transform_independent_var_values_only(res, trans_marker.pivot_id);
		const auto result_final = invert_unknowns_values_only(independent_var_transform_result, trans_marker.inversion_map);
		return { result_final, trans_marker };
	}

	/// <summary>
	/// Subroutine to perform transformation of a minimal result of evaluation (only values of right hand side functions)
	/// </summary>
	static eval_result_transformed_minimal transform_values_only(const typename ode_system<R, eqCnt>::eval_result_minimal& res, const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		const auto trans_marker = get_transform_marker(res, trans_restrict);
		return transform_values_only(res, trans_marker);
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
	/// Calculates "agreement factor" that is aimed to handle the "transition" situatuion, when inversion (the "flip") transformation is applied on the current interval
	/// but was not applied on the previous interval (or vise versa). This transitions state requires us to still work with a not transformed variable 
	/// despite of the fact that equations that we built correspond to the transformed one.
	/// In the other words, we have calculated partial derivative that corresponds to the transformed variable `w = 1/v`, however we need to have 
	/// the partial derivative with respect to `v`. This can be handled by multiplying the partial derivative with respect to `w` by `-1/v^2`.
	/// The latter is actually the "agreement" factor that the function calculates with respect to each variable
	/// </summary>
	static std::array<R, eqCnt> calc_inversion_agreement_factor(const int pivot_id_prev, const std::array<bool, eqCnt>& flip_map_prev,
		const std::array<bool, eqCnt>& flip_map_next,
		const mesh_point<R, eqCnt + 1>& current_pt)
	{
		std::array<R, eqCnt> result{};

		for (int unknown_id = 0; unknown_id < eqCnt; unknown_id++)
		{
			//An obvious case, when the "agreement factor" is trivial (i.e. equal to "1") is when 
			//the "flip" flags for this variable coinside for the current and previious steps
			//Another (more subtle) case, when the agreement factor shoud be trivial is when 
			//the variable with the given index was an independent variable on the previous step.
			//In the latter case, the "flip" flag on the previous step cannot be set by definition (the "flipping" cannot be performed on an independent variable)
			//and the "flip" flg on the current ("next") step, in case it is set, in fact, corresponds to another variable variale and must be ignored
 			if (unknown_id == pivot_id_prev || flip_map_prev[unknown_id] == flip_map_next[unknown_id])
			{
				result[unknown_id] = R(1);
				continue;
			}

			if (flip_map_prev[unknown_id])
				result[unknown_id] = -current_pt[unknown_id] * current_pt[unknown_id];
			else {
				const auto temp = R(1) / current_pt[unknown_id];
				result[unknown_id] = -temp * temp;
			}
		}

		return result;
	}

	/// <summary>
	/// Returns the Jacobian matrix of the vector function F(x_1, x_2, ...,x_{eqCnt+1}) = ([y_1, y_2, ...,y_{eqCnt + 1}] - [x_1, x_2, ...,x_{eqCnt + 1}])/(y_{independent_var_id} - x_{independent_var_id})
	/// multiplied by (y_{independent_var_id} - x_{independent_var_id}) (i.e., the Jacobian matrix is multiplied by the mentioned factor),
	/// where X and Y vectors are defined with their "divided difference"
	/// </summary>
	static LinAlg::Matrix<R, eqCnt, eqCnt + 1> CalcJacobianOfDividedDifference(const mesh_point<R, eqCnt + 1>& divided_difference, const int independent_var_id)
	{
		LinAlg::Matrix<R, eqCnt, eqCnt + 1> result{};

		for (int row_id = 0; row_id < eqCnt; row_id++)
		{
			//Since the variable "swapping" mechanism is involved, the row index is not necessarly corresponds to the "variaable" index
			//however, knowing the index of the independent variable, we can deduce the actual variable index according to the following formula
			const auto var_id = row_id != independent_var_id ? row_id : eqCnt;

			for (int col_id = 0; col_id < eqCnt + 1; col_id++)
			{
				//The case when column index is equal to the variable index, 
				//which means that we need to take derivative of (y_{col_id} - x_{col_id})/(y_{independent_var_id} - x_{independent_var_id}) with respect to x_{col_id}
				//and mulltily the resutl by (y_{independent_var_id} - x_{independent_var_id}) (see the summary of the method)
				if (col_id == var_id)
					result[row_id][col_id] = - R(1);
				//The case when we need to take derivative of (y_{col_id} - x_{col_id})/(y_{independent_var_id} - x_{independent_var_id})  
				//but this time with respect to x_{independent_var_id}) and again multiply the result by (y_{independent_var_id} - x_{independent_var_id})
				else if (col_id == independent_var_id)
					result[row_id][col_id] = divided_difference[var_id];
				else
			    //The case when the expression that we need to differentiate is not explicitely depentent of x_{col_id}, so the derivative if zero
					result[row_id][col_id] = R(0);
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
	static block_matrix construct_extended_gradient_matrix(const sys& system,
		const std::vector<mesh_point<R, eqCnt + 1>>& init_guess,
		const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		if (init_guess.size() <= 1)
			throw std::exception("Invalid input");

		block_matrix result(init_guess.size() - 1);

		auto res_prev = system.evaluate(init_guess[0]);
		transformation_maker<eqCnt> trans_marker_prev{ -1 };
		for (auto pt_id = 0; pt_id < init_guess.size() - 1; pt_id++)
		{
			const auto res_prev_transformed = transform(res_prev, trans_restrict);
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
			const auto divided_diff_J = CalcJacobianOfDividedDifference(divided_diff, independent_var_id_next);


			const auto agreement_factor = calc_inversion_agreement_factor(independent_var_id_prev, inversion_map_prev, inversion_map_next, init_guess[pt_id]);

			auto& current_stripe = result[pt_id];
			current_stripe.trans_marker = trans_marker_next;
			LinAlg::Matrix<R, eqCnt, eqCnt> temp;

			for (int row_id = 0; row_id < eqCnt; row_id++)
			{
				for (int col_id = 0; col_id < eqCnt; col_id++)
				{
					const auto var_id_prev = col_id != independent_var_id_prev ? col_id : eqCnt;
					current_stripe.m[row_id][col_id] = -OneHalf * agreement_factor[col_id] * res_prev_transformed.res[row_id].grad[var_id_prev] +
						agreement_factor[col_id] * divided_diff_J[row_id][var_id_prev] * one_over_h;
					temp[row_id][col_id] = -OneHalf * res_next_transformed.res[row_id].grad[col_id != independent_var_id_next ? col_id : eqCnt];
				}

				temp[row_id][row_id] += one_over_h;

				const int actual_var_id = row_id != independent_var_id_next ? row_id : eqCnt;
				current_stripe.b[row_id][0] = divided_diff[actual_var_id] - OneHalf * (res_prev_transformed.res[row_id].v + res_next_transformed.res[row_id].v);
			}

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
		transformation_maker<eqCnt> trans_marker{};

		/// <summary>
		/// Correction magnitude
		/// </summary>
		R magnitude() const
		{
			return point.max_abs();
		}

		/// <summary>
		/// Returns string representation of the current instance
		/// </summary>
		std::string to_string() const
		{
			return point.to_string() + " " + trans_marker.to_string();
		}
	};

	/// <summary>
	/// A helper method to copy result from vector-column to the mesh point
	/// </summary>
	static void copy(correction& dest, const LinAlg::Matrix<R, eqCnt, 1>& source, const transformation_maker<eqCnt>& trans_marker)
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

		const auto solution = - temp.Inverse() * final_block.b;

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
	/// All the data needed to perform "smart" mesh refinement
	/// </summary>
	struct refinement_data
	{
		/// <summary>
		/// The right hand side ('rhs') part of the system evaluated at some point
		/// </summary>
		typename ode_system<R, eqCnt>::eval_result_minimal rhs;

		/// <summary>
		/// The transformed version of the right hand side ('rhs') part of the system evaluated at some point
		/// </summary>
		eval_result_transformed_minimal rhs_transformed;

		/// <summary>
		/// Returns index of the independent variable that was chosen when calculating the transformed right hand side values
		/// </summary>
		int independent_var_id() const
		{
			return rhs_transformed.trans_marker.pivot_id;
		}

		/// <summary>
		/// Returns transformend right hand side evaluation result as a mesh point of the corresponding dimension
		/// </summary>
		mesh_point<R, eqCnt + 1> transformed_to_point() const
		{
			return rhs_transformed.values_to_mesh_point();
		}

		/// <summary>
		/// Transforms the right hand side evaluation result according to the given transformation marker and returns the result 
		/// in a form of a mesh point of the corresponding dimension
		/// For the optomization purposes, the inpur marker is compared to the local transformation marker and if they are equal
		/// the existing transformed result is used (no additional calculations)
		/// </summary>
		mesh_point<R, eqCnt + 1> transformed_to_point(const transformation_maker<eqCnt>& trans_marker) const
		{
			if (trans_marker == rhs_transformed.trans_marker)
				return rhs_transformed.values_to_mesh_point();

			return transform_values_only(rhs, trans_marker).values_to_mesh_point();
		}
	};

	/// <summary>
	/// Returns the refinement data structure calculated for the given mesh point
	/// </summary>
	static refinement_data calc_refinement_data(const mesh_point<R, eqCnt + 1>& pt,
		const sys& system, const transform_restrictions<R, eqCnt>& trans_restrict)
	{
		auto eval_res = system.evaluate_minimal(pt);
		auto eval_res_transformed = transform_values_only(eval_res, trans_restrict);

		return { eval_res, eval_res_transformed };
	}

	/// <summary>
	/// Performs clean up and refinement of the solution
	/// </summary>
	static bool cleanup_and_refine_solution(const sys& system, std::vector<mesh_point<R, eqCnt + 1>>& solution,
		const transform_restrictions<R, eqCnt>& trans_restrict, const R& desired_step_size, const R& second_deriv_threshold,
		const R& min_step_size, const bool optimize_step_size, const R remove_threshold = R(0.1))
	{
		if (solution.size() <= 2)
			return false; //There is literally nothing to refine

		const auto  solution_input = solution;

		const auto argument_min = solution[0][eqCnt];
		const auto argument_max = solution[solution.size() - 1][eqCnt];
		const auto independent_var_scale_factor = second_deriv_threshold / desired_step_size;

		solution.clear();
		auto pt_prev = solution_input[0];
		auto refine_data_prev = calc_refinement_data(pt_prev, system, trans_restrict);
		
		auto pt_next = solution_input[1];
		auto refine_data_next = calc_refinement_data(pt_next, system, trans_restrict);

		solution.push_back(pt_prev);
		bool solution_altered = false;

		int pt_id = 1;
		while (true)
		{
			const auto independent_var_id = refine_data_prev.independent_var_id();
			const auto prev_trans_pt = refine_data_prev.transformed_to_point();
			const auto next_trans_pt = refine_data_next.transformed_to_point(refine_data_prev.rhs_transformed.trans_marker);
			auto diff = prev_trans_pt - next_trans_pt;
			diff[eqCnt] *= independent_var_scale_factor* (optimize_step_size ? auxutils::Sqrt(R(0.5) * (prev_trans_pt.max_abs(eqCnt) + next_trans_pt.max_abs(eqCnt))) : R(1));

			const auto actual_step_size = diff.max_abs();

			//Check for "zigzags"
			if ((pt_next[eqCnt] < argument_min || pt_next[eqCnt] > argument_max || 
				pt_next[eqCnt] <= pt_prev[eqCnt] ||
				actual_step_size < second_deriv_threshold * remove_threshold ||
				auxutils::Abs(pt_next[independent_var_id] - pt_prev[independent_var_id]) < min_step_size) && // minimal allowed step size check
				pt_id < (solution_input.size() - 1)) 
			{
				solution_altered = true;
			} else 
			{
				const int points_to_add = static_cast<int>(actual_step_size / second_deriv_threshold);
				const auto h = R(1) / (points_to_add + 1);
				const auto position_increment = (pt_next - pt_prev) * h;

				auto temp = pt_prev;
				for (int extra_pt_id = 0; extra_pt_id < points_to_add + 1; extra_pt_id++)
				{
					if (extra_pt_id != points_to_add)
						temp += position_increment;
					else
						temp = pt_next;

					const auto ind_var_id = refine_data_prev.independent_var_id();

					if (auxutils::Abs(temp[ind_var_id] - pt_prev[ind_var_id]) < min_step_size)
						continue;

					pt_prev = temp;
					refine_data_prev = calc_refinement_data(pt_prev, system, trans_restrict);
					solution.push_back(pt_prev);
					solution_altered = solution_altered || (extra_pt_id != points_to_add);//we just added a new point so the solution was altered
				}
			}

			pt_id++;

			if (pt_id == solution_input.size())
				break;
			
			pt_next = solution_input[pt_id];
			refine_data_next = calc_refinement_data(pt_next, system, trans_restrict);
		}

		return solution_altered;
	}

	/// <summary>
	/// A diagnostics tool
	/// Evaluates the right hand side part of the given system at each point of the given collection of points 
	/// thransdorms the result of avaluation according to the given tansformation restrictions 
	/// and saves it to the file with the given name (in the human readable text format)
	/// </summary>
	/// <param name="system">The given systen of ODE</param>
	/// <param name="points">The given collection of points</param>
	/// <param name="trans_restrict">The transformation restrictions</param>
	void transform_and_save_evaluation_result(const sys& system, std::vector<mesh_point<R, eqCnt + 1>>& points,
		const transform_restrictions<R, eqCnt>& trans_restrict, const char* file_name)
	{
		std::ofstream file;
		file.open(file_name, std::ofstream::out);

		for (const auto& pt : points)
		{
			const auto eval_res = system.evaluate_minimal(pt);
			const auto eval_res_transformed = transform_values_only(eval_res, trans_restrict);
			file << eval_res_transformed.values_to_mesh_point().to_string() << eval_res_transformed.trans_marker.pivot_id << " ";

			for (int eqId = 0; eqId < eqCnt; eqId++)
			{
				file << std::to_string(eval_res_transformed.trans_marker.inversion_map[eqId]).c_str() << " ";
			}

			file << pt.to_string() << std::endl;
		}

		file.close();
	}

	/// <summary>
	/// Applies correction to the given solution
	/// Returns true if actual refinement was applied to the corrected solution
	/// </summary>
	static bool apply_correction(const sys& system, std::vector<mesh_point<R, eqCnt + 1>>& solution, const std::vector<correction>& correction)
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
	/// A setting defining whether to do pre-refinement before running the Newton's procedure
	/// </summary>
	bool DoPreRefinement = false;

	/// <summary>
	/// A setting defining whether the step size should be optimized by using some step selection strategies when refining the mesh
	/// </summary>
	bool OptimizeStepSize = false;

	/// <summary>
	/// A value that defines the minimal ratio between the "actual" and "desired" step size of the mesh
	/// If the ratio is smaller the corresponding points will be removed during the refinement procedure
	/// </summary>
	R RefinementRemoveFactor = R(0.1);

	/// <summary>
	/// The main solving method
	/// </summary>
	std::vector<mesh_point<R, eqCnt + 1>> solve(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const bc_marker<eqCnt>& bcm,
		const R& precision, const R& desired_step, const bool use_swap_transform, const bool use_flip_transform, const R& derivative_threshold = R(1),
		const R& second_deriv_refinement_threshold = R(20),
		const R& min_h_threshold = R(2e-5),
		std::optional<std::array<bool, eqCnt>> swap_transform_map = std::nullopt,
		std::optional<std::array<bool, eqCnt>> flip_transform_map = std::nullopt)
	{
		if (!bcm.is_valid())
			throw std::exception("Invalid boundary condition marker");

		std::array<bool, eqCnt> swap_map{};

		if (use_swap_transform)
		{
			if (!swap_transform_map.has_value())
				for (int var_id = 0; var_id < eqCnt; var_id++)
					swap_map[var_id] = bcm.left_marker[var_id] && bcm.right_marker[var_id];
			else
				swap_map = swap_transform_map.value();
		}

		std::array<bool, eqCnt> flip_map{};

		if (use_flip_transform)
		{
			if (!flip_transform_map.has_value())
				for (int var_id = 0; var_id < eqCnt; var_id++)
					flip_map[var_id] = true;
			else
				flip_map = flip_transform_map.value();
		}


		const transform_restrictions<R, eqCnt> trans_restrict{ swap_map , flip_map, derivative_threshold };

		auto solution = init_guess;
		correction_magnitudes.clear();

		R correction_magnitude = std::numeric_limits<R>::max();

		if (DoPreRefinement)
			cleanup_and_refine_solution(system, solution, trans_restrict, desired_step,
				second_deriv_refinement_threshold, min_h_threshold, OptimizeStepSize, RefinementRemoveFactor*R(10));

		int iter_count = 0;
		while (correction_magnitude > precision && iter_count < max_iterations_count)
		{
			auto g_matrix = construct_extended_gradient_matrix(system, solution, trans_restrict);
			const auto correction = get_newton_correction(g_matrix, bcm);
			correction_magnitude = std::max_element(correction.begin(), correction.end(),
				[](const auto& a, const auto& b) { return a.magnitude() < b.magnitude(); })->magnitude();

			apply_correction(system, solution, correction);

			const auto refinement_applied = (use_swap_transform) ?
				cleanup_and_refine_solution(system, solution, trans_restrict, desired_step,
					second_deriv_refinement_threshold, min_h_threshold, OptimizeStepSize, RefinementRemoveFactor) : false;

			correction_magnitudes.emplace_back( ConvergenceInfo<R>{ correction_magnitude, refinement_applied });

			iter_count++;
		}

		succeeded = (correction_magnitude <= precision);

		return solution;
	}
};