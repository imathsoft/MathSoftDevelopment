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
	/// Constructs the extended gradient matrix, which is essentially
	/// [F'(s_{i}) | F(s_{i})],
	/// where F(s) = 0 is the system of nonliner equations that we get applying the trapezoidal scheme to the given system of ODEs 
	/// s_{i} is the "initial guess"
	/// </summary>
	static block_matrix construst_extended_gradient_matrix(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess)
	{
		if (init_guess.size() <= 1)
			throw std::exception("Invalid input");

		block_matrix result(init_guess.size() - 1);

		auto evaluate_result_current = system.Evaluate(init_guess[0]);
		for (int pt_id = 0; pt_id < init_guess.size() - 1; pt_id++)
		{
			const auto evaluate_result_next = system.Evaluate(init_guess[pt_id + 1]);

			const auto pt_diff = init_guess[pt_id + 1] - init_guess[pt_id];
			const auto one_over_h = R(1) / pt_diff.pt[eqCnt];
			const auto divided_diff = one_over_h * pt_diff;

			auto& current_stripe = result[pt_id];
			LinAlg::Matrix<R, eqCnt, eqCnt> temp;

			for (int row_id = 0; row_id < eqCnt; row_id++)
			{
				for (int col_id = 0; col_id < eqCnt; col_id++)
				{
					current_stripe.m[row_id][col_id] = -OneHalf * evaluate_result_current[row_id].grad[col_id];
					temp[row_id][col_id] = -OneHalf * evaluate_result_next[row_id].grad[col_id];
				}

				current_stripe.m[row_id][row_id] -= one_over_h;
				temp[row_id][row_id] += one_over_h;

				current_stripe.b[row_id][0] = divided_diff[row_id] - OneHalf * (evaluate_result_current[row_id].v + evaluate_result_next[row_id].v);
			}

			const auto det = temp.Determinant();

			if (std::abs<R>(det) < 1e3 * std::numeric_limits<R>::epsilon())
				throw std::exception("Singular block");

			const auto temp_inverted = temp.Inverse();

			current_stripe.m = -temp_inverted * current_stripe.m;//inverse sign so that now matrix "m" and vector "b" are on the "same side"
			current_stripe.b = temp_inverted * current_stripe.b;

			evaluate_result_current = evaluate_result_next;
		}

		return result;
	}

	/// <summary>
	/// A helper method to copy result from vector-column to the mesh point
	/// </summary>
	static void copy(mesh_point<R, eqCnt + 1>& dest, const LinAlg::Matrix<R, eqCnt, 1>& source)
	{
		for (int i = 0; i < eqCnt; i++)
			dest[i] = source[i][0];
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
	static std::vector<mesh_point<R, eqCnt + 1>> get_newton_correction(block_matrix& gradient_matrix, const bc_marker<eqCnt>& bcm)
	{
		for (int block_id = 1; block_id < gradient_matrix.size(); block_id++)
		{
			gradient_matrix[block_id].b += gradient_matrix[block_id].m * gradient_matrix[block_id - 1].b;
			gradient_matrix[block_id].m = gradient_matrix[block_id].m * gradient_matrix[block_id - 1].m;
		}

		static_assert(eqCnt == 2, "Current implementation supports only systems of ODEs with 2 equations and 2 unknown functions.");

		const auto& final_block = *gradient_matrix.rbegin();
		const auto u_0 = resolve_endpoint_corrections(final_block, bcm);

		std::vector<mesh_point<R, eqCnt + 1>> result(gradient_matrix.size() + 1);

		copy(result[0], u_0);

		for (int block_id = 0; block_id < gradient_matrix.size(); block_id++)
			copy(result[block_id + 1], gradient_matrix[block_id].m * u_0 + gradient_matrix[block_id].b);

		return result;
	}

	/// <summary>
	/// Applies correction to the given solution
	/// </summary>
	static void apply_correction(std::vector<mesh_point<R, eqCnt + 1>>& solution, const std::vector<mesh_point<R, eqCnt + 1>>& correction)
	{
		if (solution.size() != correction.size())
			throw std::exception("Invalid input");

		for (int pt_id = 0; pt_id < correction.size(); pt_id++)
			solution[pt_id] -= correction[pt_id];
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
	std::vector<mesh_point<R, eqCnt + 1>> solve(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const bc_marker<eqCnt>& bcm, const R& precision)
	{
		if (!bcm.is_valid())
			throw std::exception("Invalid boundary condition marker");

		auto solution = init_guess;
		correction_magnitudes.clear();

		R correction_magnitude = std::numeric_limits<R>::max();

		int iter_count = 0;
		while (correction_magnitude > precision && iter_count < max_iterations_count)
		{
			auto g_matrix = construst_extended_gradient_matrix(system, solution);
			const auto correction = get_newton_correction(g_matrix, bcm);
			correction_magnitude = std::max_element(correction.begin(), correction.end(),
				[](const auto& a, const auto& b) { return a.max_abs() < b.max_abs(); })->max_abs();

			correction_magnitudes.push_back(correction_magnitude);

			apply_correction(solution, correction);
			iter_count++;
		}

		succeeded = (correction_magnitude <= precision);

		return solution;

	}

};