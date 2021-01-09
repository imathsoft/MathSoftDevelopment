#pragma once
#include <vector>
#include <array>
#include <exception>
#include "../ode_system.h"
#include "../../LinearAlgebra/Matrix.h"
#include "../../Utils/AuxUtils.h"

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
	struct strip
	{
		LinAlg::Matrix<R, eqCnt, eqCnt> m;
		LinAlg::Matrix<R, eqCnt, 1> b;
	};

	/// <summary>
	/// Type to hold extended block-diagonal matrix
	/// </summary>
	typedef std::vector<strip> block_matrix;

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
					current_stripe.m[row_id][col_id] = -OneHalf*evaluate_result_current[row_id].grad[col_id];
					temp[row_id][col_id] = -OneHalf*evaluate_result_next[row_id].grad[col_id];
				}

				current_stripe.m[row_id][row_id] -= one_over_h;
				temp[row_id][row_id] += one_over_h;

				current_stripe.b[row_id][0] = divided_diff[row_id] - OneHalf * (evaluate_result_current[row_id].v + evaluate_result_next[row_id].v);
			}

			const auto det = temp.Determinant();

			if (std::abs<R>(det) < 1e3 * std::numeric_limits<R>::epsilon())
				throw std::exception("Singular block");

			const auto temp_inverted = temp.Inverse();

			current_stripe.m = - temp_inverted * current_stripe.m;//inverse sign so that now matrix "m" and vector "b" are on the "same side"
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
	/// Solves system of linear equation wiith respect to the correction of the Newton's method
	/// The correction is returned
	/// Modifies the input gradient matrix 
	/// </summary>
	static std::vector<mesh_point<R, eqCnt + 1>> get_newton_correction(block_matrix& gradient_matrix)
	{
		for (int block_id = 1; block_id < gradient_matrix.size(); block_id++)
		{
			gradient_matrix[block_id].b += gradient_matrix[block_id].m * gradient_matrix[block_id - 1].b;
			gradient_matrix[block_id].m = gradient_matrix[block_id].m * gradient_matrix[block_id - 1].m;
		}

		static_assert(eqCnt == 2, "Current implementation supports only systems of ODEs with 2 equations and 2 unknown functions.");

		//As for the boundary conditions, currently we confine ourselves to consider a particular case, when the consitions are separates,
		//and determine the value of the first unknown function at both endpoints
		//the latter means that the correction to the value of the first unknown function at the endpoints must be "0"
		//which allows us to find the value of the second unknown function at the "left" endpoint as follows:

		const auto& last_block = *gradient_matrix.rbegin();
		LinAlg::Matrix<R, 2, 1> u_0{};
		u_0[1][0] = -last_block.b[0][0] / last_block.m[0][1];

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

	/// <summary>
	/// The main solving method
	/// </summary>
	std::vector<mesh_point<R, eqCnt + 1>> solve(const sys& system, const std::vector<mesh_point<R, eqCnt + 1>>& init_guess, const R& precision)
	{
		auto solution = init_guess;
		correction_magnitudes.clear();

		R correction_magnitude = std::numeric_limits<R>::max();

		int iter_count = 0;
		while (correction_magnitude > precision && iter_count < max_iterations_count)
		{
			auto g_matrix = construst_extended_gradient_matrix(system, solution);
			const auto correction = get_newton_correction(g_matrix);
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