#pragma once

#include "trapezoidal_solver.h"
#include "transformation_strategy.h"
#include "../../Utils/SolutionUtils.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

/// <summary>
/// Diagnostics functionslity for the trapezoidal_solver class
/// The class below is supposed to be a friend to the trapezoidal_solver class
/// to det access to its inner data
/// </summary>
class ts_diagnostics
{
private:

	/// <summary>
	/// Simplifies the given curve by removing its points while keeping the length as close to the initial one as possible
	/// </summary>
	template <class R, int varCnt>
	static std::vector<mesh_point<R, varCnt>> simplify(const std::vector<mesh_point<R, varCnt>>& curve, const int pt_cnt_desired)
	{
		if (curve.size() < pt_cnt_desired)
			return curve; // we can't increase number of points, so return

		const auto percentage_points_to_keep = pt_cnt_desired * (1.0 / curve.size());
		return solutionutils::simplify_polyline<false>(curve, percentage_points_to_keep);
	}

	/// <summary>
	/// Saves the given table of data into the given file in a csv-like format
	/// </summary>
	template <class R>
	static void save_table(const std::string& file_name, const std::vector<std::vector<mesh_point<R, 2>>>& table)
	{
		std::ofstream file(file_name);

		const auto max_rows_cnt = (*std::max_element(table.begin(), table.end(),
			[](const auto& e1, const auto& e2) { return e1.size() < e2.size(); })).size();

		for (int row_id = 0; row_id < max_rows_cnt; row_id++)
		{
			for (int col_id = 0; col_id < table.size(); col_id++)
			{
				if (table[col_id].size() <= row_id)
				{
					file << " , ,";
				} else
					file << auxutils::ToString(table[col_id][row_id][0]) << ", " << auxutils::ToString(table[col_id][row_id][1]) << ", ";
			}
			file << std::endl;
		}

		file.close();
	}

public:
	/// <summary>
	/// Saves norms of the right hand side vector of the given system of ODEs
	/// evaluated at each poit of the given mesh into a csv-like table
	/// Additionally to that, for each mesh point, a transformed right hand
	/// side vector gets evaluated (as suggested by the solver) and aved
	/// </summary>
	template <class R, class S, int eqCnt>
	static void save_rhs_norms_on_mesh(const std::string& file_name, const ode_system<R, eqCnt>& system,
		const std::vector<mesh_point<R, eqCnt + 1>>& mesh, const S& solver, const int max_pt_count = 100)
	{
		std::vector<std::vector<mesh_point<R, 2>>> table(eqCnt + 2);

		for (const auto& pt : mesh)
		{
			const auto eval_res = system.evaluate_minimal(pt);
			const auto eval_res_transformed = S::transform_values_only(eval_res, solver._trans_restrict);
			const auto eval_res_pt = eval_res.values_to_mesh_point();
			table[0].push_back({ pt[eqCnt], eval_res_pt.max_abs(eqCnt) });
			const auto eval_res_transformed_pt = eval_res_transformed.values_to_mesh_point();
			const auto pivot_id = eval_res_transformed.trans_marker.pivot_id;
			table[pivot_id + 1].push_back({ eval_res_transformed_pt[eqCnt], eval_res_transformed_pt.max_abs(eqCnt) });
		}

		//remove empty columns
		const auto new_end = std::remove_if(table.begin(), table.end(), [](const auto& x) { return x.size() == 0; });
		table.resize(std::distance(table.begin(), new_end));

		for (auto& curve : table)
			curve = simplify(curve, max_pt_count);

		save_table(file_name, table);
	}

	/// <summary>
	/// Reduces the number of poitns in the given mesh to ensure the
	/// given maximal number of points and saves the result to disk in text format
	/// </summary>
	template <class R, int varCnt>
	static void reduce_and_save_mesh(const std::string& file_name, const std::vector<mesh_point<R, varCnt>>& mesh, const int max_pt_count = 100)
	{
		auto mesh_reduced = simplify(mesh, max_pt_count);
		auxutils::SaveToTextFile(mesh_reduced, file_name.c_str());
	}
};
