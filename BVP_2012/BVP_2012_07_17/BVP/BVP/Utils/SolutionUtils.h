#pragma once
#include "../FunctionApproximation/InitialCondition.h"
#include <vector>
#include <numeric>
#include <limits>
#include "AuxUtils.h"
#include <boost\heap\pairing_heap.hpp>

namespace solutionutils
{
	///A method to calculate length of the solution (Deriv == false) or its derivatve (Deriv == true)
	template <bool Deriv, class T>
	T length(const std::vector<InitCondition<T>>& polyline)
	{
		T result = T(0);

		for (int pt_id = 0; pt_id < polyline.size() - 1; pt_id++)
		{
			result += auxutils::Sqrt((polyline[pt_id + 1] - polyline[pt_id]).NormSquaredPartial<Deriv>());
		}

		return result;
	}

	template <bool Deriv, class T>
	T calculate_cost(const std::vector<InitCondition<T>>& polyline, const std::vector<std::array<int, 2>>& neighbors, const int pt_id)
	{
		const auto a = auxutils::Sqrt((polyline[pt_id] - polyline[neighbors[pt_id][0]]).NormSquaredPartial<Deriv>());
		const auto b = auxutils::Sqrt((polyline[neighbors[pt_id][1]] - polyline[pt_id]).NormSquaredPartial<Deriv>());
		const auto c = auxutils::Sqrt((polyline[neighbors[pt_id][1]] - polyline[neighbors[pt_id][0]]).NormSquaredPartial<Deriv>());

		return  a + b - c;
	}

	template <class T>
	struct comparer
	{
		bool operator () (const std::pair<int, T>& a, const std::pair<int, T>& b) const
		{
			return a.second > b.second;
		}
	};

	//A method to simplify the solution (Deriv == false) or its derivative (Deriv = true), minimizing the change of the curve's length
	template <bool Deriv, class T>
	std::vector<InitCondition<T>> simplify_polyline(const std::vector<InitCondition<T>>& solution, const float percentagePointsToKeep)
	{
		if (solution.size() <= 2)
			return solution;

		const int pointsToKeep = std::max<int>(2, solution.size() * percentagePointsToKeep - 2);

		if (pointsToKeep >= solution.size())
			return solution;
		 
		std::vector<std::array<int, 2>> neighbors(solution.size());
		boost::heap::pairing_heap<std::pair<int,T>, boost::heap::compare<comparer<T>>> pair_heap;
		std::vector<boost::heap::pairing_heap<std::pair<int,T>, boost::heap::compare<comparer<T>>>::handle_type> heap_handles(solution.size()); 

		neighbors[0][0] = -1;
		neighbors[0][1] = 1;
		neighbors[neighbors.size() - 1][0] = neighbors.size()-2;
		neighbors[neighbors.size() - 1][1] = -1;

		for (int pt_id = 1; pt_id < solution.size() - 1; pt_id++)
		{
			neighbors[pt_id][0] = pt_id - 1;
			neighbors[pt_id][1] = pt_id + 1;
			const auto cost = calculate_cost<Deriv>(solution, neighbors, pt_id);
			heap_handles[pt_id] = pair_heap.push(std::pair<int, T>(pt_id, cost));
		}

		while (pair_heap.size() > pointsToKeep)
		{
		    const int pt_id = pair_heap.top().first;
			pair_heap.pop(); 

			neighbors[neighbors[pt_id][0]][1] = neighbors[pt_id][1];
			neighbors[neighbors[pt_id][1]][0] = neighbors[pt_id][0];

			if (neighbors[pt_id][0] != 0) //do not update the heap if the left neighbor is the first point in the polyline
			{
				const auto cost = calculate_cost<Deriv>(solution, neighbors, neighbors[pt_id][0]);
				pair_heap.update(heap_handles[neighbors[pt_id][0]], std::pair<int, T>(neighbors[pt_id][0], cost));
			}

			if (neighbors[pt_id][1] != solution.size() - 1)//do not update the heap if the left neighbor is the last point in the polyline
			{
				const auto cost = calculate_cost<Deriv>(solution, neighbors, neighbors[pt_id][1]);
				pair_heap.update(heap_handles[neighbors[pt_id][1]], std::pair<int, T>(neighbors[pt_id][1], cost));
			}

			neighbors[pt_id][1] = -1;
			neighbors[pt_id][0] = -1;
		}

		std::vector<InitCondition<T>> result;

		result.push_back(solution[0]);

		int next_pt_id = neighbors[0][1];

		while (next_pt_id != -1)
		{
			result.push_back(solution[next_pt_id]);
			if (next_pt_id >= neighbors[next_pt_id][1])
				int i = 7;

			next_pt_id = neighbors[next_pt_id][1];
		}

		return result;
	}
}