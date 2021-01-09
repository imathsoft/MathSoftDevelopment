#pragma once
#include "simple_bvp.h"

/// <summary>
/// Factory that creates different boundary value problems for systems of differential equations
/// </summary>
template<class R>
class bvp_sys_factory
{
public:

	/// <summary>
	/// Returns the Troesch's problem in the "system" formulation
	/// </summary>
	static simple_bvp<R, 2> Troesch(const R lambda)
	{
		return simple_bvp<R, 2>(
			{ {
				[](const auto& args) { return args[1]; }, // v_{0}^{'}(t) = v_{i}(t)
				[lambda](const auto& args) { return lambda * Sinh(lambda * args[0]); }, // v_{1}^{'} = \lambda Sinh(lambda v_{0})
			} },
			{ {0, 0}, {1, 1} });
	}

};
