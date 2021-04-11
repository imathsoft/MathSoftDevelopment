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
				[](const auto& args) { return args[1]; }, // v_{0}^{'}(t) = v_{1}(t)
				[lambda](const auto& args) { return lambda * Sinh(lambda * args[0]); }, // v_{1}^{'} = \lambda Sinh(lambda v_{0})
			} },
			{ {0, 0}, {1, 1} });
	}

	/// <summary>
	/// The Troesch's problem with inverted derivative of the unknown function
	/// </summary>
	static simple_bvp<R, 2> Troesch_deriv_inverted(const R lambda)
	{
		return simple_bvp<R, 2>(
			{ {
				[](const auto& args) { return R(1)/args[1]; }, // v_{0}^{'}(t) = 1/v_{1}(t)
				[lambda](const auto& args) { return  - lambda * Sinh(lambda * args[0]) * args[1] * args[1]; }, // v_{1}^{'} = = \lambda Sinh(lambda v_{0}) * v_{1}^{2}
			} },
			{ {0, 0}, {1, 1} });
	}

	/// <summary>
	/// Test boundary value problem #1
	/// Exact solution [4*sin(t), 4*cos(t)]
	/// </summary>
	static simple_bvp<R, 2> BVP_1()
	{
		return simple_bvp<R, 2>(
			{ {
				[](const auto& args) { return args[0]* args[1] + args[1] - R(16) * Cos(args[2])*Sin(args[2]); },
				[](const auto& args) { return args[0] * args[0] + args[1] * args[1] + args[0] - R(16) - R(8) * Sin(args[2]); },
			} },
			{ {0, R(4) * Sin(R(0))}, {1, R(4) * Sin(R(1))} },
			{ [](const auto t) { return R(4) * Sin(t); },
			  [](const auto t) { return R(4) * Cos(t); } });
	}
};
