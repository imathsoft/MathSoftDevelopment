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
				ode_system<R, 2>::create_func([](const auto& args) { return args[1]; }), // v_{0}^{'}(t) = v_{1}(t)
				ode_system<R, 2>::create_func([lambda](const auto& args) { return lambda * Sinh(lambda * args[0]); }), // v_{1}^{'} = \lambda Sinh(lambda v_{0})
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
				// v_{0}^{'}(t) = 1/v_{1}(t)
				ode_system<R, 2>::create_func([](const auto& args) { return R(1)/args[1]; }),
				// v_{1}^{'} = = \lambda Sinh(lambda v_{0}) * v_{1}^{2}
				ode_system<R, 2>::create_func([lambda](const auto& args) { return  - lambda * Sinh(lambda * args[0]) * args[1] * args[1]; }),
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
				ode_system<R, 2>::create_func([](const auto& args) { return args[0]* args[1] + args[1] - R(16) * Cos(args[2])*Sin(args[2]); }),
				ode_system<R, 2>::create_func([](const auto& args) { return args[0] * args[0] + args[1] * args[1] + args[0] - R(16) - R(8) * Sin(args[2]); }),
			} },
			{ {0, R(4) * Sin(R(0))}, {1, R(4) * Sin(R(1))} },
			{ [](const auto t) { return R(4) * Sin(t); },
			  [](const auto t) { return R(4) * Cos(t); } });
	}

	/// <summary>
	/// Test boundary value problem #1
	/// Exact solution [t^{3}, t^{2} + 1]
	/// </summary>
	static simple_bvp<R, 2> BVP_2()
	{
		return simple_bvp<R, 2>(
			{ {
				ode_system<R, 2>::create_func([](const auto& args) { return R(3) * args[1] - R(3) + Sin(args[0] + args[1]) - Sin(args[2] * args[2] * (args[2] + R(1)) + R(1)); }),
				ode_system<R, 2>::create_func([](const auto& args) { return R(2) * (args[0] + args[2]) / args[1] + Cos(args[0] / args[1]) - Cos(args[2] * args[2] * args[2] / (args[2] * args[2] + R(1))); }),
			} },
			{ {0, R(0) }, {1, R(1) } },
			{ [](const auto t) { return t * t * t; },
			  [](const auto t) { return t * t + R(1); } });
	}

	/// <summary>
	/// Test boundary value problem bvp_t_30
	/// </summary>
	static simple_bvp<R, 2> BVP_T30(const R lambda)
	{
		return simple_bvp<R, 2>(
			{ {
				ode_system<R, 2>::create_func([lambda](const auto& args) { return args[1]; }),
				ode_system<R, 2>::create_func([lambda](const auto& args) { return lambda * args[0] * (R(1) - args[1]); }),
			} },
			{ {0, -R(7)/R(6) }, {1, R(3)/R(2) } },
			{ [](const auto t) { return std::numeric_limits<R>::quiet_NaN(); },
			  [](const auto t) { return std::numeric_limits<R>::quiet_NaN(); } });
	}

	/// <summary>
	/// Test boundary value problem bvp_t_30 (modified version with derivatives tending to zero at both boundary points as lambda tends to infinity)
	/// </summary>
	static simple_bvp<R, 2> BVP_T30_1(const R lambda)
	{
		return simple_bvp<R, 2>(
			{ {
				ode_system<R, 2>::create_func([lambda](const auto& args) { return args[1]; }),
				ode_system<R, 2>::create_func([lambda](const auto& args) { return - lambda * (args[0] + args[2]) * args[1]; }),
			} },
			{ {0, -R(7) / R(6) }, {1, R(3) / R(2) - R(1) } },
			{ [](const auto t) { return std::numeric_limits<R>::quiet_NaN(); },
			  [](const auto t) { return std::numeric_limits<R>::quiet_NaN(); } });
	}


	/// <summary>
	/// Test boundary value problem bvp_t20, possesing a corner layer at point "alpha"
	/// </summary>
	static simple_bvp<R, 2> BVP_T20(const R lambda, const R alpha = R(0.5))
	{
		return simple_bvp<R, 2>(
			{ {
				ode_system<R, 2>::create_func([lambda](const auto& args) { return args[1]; }),
				ode_system<R, 2>::create_func([lambda](const auto& args) { return -lambda * (args[1] * args[1] - 1); }),
			} },
			{ {0, R(1) + auxutils::Log(auxutils::Cosh(-lambda*alpha))/lambda}, {1, R(1) + auxutils::Log(auxutils::Cosh(lambda*(R(1) - alpha)))/lambda } },
			{ [lambda, alpha](const auto t) { return R(1) + auxutils::Log(auxutils::Cosh(lambda * (t - alpha)))/lambda; },
			  [lambda, alpha](const auto t) { return auxutils::Tanh(lambda * (t - alpha)); } });
	}


};
