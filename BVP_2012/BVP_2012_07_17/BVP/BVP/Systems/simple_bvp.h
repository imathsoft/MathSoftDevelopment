#pragma once
#include "ode_system.h"

/// <summary>
/// A structure to represet an argument-value pair of a scalar function of scalar argument
/// </summary>
template <class R>
struct func_value
{
	/// <summary>
	/// Function argument
	/// </summary>
	R x;

	/// <summary>
	/// Function value
	/// </summary>
	R v;
};

/// <summary>
/// A structure to holde the "simplest" boundary conditions
/// </summary>
template <class R>
struct bnd_cond_simple
{
	/// <summary>
	/// Boundary condition at the "left" argument point
	/// </summary>
	func_value<R> left;

	/// <summary>
	/// Boundary condition ar the "right" argument point
	/// </summary>
	func_value<R> right;
};

/// <summary>
/// A class representing a "simple" boundary value problem
/// for system of ordinary differential equations of the first order
/// The "simplicity" comes from the fact that boundary conditions are decoupled and imposed only on the value of the 
/// first unknown function
/// </summary>
template <class R, int eqCnt>
class simple_bvp
{
private:
	/// <summary>
	/// Boundary conditions
	/// </summary>
	bnd_cond_simple<R> boundary_conditions;

	/// <summary>
	/// System of differential equations
	/// </summary>
	ode_system<R, eqCnt> system;

public:

	/// <summary>
	/// Constructor
	/// </summary>
	simple_bvp(ode_system<R, eqCnt> sys, const bnd_cond_simple<R> bc) : system{ sys }, boundary_conditions{bc}
	{}

	/// <summary>
	/// Getter for the boundary conditions field
	/// </summary>
	bnd_cond_simple<R> get_bc() const
	{
		return boundary_conditions;
	}

	/// <summary>
	/// Getter for the "system" field
	/// </summary>
	const ode_system<R, eqCnt>& get_system() const
	{
		return system;
	}
};
