#ifndef GUARD_X_CANNON
#define GUARD_X_CANNON

#include "XCannonAbstract.h"
#include "../FunctionApproximation/X_Function.h"

template <class T>
class XCannon : public XCannonAbstract<T>
{
protected:

	/// <summary>
	/// Gets the next knot.
	/// </summary>
	/// <param name="prevKnot">The previous knot.</param>
	/// <param name="argFinish">The argument finish.</param>
	/// <returns></returns>
	virtual inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish) override
	{
		T A = _problem->GetACoeff(prevKnot.Derivative, prevKnot.Value, prevKnot.Argument);
		T B = _problem->GetBCoeff(prevKnot.Derivative, prevKnot.Value, prevKnot.Argument);
		T C = prevKnot.Derivative;
		T D = prevKnot.Value;
		T E = _problem->GetECoeff(prevKnot.Argument);
		T F = _problem->GetFCoeff(prevKnot.Argument);

		T hOpt = _problem->get_optimal_step_straight(A, B, _h);

		T H = sgn(_h) * min(hOpt, abs(argFinish - prevKnot.Argument));

		InitCondition<T> knot = _problem->step_streight(A, B, C, D, E, F, H, _precision);
		knot.Argument = prevKnot.Argument + H;

		return knot;
	};

public:	
	/// <summary>
	/// Initializes a new instance of the <see cref="XCannon{T}"/> class.
	/// </summary>
	/// <param name="N">The n.</param>
	/// <param name="dN">The d n.</param>
	/// <param name="defaultStepSize">Default size of the step.</param>
	/// <param name="hFunc">The h function.</param>
	/// <param name="checkFunc">The check function.</param>
	/// <param name="precision">The precision.</param>
	XCannon(const ProblemAbstract<T>& problem, const T defaultStepSize, const T precision, 
		std::function<bool(InitCondition<T>&)> checkFunc = [](InitCondition<T>& ic){ return true; }, 
		std::function<T(const int, const int)> hFunc = [](const int a, const int b){ return 1; }) : 
		XCannonAbstract(problem, defaultStepSize, precision, checkFunc, hFunc)
	{};

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="fileName">Name of the file.</param>
	virtual void SaveToFile(const char* fileName) override
	{
		SaveFunctionToFile(fileName, false /*invertMapping*/);
	};

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="saveFileStream">The save file stream.</param>
	virtual void SaveToFile(ofstream& saveFileStream) override
	{
		SaveFunctionToFile(saveFileStream, false /*invertMapping*/);
	};
};

#endif