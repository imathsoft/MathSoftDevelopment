#ifndef GUARD_XCANNON_INVERSE
#define GUARD_XCANNON_INVERSE

#include "XCannonAbstract.h"
#include "../FunctionApproximation/X_Function.h"

template <class T>
class XCannonInverse : public XCannonAbstract<T>
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
		T A = _aCoeff(prevKnot.Derivative, prevKnot.Value, prevKnot.Argument);
		T B = _bCoeff(prevKnot.Derivative, prevKnot.Value, prevKnot.Argument);
    	T C = prevKnot.Derivative;
		T D = prevKnot.Value;

		T hOpt = (abs(A) > 0) ? 1/abs(A) : abs(_h);

		T H = sgn(_h) * min(min(hOpt, abs(_h)), abs(argFinish - prevKnot.Argument));

		InitCondition<T> knot = XI_Func(A, B, C, D, H, _precision);
		knot.Argument = prevKnot.Argument + H;

		return knot;
	};

public:	
	/// <summary>
	/// Initializes a new instance of the <see cref="XCannonInverse{T}"/> class.
	/// </summary>
	/// <param name="N">The n.</param>
	/// <param name="dN">The d n.</param>
	/// <param name="defaultStepSize">Default size of the step.</param>
	/// <param name="hFunc">The h function.</param>
	/// <param name="checkFunc">The check function.</param>
	/// <param name="precision">The precision.</param>
	XCannonInverse(ProblemAbstract<T>& problem, const T defaultStepSize, double precision, 
		std::function<bool(InitCondition<T>&)> checkFunc = [](InitCondition<T>& ic){ return true; }, 
		std::function<T(const int, const int)> hFunc = [](const int a, const int b){ return 1; }) : 
		XCannonAbstract(problem, defaultStepSize, precision, checkFunc, hFunc)
	{
		_aCoeff = _problem->GetACoeffInverse();
		_bCoeff = _problem->GetBCoeffInverse();
	};

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="fileName">Name of the file.</param>
	virtual void SaveToFile(const char* fileName) override
	{
		SaveFunctionToFile(fileName, true /*invertMapping*/);
	};

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="saveFileStream">The save file stream.</param>
	virtual void SaveToFile(ofstream& saveFileStream) override
	{
		SaveFunctionToFile(saveFileStream, true /*invertMapping*/);
	};
};

#endif