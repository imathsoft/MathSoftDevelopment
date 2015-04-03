#ifndef GUARD_XCANNON_INVERSE
#define GUARD_XCANNON_INVERSE

#include "XCannonAbstract.h"

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
		T F = _N(prevKnot.Argument);
		T derivSquared = auxutils::sqr(prevKnot.Derivative);
		T A = - (_dN(prevKnot.Argument)*prevKnot.Argument + F - 2*auxutils::sqr(F*prevKnot.Argument)*derivSquared)*derivSquared;
		T B = - F*prevKnot.Argument*derivSquared;
		T C = prevKnot.Derivative;
		T D = prevKnot.Value;

		T hOpt = (abs(A) > 0) ? 1/abs(A) : abs(_h);

		T H = sgn(_h) * min(min(hOpt, abs(_h)), abs(argFinish - prevKnot.Argument));

		InitCondition<T> knot = XI_Func(A, B, C, D, H, _precision);
		knot.Argument = prevKnot.Argument + H;

		return knot;
	}

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
	XCannonInverse(std::function<T(const T&)> N, std::function<T(const T&)> dN, const T defaultStepSize, 
		std::function<T(const int, const int)> hFunc, std::function<bool(InitCondition<T>&)> checkFunc,
		double precision) : 
		XCannonAbstract(N, dN, defaultStepSize, 
		hFunc, checkFunc, precision)
	{
	}

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="fileName">Name of the file.</param>
	virtual void SaveToFile(const char* fileName) override
	{
		SaveFunctionToFile(fileName, true /*invertMapping*/);
	}

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="saveFileStream">The save file stream.</param>
	virtual void SaveToFile(ofstream& saveFileStream) override
	{
		SaveFunctionToFile(saveFileStream, true /*invertMapping*/);
	}
};

#endif