#ifndef GUARD_XCANNON_ABSTRACT
#define GUARD_XCANNON_ABSTRACT

#include "..\Problems\ProblemAbstract.h"
#include "..\FunctionApproximation\InitialCondition.h"
#include "..\Utils\AuxUtils.h"
#include <vector>
#include <list>

using namespace auxutils;

template <class T>
class BisectionComponent;

template <class T>
class XCannonAbstract
{
protected:
	typedef std::list<InitCondition<T>> ICLIST;
	ICLIST _listIC;
	typename ICLIST::iterator iLst;
	ProblemAbstract<T>* _problem;
	std::function<T(const T&, const T&, const T&)> _aCoeff;
	std::function<T(const T&, const T&, const T&)> _bCoeff;
	std::function<T(const int, const int)> _hFunc;
	std::function<bool(InitCondition<T>&)> _checkFunc;
	T _h; // default step size
	T _precision;

	virtual inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish) = 0;

	/// <summary>
	/// Saves the function to file.
	/// </summary>
	/// <param name="saveFile">The save file.</param>
	/// <param name="invertMapping">if set to <c>true</c> [invert mapping].</param>
	/// A method to save function values to the specified file
	void SaveFunctionToFile(std::ofstream& saveFile, bool invertMapping)
	{
		if (invertMapping)
		{
			for (iLst = _listIC.begin(); iLst != _listIC.end(); ++iLst)
			{
				auxutils::WriteToStream(saveFile, (*iLst).Value);
				auxutils::WriteToStream(saveFile, (*iLst).Argument);
				auxutils::WriteToStream(saveFile, 1/(*iLst).Derivative);
				saveFile <<  endl;
			}
		}
		else
		{
			for (iLst = _listIC.begin(); iLst != _listIC.end(); ++iLst)
			{
				auxutils::WriteToStream(saveFile, (*iLst).Argument);
				auxutils::WriteToStream(saveFile, (*iLst).Value);
				auxutils::WriteToStream(saveFile, (*iLst).Derivative);
				saveFile <<  endl;
			}
		}
	};

	/// <summary>
	/// Saves the function to file.
	/// </summary>
	/// <param name="fileName">Name of the file.</param>
	/// <param name="invertMapping">if set to <c>true</c> [invert mapping].</param>
	void SaveFunctionToFile(const char* fileName, bool invertMapping)
	{
		ofstream saveFile;
		saveFile.open(fileName);
     
		SaveFunctionToFile(saveFile, invertMapping);
		saveFile.close();
	};

	///Method to validate start and finish arguments with respect to step value
	///Returns false if it is impossible to get to the finish argument from the start argument 
	///using the given step size value
	bool ValidateArgumentSegmentAndStep(const T& argStart, const T& argFinish)
	{
		return ValidateArgumentSegmentAndStep(argStart, argFinish, _h);
	}

	bool ValidateArgumentSegmentAndStep(const T& argStart, const T& argFinish, const T& step)
	{
		return abs(step) > 0 &&  sgn(argFinish - argStart) == sgn(step);;
	}

	///Returns true if the two given arguments are equal with respect to the given precision
	bool ArgumentsAreEqualWithPrecision(const T& arg1, const T& arg2)
	{
		return abs(arg1 - arg2) <= 10*_precision;
	}

	bool _destinationPointIsHit;

	/// <summary>
	/// Shoots the specified argument start.
	/// </summary>
	/// <param name="argStart">The argument start.</param>
	/// <param name="argFinish">The argument finish.</param>
	/// <param name="funcStart">The function start.</param>
	/// <param name="funcTarget">The function target.</param>
	/// <param name="dFuncStart">The d function start.</param>
	/// <returns></returns>
	virtual InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& funcTarget, const T& dFuncStart)
	{
		_destinationPointIsHit = false;
		_listIC.clear();
		InitCondition<T> ic;
		ic.Value = funcStart;
		ic.Derivative = dFuncStart;
		ic.Argument = argStart;
		_listIC.insert(_listIC.end(), ic);

		if (!ValidateArgumentSegmentAndStep(argStart, argFinish))
			return ic;
					
		do
		{
			ic = GetNextKnot(ic, argFinish);
			_listIC.insert(_listIC.end(), ic);
		} while (!ArgumentsAreEqualWithPrecision(ic.Argument, argFinish) && _checkFunc(ic)); 

		_destinationPointIsHit = ArgumentsAreEqualWithPrecision(ic.Argument, argFinish);

		return ic;
	};

public:	
	/// <summary>
	/// Initializes a new instance of the <see cref="XCannonAbstract{T}"/> class.
	/// </summary>
	/// <param name="N">The n.</param>
	/// <param name="dN">The d n.</param>
	/// <param name="defaultStepSize">Default size of the step.</param>
	/// <param name="hFunc">The h function.</param>
	/// <param name="checkFunc">The check function.</param>
	/// <param name="precision">The precision.</param>
	XCannonAbstract(ProblemAbstract<T>& problem, const T defaultStepSize, const T precision, 
		std::function<bool(InitCondition<T>&)> checkFunc = [](InitCondition<T>& ic){ return true; }, 
		std::function<T(const int, const int)> hFunc = [](const int a, const int b){ return 1; })
	{
		_problem = &problem;
		_h = defaultStepSize;
		_hFunc = hFunc;
		_checkFunc = checkFunc;
		_precision = precision;
		_destinationPointIsHit = false;
	};

	/// <summary>
	/// Gets the precision.
	/// </summary>
	/// <returns></returns>
	T GetPrecision()
	{
		return _precision;
	};

	///Returns true if after shooting we hit the end of interval
	bool DestinationPointIsHit()
	{
		return _destinationPointIsHit;
	}

	//Returns vector of knots
	/// <summary>
	/// Gets the knot vector.
	/// </summary>
	/// <returns></returns>
	virtual std::vector<InitCondition<T>> GetKnotVector()
	{
		std::vector<InitCondition<T>> result(std::begin(_listIC), std::end(_listIC));
		return result;
	};

	///Returns knot vector with respect to straight solution
	virtual std::vector<InitCondition<T>> GetKnotVectorStreight()
	{
		return GetKnotVector();
	}

	/// <summary>
	/// Shoots the specified argument start.
	/// </summary>
	/// <param name="argStart">The argument start.</param>
	/// <param name="argFinish">The argument finish.</param>
	/// <param name="funcStart">The function start.</param>
	/// <param name="dFuncStart">The d function start.</param>
	/// <returns></returns>
	virtual InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& dFuncStart)
	{
		return Shoot(argStart, argFinish, funcStart, 0 /*does not matter*/, dFuncStart);
	};

	/// <summary>
	/// Finalizes an instance of the <see cref="XCannonAbstract{T}"/> class.
	/// </summary>
	virtual ~XCannonAbstract()
	{
		_listIC.clear();
	};

	virtual void SaveToFile(const char* fileName) = 0;

	virtual void SaveToFile(std::ofstream& saveFileStream) = 0;

	friend class BisectionComponent<T>;
};

#endif