#ifndef GUARD_XCANNON_ABSTRACT
#define GUARD_XCANNON_ABSTRACT

#include "..\FunctionApproximation\X_Function.h"
#include <vector>

template <class T>
class BisectionComponent;

template <class T>
class XCannonAbstract
{
protected:
	typedef list<InitCondition<T>> ICLIST;
	ICLIST _listIC;
	typename ICLIST::iterator iLst;
	std::function<T(const T&)> _N;
	std::function<T(const T&)> _dN;
	std::function<T(const int, const int)> _hFunc;
	std::function<bool(InitCondition<T>&)> _checkFunc;
	T _h; // default step size
	double _precision;

	virtual inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish) = 0;

	/// <summary>
	/// Saves the function to file.
	/// </summary>
	/// <param name="saveFile">The save file.</param>
	/// <param name="invertMapping">if set to <c>true</c> [invert mapping].</param>
	/// A method to save function values to the specified file
	void SaveFunctionToFile(ofstream& saveFile, bool invertMapping)
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
	}

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
	}

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
		_listIC.clear();
		InitCondition<T> ic;
		ic.Value = funcStart;
		ic.Derivative = dFuncStart;
		ic.Argument = argStart;
		_listIC.insert(_listIC.end(), ic);
					
		do
		{
			ic = GetNextKnot(ic, argFinish);
			_listIC.insert(_listIC.end(), ic);
		} while (abs(ic.Argument- argFinish) > 10*_precision && _checkFunc(ic)); 

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
	XCannonAbstract(std::function<T(const T&)> N, std::function<T(const T&)> dN, const T defaultStepSize, 
		std::function<T(const int, const int)> hFunc, std::function<bool(InitCondition<T>&)> checkFunc, double precision)
	{
		_h = defaultStepSize;
		_N = N;
		_dN = dN;
		_hFunc = hFunc;
		_checkFunc = checkFunc;
		_precision = precision;
	}

	/// <summary>
	/// Gets the precision.
	/// </summary>
	/// <returns></returns>
	double GetPrecision()
	{
		return _precision;
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
	}

	/// <summary>
	/// Finalizes an instance of the <see cref="XCannonAbstract{T}"/> class.
	/// </summary>
	virtual ~XCannonAbstract()
	{
		_listIC.clear();
	}

	virtual void SaveToFile(const char* fileName) = 0;

	virtual void SaveToFile(ofstream& saveFileStream) = 0;

	friend class BisectionComponent<T>;
};

#endif