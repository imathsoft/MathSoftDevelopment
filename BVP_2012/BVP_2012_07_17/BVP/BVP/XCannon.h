#ifndef GUARD_X_CANNON
#define GUARD_X_CANNON

#include "X_Function.h"
#include "AuxUtils.h"
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
	T (*_N)(const T&); 
	T (*_dN)(const T&);
	T (*_hFunc)(const int, const int);
	bool (*_checkFunc)(InitCondition<T>&);
	T _h; // default step size
	double _precision;

	virtual inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish) = 0;

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

	void SaveFunctionToFile(const char* fileName, bool invertMapping)
	{
		ofstream saveFile;
		saveFile.open(fileName);
     
		SaveFunctionToFile(saveFile, invertMapping);
		saveFile.close();
	}

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
	XCannonAbstract(T (*N)(const T&), T (*dN)(const T&), const T defaultStepSize, 
		T (*hFunc)(const int, const int), bool (*checkFunc)(InitCondition<T>&), double precision)
	{
		_h = defaultStepSize;
		_N = N;
		_dN = dN;
		_hFunc = hFunc;
		_checkFunc = checkFunc;
		_precision = precision;
	}

	double GetPrecision()
	{
		return _precision;
	}

	//Returns vector of knots
	virtual std::vector<InitCondition<T>> GetKnotVector()
	{
		std::vector<InitCondition<T>> result(std::begin(_listIC), std::end(_listIC));
		return result;
	};

    virtual InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& dFuncStart)
	{
		return Shoot(argStart, argFinish, funcStart, 0 /*does not matter*/, dFuncStart);
	}

	virtual ~XCannonAbstract()
	{
		_listIC.clear();
	}

	virtual void SaveToFile(const char* fileName) = 0;

	virtual void SaveToFile(ofstream& saveFile) = 0;

	friend class BisectionComponent<T>;
};

template <class T>
class XCannon : public XCannonAbstract<T>
{
protected:

	virtual inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish) override
	{
		T A = _dN(prevKnot.Value)*prevKnot.Derivative;
		T B = _N(prevKnot.Value);
		T C = prevKnot.Derivative;
		T D = prevKnot.Value;

		T hOpt = (abs(A) > 1) ? 1/auxutils::RoughSqrt(abs(A)) : abs(_h);

		T H = sgn(_h) * min(min(hOpt, abs(_h)), abs(argFinish - prevKnot.Argument));

		InitCondition<T> knot = X3_Func(A, B, C, D, H, _precision);
		knot.Argument = prevKnot.Argument + H;

		return knot;
	}

public:	
	XCannon(T (*N)(const T&), T (*dN)(const T&), const T defaultStepSize, 
		T (*hFunc)(const int, const int), bool (*checkFunc)(InitCondition<T>&), double precision) : 
		XCannonAbstract(N, dN, defaultStepSize, 
		hFunc, checkFunc, precision)
	{
	}

	virtual void SaveToFile(const char* fileName) override
	{
		SaveFunctionToFile(fileName, false /*invertMapping*/);
	}

	virtual void SaveToFile(ofstream& saveFile) override
	{
		SaveFunctionToFile(saveFile, false /*invertMapping*/);
	}
};

template <class T>
class XCannonInverse : public XCannonAbstract<T>
{
protected:

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
	XCannonInverse(T (*N)(const T&), T (*dN)(const T&), const T defaultStepSize, 
		T (*hFunc)(const int, const int), bool (*checkFunc)(InitCondition<T>&), double precision) : 
		XCannonAbstract(N, dN, defaultStepSize, 
		hFunc, checkFunc, precision)
	{
	}

	virtual void SaveToFile(const char* fileName) override
	{
		SaveFunctionToFile(fileName, true /*invertMapping*/);
	}

	virtual void SaveToFile(ofstream& saveFile) override
	{
		SaveFunctionToFile(saveFile, true /*invertMapping*/);
	}
};

#endif