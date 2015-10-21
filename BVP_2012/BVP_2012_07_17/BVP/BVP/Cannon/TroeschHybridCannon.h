#ifndef GUARD_TROESCH_HYBRID_CANNON
#define GUARD_TROESCH_HYBRID_CANNON

#include "XCannonInverse.h"
#include "XCannon.h"
#include "../Utils/AuxUtils.h"

template <class T>
class TroeschHybridCannon : public XCannonAbstract<T>
{
private:
	XCannon<T>*        _straightCannon;
	XCannonInverse<T>* _inverseCannon;

	static bool StraightCannonCheckFunc(InitCondition<T>& ic)
	{
		return (abs(ic.Derivative) <= 1);
	}

protected:
	inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish)
	{
		InitCondition<T> result = {0,0,0,0};
		return result;
	}

public:
	///Constructor
	TroeschHybridCannon(ProblemAbstract<T>& problem, const T defaultStepSize, 
		double precision, 
		std::function<bool(InitCondition<T>&)> checkFunc = [](InitCondition<T>& ic){ return true; }, 
		std::function<T(const int, const int)> hFunc = [](const int a, const int b){ return 1; }) : 
		XCannonAbstract(problem, defaultStepSize, precision, checkFunc, 
		hFunc)
	{
		_straightCannon = new XCannon<T>(problem, defaultStepSize, precision, TroeschHybridCannon<T>::StraightCannonCheckFunc);
		_inverseCannon = new XCannonInverse<T>(problem, defaultStepSize, precision);
	}

	~TroeschHybridCannon()
	{
		delete _straightCannon;
		delete _inverseCannon;
	}

	virtual InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& funcTarget, const T& dFuncStart) override
	{
		InitCondition<T> result = _straightCannon->Shoot(argStart, argFinish, funcStart,  dFuncStart);
		if (abs(argFinish - result.Argument) > 10*_precision && abs(result.Derivative) > 10*_precision )
		{
			InitCondition<T> tempResult = _inverseCannon->Shoot(result.Value, funcTarget, result.Argument, 1/result.Derivative);
			result.Argument = tempResult.Value;
			result.Value    = tempResult.Argument;
			result.Derivative = 1/tempResult.Derivative;
		}
		return result;
	}

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="fileName">Name of the file.</param>
	virtual void SaveToFile(const char* fileName) override
	{
		ofstream saveFile;
		saveFile.open(fileName);
     
		SaveToFile(saveFile);
		saveFile.close();
	}

	/// <summary>
	/// Saves to file.
	/// </summary>
	/// <param name="saveFile">Output file stream.</param>
	virtual void SaveToFile(ofstream& saveFile) override
	{
		auto knots = GetKnotVector();

		for each (auto knot in knots)
			{
				auxutils::WriteToStream(saveFile, knot.Argument);
				auxutils::WriteToStream(saveFile, knot.Value);
				auxutils::WriteToStream(saveFile, knot.Derivative);
				saveFile <<  endl;
			}
	}

	//Returns vector of knots
	virtual std::vector<InitCondition<T>> GetKnotVector() override
	{
		std::vector<InitCondition<T>> result = _straightCannon->GetKnotVector();
		std::vector<InitCondition<T>> resultInv = _inverseCannon->GetKnotVector();
		if (resultInv.size() > 0)
		{
			std::transform(++resultInv.begin(), resultInv.end(), ++resultInv.begin(), [](InitCondition<T> ic) ->
				InitCondition<T>{ 
					InitCondition<T> result;
					result.Value = ic.Argument;
					result.Argument = ic.Value;
					result.Derivative = 1/ic.Derivative;
					result.SecDerivative = - ic.SecDerivative/(auxutils::sqr(ic.Derivative) * ic.Derivative);
					return result;
			});
			result.insert(result.end(), ++resultInv.begin(), resultInv.end());
		}
		return result;
	};
};

#endif