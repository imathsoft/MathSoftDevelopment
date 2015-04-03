#ifndef GUARD_TROESCH_HYBRID_CANNON
#define GUARD_TROESCH_HYBRID_CANNON

#include "XCannonInverse.h"
#include "XCannon.h"

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
		InitCondition<T> result;
		return result;
	}

public:
	///Constructor
	TroeschHybridCannon(std::function<T(const T&)> N, std::function<T(const T&)> dN, const T defaultStepSize, 
		std::function<T(const int, const int)> hFunc, std::function<bool(InitCondition<T>&)> checkFunc, 
		double precision) : 
		XCannonAbstract(N, dN, defaultStepSize, 
		hFunc, checkFunc, precision)
	{
		_straightCannon = new XCannon<T>(N, dN, defaultStepSize, hFunc, TroeschHybridCannon<T>::StraightCannonCheckFunc, precision);
		_inverseCannon = new XCannonInverse<T>(N, dN, defaultStepSize, hFunc, checkFunc, precision);
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
		std::transform(++resultInv.begin(), resultInv.end(), ++resultInv.begin(), [](InitCondition<T> ic) ->
			InitCondition<T>{ 
				InitCondition<T> result;
				result.Value = ic.Argument;
				result.Argument = ic.Value;
				result.Derivative = 1/ic.Derivative;
				return result;
		});
		result.insert(result.end(), ++resultInv.begin(), resultInv.end());
		return result;
	};
};

#endif