#ifndef GUARD_TROESCH_HYBRID_CANNON
#define GUARD_TROESCH_HYBRID_CANNON

#include "XCannonInverse.h"
#include "XCannon.h"
#include "memory"
#include "../Utils/AuxUtils.h"

using namespace auxutils;

template <class T>
class HybridCannon : public XCannonAbstract<T>
{
private:
	typedef std::list<std::unique_ptr<XCannonAbstract<T>>> CANNON_LIST;

	CANNON_LIST _cannonList;

protected:
	inline InitCondition<T> GetNextKnot(const InitCondition<T> prevKnot, const T& argFinish)
	{
		InitCondition<T> result = {0,0,0,0};
		return result;
	};

public:
	///Constructor
	HybridCannon(ProblemAbstract<T>& problem, const T defaultStepSize, 
		T precision, 
		std::function<bool(InitCondition<T>&)> checkFunc = [](InitCondition<T>& ic){ return true; }, 
		std::function<T(const int, const int)> hFunc = [](const int a, const int b){ return 1; }) : 
		XCannonAbstract(problem, defaultStepSize, precision, checkFunc, 
		hFunc)
	{};

	~HybridCannon()
	{};

	virtual InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& funcTarget, const T& dFuncStart) override
	{
		_destinationPointIsHit = false;

		InitCondition<T> result = {funcStart, dFuncStart, 0, argStart};
		_cannonList.clear();

		if (!ValidateArgumentSegmentAndStep(argStart, argFinish))
			return result;

		auto subCheckFunction = [=](InitCondition<T>& ic){ 
			return abs(ic.Derivative) <= 1 && _checkFunc(ic); };

		while (_checkFunc(result) && 
			!_destinationPointIsHit)
		{
			if (abs(result.Derivative) < 1)
			{
				auto cannon = std::unique_ptr<XCannon<T>>(
					new XCannon<T>(*_problem, _h, _precision, subCheckFunction ));

				result = cannon->Shoot(result.Argument, argFinish, result.Value, result.Derivative);

				_destinationPointIsHit = cannon->DestinationPointIsHit();

				_cannonList.push_back(std::move(cannon));
			}
			else
			{
				T inverseStepSize;
				if (result.Derivative > 0)
					inverseStepSize = _h;
				else
					inverseStepSize = -_h;

				auto cannon = std::unique_ptr<XCannonInverse<T>>(
					new XCannonInverse<T>(*_problem, inverseStepSize, _precision, subCheckFunction));

				T inverseArgFinish;
				if (ValidateArgumentSegmentAndStep(result.Value, funcTarget, inverseStepSize))
				   inverseArgFinish = funcTarget;
				else 
					inverseArgFinish = sgn(inverseStepSize) * 1e6;

				InitCondition<T> tempResult = cannon->Shoot(result.Value, inverseArgFinish, result.Argument, 1/result.Derivative);
				result.Argument = tempResult.Value;
				result.Value    = tempResult.Argument;
				result.Derivative = 1/tempResult.Derivative;

				_destinationPointIsHit = cannon->DestinationPointIsHit();

				_cannonList.push_back(std::move(cannon));
			}

		}

		return result;
	};

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
	};

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
	};

	//Returns vector of knots
	virtual std::vector<InitCondition<T>> GetKnotVector() override
	{
		std::vector<InitCondition<T>> result;
		if (_cannonList.empty())
			return result;

		result = _cannonList.front()->GetKnotVectorStreight();
		for (auto iter = ++_cannonList.begin(); iter != _cannonList.end(); ++iter)
		{
			auto vector = (*iter)->GetKnotVectorStreight();
			result.insert(result.end(), ++vector.begin(), vector.end());
		}

		return result;
	};
};

#endif