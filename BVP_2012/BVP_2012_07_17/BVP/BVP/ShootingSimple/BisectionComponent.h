#ifndef GUARD_BISECTION_COMPONENT
#define GUARD_BISECTION_COMPONENT

#include "..\Cannon\XCannon.h"

template <class T>
class BisectionComponent
{
private:
	XCannonAbstract<T>* _cannon;
	double _precision;

	InitCondition<T> Shoot(const T& argStart, const T& argFinish, const T& funcStart, const T& uTarget, const T& dFuncStart)
	{
		return _cannon->Shoot(argStart, argFinish, funcStart, uTarget, dFuncStart);
	}

public:
	BisectionComponent(XCannonAbstract<T>& cannon)
	{
		_cannon = &cannon;
		_precision = _cannon->GetPrecision();
	}

	bool DerivativeBisection(const T& argStart, const T& argFinish, const T& uStart, const T& uTarget, const T& duLeft, const T& duRight)
	{
		std::function<int(const InitCondition<T>&)> evalFunc = [&](const InitCondition<T>& ic) { return sgn(ic.Value - uTarget); };

		return DerivativeBisectionGen(argStart, argFinish, uStart, uTarget, duLeft, duRight, evalFunc);
	}

	bool DerivativeBisectionGen(const T& argStart, const T& argFinish, const T& uStart, const T& uTarget, 
		const T& duLeft, const T& duRight, std::function<int(const InitCondition<T>&)>& evaluateResultFunc)
	{
		InitCondition<T> icLeft, icRight, icCurrent;
		T derLeft, derRight, derCurrent;

		icLeft  = Shoot(argStart, argFinish, uStart, uTarget, duLeft);
		icRight = Shoot(argStart, argFinish, uStart, uTarget, duRight);

		if (evaluateResultFunc(icLeft)*evaluateResultFunc(icRight) > 0)
			return false;

		derRight = duRight;
		derLeft = duLeft;

		while (abs(derRight - derLeft) > abs(derLeft)*_precision)
		{
			derCurrent = (derRight + derLeft)/2;
     		icCurrent  = Shoot(argStart, argFinish, uStart, uTarget, derCurrent);

			if (evaluateResultFunc(icCurrent)*evaluateResultFunc(icRight) > 0)
			{
				derRight = derCurrent;
				icRight  = icCurrent;
			}
			else
			{
				derLeft = derCurrent;
				icLeft  = icCurrent;
			}
		}

		return true;
	}

};

#endif