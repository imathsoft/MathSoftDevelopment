#ifndef GUARD_X_FUNCTION
#define GUARD_X_FUNCTION
#include <iostream>
#include <fstream>
#include <list>
#include "GradientVector.h"
#include "InitialCondition.h"
#include "../Utils/AuxUtils.h"

using namespace std;

//Estimates of iterations for algorithms with double integration
template <class T>
inline int EstimateIterationsDouble(const T& numerator, const T& precision)
{
	int i = 2;
	T t = numerator/2;

	while (t > precision)
	{
		t = t*numerator/((i+1)*(i+2));
		i=i+2;
	};
	return i/2;
};

//Estimates of iterations for algorithms with single integration
template <class T>
inline int EstimateIterationsSingle(const T& numerator, const T& precision)
{
	int i = 1;
	T t = numerator;

	while (t > precision)
	{
		t = t*numerator/(i+1);
		i=i++;
	};
	return i;
};

template <class T, class P>
inline T X_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const P& precision){
	return X(A, B, C, D, h, EstimateIterations((abs(A)*abs(h)+abs(B))*abs(h*h), precision/(abs(C*h)+abs(D))));
};

//This function calculates value of u(h), where
//diff(u(x), x, x) = (A*x+B)*u(x);   u'(0) = C; u(0)=D;
//Fixed point iteration method is used
//n -- number of iterations

template <class T>
inline T X(const T& A, const T& B, const T& C, const T& D, const T& h, int n){
	T* mas = new T[3*n+2];
	T result;
    int k;
    mas[0] = D;
    mas[1] = C;
    mas[2] = B*D/2;
    mas[3] = (B*C+A*D)/6;
    mas[4] = A*C/12;

    for (int i = 2; i <= n; i++) {
        k = 3 * i - 2;
        mas[k+3] = mas[k]*A/((k+3)*(k+2));
        for (int j = k; j >= 2 * i - 2; j--) { 
            mas[j+2] = (mas[j]*B+mas[j-1]*A)/((j+1)*(j+2));
        }
    }

    result =  mas[3 * n + 1] * h;
    for (int i = 3 * n; i > 0; i--) {
        result = result + mas[i];
		result = result*h;
    }
	result += D;
	delete [] mas;
    return result;
};

template <class T>
inline InitCondition<T> X3(const T& A, const T& B, const T& C, const T& D, const T& h, int n){
	std::vector<T> baseCoefficientBuffer(3*n+3);
	T* mas = &baseCoefficientBuffer[1];
	mas[-1] = T(0);
    mas[0] = D;
    mas[1] = C;

    int k;
    for (int i = 1; i <= n; i++) {
        k = 3 * i - 2;
        mas[k+3] = mas[k]*A/((k+3)*(k+2));
        for (int j = k; j >= 2 * i - 2; j--) { 
            mas[j+2] = (mas[j]*B+mas[j-1]*A)/((j+1)*(j+2));
        }
    }

	InitCondition<T> result;
	result.Value =  mas[3 * n + 1] * h;
    for (int i = 3 * n; i > 0; i--) {
		result.Value = (result.Value + mas[i]) * h;
    }
	result.Value += D;

	k = 3 * n + 1;
    result.Derivative = mas[k]*A*h/(k+2);
    for (int j = k; j >= 1; j--) { 
        result.Derivative = (result.Derivative + (mas[j]*B+mas[j-1]*A)/(j+1)) * h;
    }

	result.Derivative = (result.Derivative + mas[0] * B) * h + C;
	result.SecDerivative = (A * h + B) * result.Value;

    return result;
};

template <class T, class P>
inline InitCondition<T> X4_Func(const T& A, const T& B, const T& C, const T& D, const T& E, const T& F, const T& h, const P& precision){
	const int maxPower = 40;
	const int minPower = 10;
	std::vector<T> baseCoefficientBuffer(maxPower+2);
	T* mas = &baseCoefficientBuffer[1];
	mas[-1] = T(0);
    mas[0] = D;
    mas[1] = C;
	mas[2] = (B*D + F)/2;
	mas[3] = (A*D + B*C + E)/6;

	int currentPower = 3;
	const auto hAbs = abs(h);
	auto hPowerAbs = hAbs;

	int lengthOfTailBelowThreshold = 0;

	do
	{
		currentPower++;

		if (baseCoefficientBuffer.size() <= currentPower + 1)
		{
			baseCoefficientBuffer.resize(2*baseCoefficientBuffer.size());
			mas = &baseCoefficientBuffer[1];
		}

		mas[currentPower] = (mas[currentPower-2]*B+mas[currentPower-3]*A)/((currentPower-1)*(currentPower));
		hPowerAbs *= hAbs;

		if (abs(mas[currentPower])*hPowerAbs <= precision)
			lengthOfTailBelowThreshold++;
		else 
			lengthOfTailBelowThreshold = 0;

	} while(currentPower < minPower || lengthOfTailBelowThreshold < 4);

	InitCondition<T> result;
	result.Value =  mas[currentPower];
    for (int i = currentPower - 1; i >= 0; i--) {
		result.Value = result.Value * h + mas[i];
    }

    result.Derivative = mas[currentPower]*A*h/(currentPower+2);
    for (int j = currentPower; j >= 0; j--) { 
        result.Derivative = (result.Derivative + (mas[j]*B+mas[j-1]*A)/(j+1)) * h;
    }

	result.Derivative += (E * h / 2 + F) * h  + C;
	result.SecDerivative = (A * h + B) * result.Value + E*h + F;

    return result;
};

template <class T>
inline InitCondition<T> XI(const T& A, const T& B, const T& C, const T& D, const T& h, int n){
	T* mas = new T[2*n+1];
	InitCondition<T> result;
    int k;
    mas[0] = C;
    mas[1] = B*C;
    mas[2] = A*C/2;

    for (int i = 2; i <= n; i++) {
        k = 2 * i - 2;
        mas[k+2] = mas[k]*A/(k+2);
        for (int j = k; j >= i - 1; j--) { 
            mas[j+1] = (mas[j]*B+mas[j-1]*A)/(j+1);
        }
    }

	result.Derivative =  mas[2 * n] * h;
	result.Value = mas[2 * n] * h / (2 * n + 1);
    for (int i = 2 * n - 1; i > 0; i--) {
		result.Derivative = (result.Derivative + mas[i]) *h;
		result.Value = (result.Value + mas[i]/(i + 1)) *h;
    }
	result.Derivative = result.Derivative + C;
	result.Value = (result.Value  + C) * h + D;
	result.SecDerivative = (A * h + B) * result.Derivative;

	delete [] mas;
    return result;
};

///Returns value of the function x(u) at point x = h, where x(u) is the solution to the Cauchy problem
/// x''(u) = (A*u+B)*(x')^3,  x'(0) = C,  x(0) = D;
///The value is calculated with the given "precision"
template <class T, class P>
inline InitCondition<T> XI_Bernoulli(const T& A, const T& B, const T& C, const T& D, const T& h, const P& precision)
{
	const int maxPower = 40;
	const int minPower = 10;
	std::vector<T> sol_series(maxPower + 1);
	std::vector<T> sol_series_sqr(maxPower + 1);

	sol_series[0] = C;
	T prev_cube_coeff = T(0);
	T next_cube_coeff = T(0);

	int lengthOfTailBelowThreshold = 0;
	int current_index = 0;

	do
	{
		sol_series_sqr[current_index] = T(0);
		for (int i = 0; i <= current_index; i++)
		{
			sol_series_sqr[current_index] += sol_series[i]*sol_series[current_index - i];
		}

		prev_cube_coeff = next_cube_coeff;

		next_cube_coeff = T(0);

		for (int i = 0; i<= current_index; i++)
		{
			next_cube_coeff += sol_series[i]*sol_series_sqr[current_index - i];
		}

		current_index++;

		if (current_index >= 1000)
			throw "Series does not converge";

		if (sol_series.size() <= current_index)
		{
			sol_series.resize(2*sol_series.size());
			sol_series_sqr.resize(2*sol_series.size());
		}

		sol_series[current_index] = h*(A*h*prev_cube_coeff + B*next_cube_coeff)/T(current_index);

		if (abs(sol_series[current_index]) < precision)
		{
			lengthOfTailBelowThreshold++;
		}
		else
		{
			lengthOfTailBelowThreshold = 0;
		}

	} while (lengthOfTailBelowThreshold < 4 || current_index < minPower);


	InitCondition<T> result;
	result.Value = T(0);
	result.Derivative = T(0);

	for (int index = current_index; index >= 0; index--)
	{
		result.Derivative += sol_series[index];
		result.Value += h*sol_series[index]/T(index+1);
	}

	result.Value += D;

	result.SecDerivative = (A*h+B)*result.Derivative*result.Derivative*result.Derivative;

	return result;
}

// A method that computes u(h) where u(x) stisfies the following Cauchy problem
// u''(x) = (A*x+B)*u(x);
// u(0) = D;
// u'(0) = C;
template <class T, class P>
inline InitCondition<T> X3_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const P& precision){
	return X3(A, B, C, D, h, EstimateIterationsDouble((abs(A)*abs(h)+abs(B))*abs(h*h), precision/(abs(C*h)+abs(D))));
};

// A method that computes u(h) where u(x) stisfies the following Cauchy problem
// u''(x) = (A*x+B)*u'(x);
// u(0) = D;
// u'(0) = C;
template <class T, class P>
inline InitCondition<T> XI_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const P& precision){
	return XI(A, B, C, D, h, EstimateIterationsSingle((abs(A)*abs(h)+abs(B))*abs(h),precision/(abs(C))));
};

template <class T>
class X_Func_Gradient 
{
public :
	T X;
	T dA;
	T dB;
	T dC;
	T dD;
	T dE;
	T dF;
	T dh;
	T dhdh;
	T dhdA;
	T dhdB;
	T dhdC;
	T dhdD;
	T dhdE;
	T dhdF;

	///Overloaded assignment operator
	inline X_Func_Gradient<T>& operator=(const InitCondition<typename GVTypes<T>::OptimalGradientVector>& fullGradient)
	{
		X = fullGradient.Value[0];
		dA = fullGradient.Value[1];
		dB = fullGradient.Value[2];
		dC = fullGradient.Value[3];
		dD = fullGradient.Value[4];
		dE = T(0);
		dF = T(0);

		dh = fullGradient.Derivative[0];
		dhdA = fullGradient.Derivative[1];
		dhdB = fullGradient.Derivative[2];
		dhdC = fullGradient.Derivative[3];
		dhdD = fullGradient.Derivative[4];
		dhdE = T(0);
		dhdF = T(0);

		dhdh = fullGradient.SecDerivative[0];

		return *this;
	}

	///Overloaded assignment operator
	inline X_Func_Gradient<T>& operator=(const InitCondition<typename GVTypes<T>::OptimalGradientVectorExtended>& fullGradient)
	{
		X = fullGradient.Value[0];
		dA = fullGradient.Value[1];
		dB = fullGradient.Value[2];
		dC = fullGradient.Value[3];
		dD = fullGradient.Value[4];
		dE = fullGradient.Value[5];
		dF = fullGradient.Value[6];

		dh = fullGradient.Derivative[0];
		dhdA = fullGradient.Derivative[1];
		dhdB = fullGradient.Derivative[2];
		dhdC = fullGradient.Derivative[3];
		dhdD = fullGradient.Derivative[4];
		dhdE = fullGradient.Derivative[5];
		dhdF = fullGradient.Derivative[6];

		dhdh = fullGradient.SecDerivative[0];

		return *this;
	}

    ///Returns gradient of X3_Func(...)
	static inline typename X_Func_Gradient<T> X3_Func_Gradient(const T& A, const T& B, const T& C, const T& D, const T& E, const T& F, const T& h, const T& precision)
	{
		GVTypes<T>::OptimalGradientVectorExtended EH; EH << h << 0 << 0 << 0 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVectorExtended EA; EA << A << 1 << 0 << 0 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVectorExtended EB; EB << B << 0 << 1 << 0 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVectorExtended EC; EC << C << 0 << 0 << 1 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVectorExtended ED; ED << D << 0 << 0 << 0 << 1 << 0 << 0;
		GVTypes<T>::OptimalGradientVectorExtended EE; EE << E << 0 << 0 << 0 << 0 << 1 << 0;
		GVTypes<T>::OptimalGradientVectorExtended EF; EF << F << 0 << 0 << 0 << 0 << 0 << 1;

		X_Func_Gradient<T> result;

		return result = X4_Func<typename GVTypes<T>::OptimalGradientVectorExtended, T>(EA, EB, EC, ED, EE, EF, EH, precision);
	}

	///Returns gradient of XI_Func_Gradient
	static inline typename X_Func_Gradient<T> XI_Func_Gradient(const T& A, const T& B, const T& C, const T& D, const T& h, const T& precision)
	{
		GVTypes<T>::OptimalGradientVector EH; EH << h << 0 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVector EA; EA << A << 1 << 0 << 0 << 0;
		GVTypes<T>::OptimalGradientVector EB; EB << B << 0 << 1 << 0 << 0;
		GVTypes<T>::OptimalGradientVector EC; EC << C << 0 << 0 << 1 << 0;
		GVTypes<T>::OptimalGradientVector ED; ED << D << 0 << 0 << 0 << 1;

		X_Func_Gradient<T> result;

		return result = XI_Func<typename GVTypes<T>::OptimalGradientVector, T>(EA, EB, EC, ED, EH, precision);
	}
};

#endif