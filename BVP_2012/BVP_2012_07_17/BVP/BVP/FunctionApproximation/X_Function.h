#ifndef GUARD_X_FUNCTION
#define GUARD_X_FUNCTION
#include <iostream>
#include <fstream>
#include <conio.h>
#include <Eigen/Dense>
#include <Eigen/MPRealSupport>
#include <array>
#include <list>
#include "GradientVector.h"
#include "InitialCondition.h"

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

template <class T>
inline T X_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const double precision){
	return X(A, B, C, D, h, EstimateIterations((abs(A)*abs(h)+abs(B))*abs(h*h), precision/(abs(C*h)+abs(D))));
}

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
}

template <class T>
inline InitCondition<T> X3(const T& A, const T& B, const T& C, const T& D, const T& h, int n){
	T* mas = new T[3*n+2];
	InitCondition<T> result;
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

	delete [] mas;
    return result;
}

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
}

// A method that computes u(h) where u(x) stisfies the following Cauchy problem
// u''(x) = (A*x+B)*u(x);
// u(0) = D;
// u'(0) = C;
template <class T>
inline InitCondition<T> X3_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const double precision){
	return X3(A, B, C, D, h, EstimateIterationsDouble((abs(A)*abs(h)+abs(B))*abs(h*h), precision/(abs(C*h)+abs(D))));
}

// A method that computes u(h) where u(x) stisfies the following Cauchy problem
// u''(x) = (A*x+B)*u'(x);
// u(0) = D;
// u'(0) = C;
template <class T>
inline InitCondition<T> XI_Func(const T& A, const T& B, const T& C, const T& D, const T& h, const double precision){
	return XI(A, B, C, D, h, EstimateIterationsSingle((abs(A)*abs(h)+abs(B))*abs(h),precision/(abs(C))));
}

#endif