#ifndef GUARD_BVP_PC_T
#define GUARD_BVP_PC_T
#include <iostream>
#include <conio.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/MPRealSupport>
#include <string>
#include <fstream>
//#include "../../mpreal_parser/fparser.hh"
//#include "../../mpreal_parser/fparser_mpfr.hh"

const int Cheb_Intervals_N=12;

typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic,1>        VectorXmp; // This is a dynamic vector type
typedef Eigen::Matrix<mpfr::mpreal, 2,1>        Vector2mp; // This is a 2-dimensional vector type
typedef Eigen::Matrix<mpfr::mpreal, 3,1>        Vector3mp; // This is a 3-dimensional vector type
typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic,Eigen::Dynamic>  MatrixXmp; // This is a dynamic matrix type
typedef Eigen::Matrix<mpfr::mpreal, 2,2>  Matrix2x2mp; // This is a 2x2 matrix type
typedef Eigen::SparseMatrix<mpfr::mpreal>  SMatrixXmp;// This is a sparse matrix type for the multiple shooting method

struct IntervalData {
   Vector2mp U;
   Vector2mp C;
   mpfr::mpreal a;
   mpfr::mpreal b;
   mpfr::mpreal x;
};
struct mp_pair {
	mpfr::mpreal value;
	mpfr::mpreal gamma;
};

class T_BVP_PC_T{
private:
  mpfr::mpreal lambda;//The value of Lambda in the Troesch test
  mpfr::mpreal u0;// The left boundary value
  mpfr::mpreal u1;// The right boundary value
  std::string Nonlin;// The expression for nonlinearity;
//  FunctionParser_mpfr PNonlin;//Parser for the nonlinearity;
  mpfr::mpreal h;// Discretization step;
  mpfr::mpreal* h_array; // The array for storing step sizes of a non-uniform mesh 
  mpfr::mpreal x0;// Left point  
  mpfr::mpreal x1;// Right point
  int n;// The number of discretization intervals
  mpfr::mpreal** bpU; // Coefficients for the solution of the basic problem 
  IntervalData* bpUD; // The Data about basic problem with adaptive algorithm
  VectorXmp C;// The vector of coefficients for the Multiple Shooting Ptocedure 
  mpfr::mpreal (*F)(const mpfr::mpreal&);// Pointer to function that represents the nonlinearity
  mpfr::mpreal (*dF)(const mpfr::mpreal&);// Pointer to function that represents the derivative of nonlinearity
  mpfr::mpreal (*hF)(const mpfr::mpreal&);// Pointer to a strictly increasing (!!!) function [0,1]<->[0,1] that defines a non-uniform mesh
//==========================================================================================================
  //Tis is A and B functions which we will use in the adaptive algorithm for solving the basic problem
//==========================================================================================================
  mp_pair _1div3, _2div3, _4div3, _5div3;
  MatrixXmp MCH_1div3[Cheb_Intervals_N];
  MatrixXmp MCH_2div3[Cheb_Intervals_N];
  MatrixXmp MCH_4div3[Cheb_Intervals_N];
  MatrixXmp MCH_5div3[Cheb_Intervals_N];
  static int Cheb_Cof_N;
  static mpfr::mpreal Tan_threshold; // The value that defines a threshold between pice-wise constant and pice-wise linear approximations of the ninlinearity
  int LoadChebPadeCoeffs();
  mpfr::mpreal F01(const mp_pair& a1, const mpfr::mpreal& z); //This is the hypergeometric function 0F1
  mpfr::mpreal F01_1div3(const mpfr::mpreal& z);// This is the hypergeometric function 0F1(;1/3;z);
  mpfr::mpreal F01_2div3(const mpfr::mpreal& z);// This is the hypergeometric function 0F1(;2/3;z);
  mpfr::mpreal F01_4div3(const mpfr::mpreal& z);// This is the hypergeometric function 0F1(;4/3;z);
  mpfr::mpreal F01_5div3(const mpfr::mpreal& z);// This is the hypergeometric function 0F1(;5/3;z);
  mpfr::mpreal A(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z); // This is the A-function (see description below) 
  mpfr::mpreal B(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z); // This is the B-function (see description below) 
  mpfr::mpreal dA(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z); // This is the derivative of A-function (see description below) 
  mpfr::mpreal dB(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z); // This is the derivative of B-function (see description below) 
  int GetMatchMatrix(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z, Matrix2x2mp& M); // The procedure that returns the "match matrix" 
//==========================================================================================================
public:
  T_BVP_PC_T(const std::string& str_l, const std::string& str_u0, const std::string& str_u1, const int& _n, 
	  mpfr::mpreal (*MF)(const mpfr::mpreal&), mpfr::mpreal (*dMF)(const mpfr::mpreal&), mpfr::mpreal (*MhF)(const mpfr::mpreal&));// Constructor
  ~T_BVP_PC_T();//Destructor
  mpfr::mpreal ShootFunc(const mpfr::mpreal& du0); //The shooting function for solving the basic problem
  mpfr::mpreal AB_ShootFunc(const mpfr::mpreal& du0); //The shooting function for solving the basic problem with the AB-algorithm;
  int BisectionProcedure(const mpfr::mpreal& du0_left, const mpfr::mpreal& du0_right); //This function implements a bisection procedure for solving the basic problem
                                                   //the only argument it takes is a tangent of the unknown solution u^{(0)}(x) 
												   //at the begining of the interval
  int AdjustmentOfFireProcedure(const mpfr::mpreal& du0_left, const mpfr::mpreal& du0_right); 
												   //This function implements a bisection procedure for solving the basic problem
                                                   //but in opposite to the BisectionProcedure its purpose is to obtain a 0-approximation  C0
												   //for the MultipleShootingProcedure presented below
  int MultipleShootingProcedure(VectorXmp& C0); // This function implements a multiple shooting procedure for solving the basic problem
                                                      // Vector C0 is an initial approximation for the Newton method used to solve a nonlinear 

  int BasicProblemSolve(); // The procedure that solves the basic problem using Single and Multiple Shooting Procedures
  int BasicProblemSolveSimple(); // The procedure that solves the basic problem using only Single Shooting Procedure
  class BadInitData { }; //Exception class 
  class MSPBadInitData { }; //Exception class Bad initial data in the Multiple Shooting procedure
  int SaveToFile(const std::string& filename); 
  int SaveToFile(const std::string& filename, const int& substeps); 
};

#endif