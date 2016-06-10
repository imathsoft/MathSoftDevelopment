# MathSoftDevelopment
The project is focused on developing algorithms for solving stiff boundary value problems (BVP) and initial value problems (IVP) for ordinary differential equations (ODE).

Currently it containes an implementation of a new algorithm for solving stiff BVPs for the second order ODE's, such like the well known Troesch's problem (which is an acknowledged banchmark for testing BVP solvers).
The capabilities of the algorithm (lets call it SI-algorithm, which means "Straight-Inverse", although the name might change in future)
allows it to solve the Troesch's problem for 'lambda = 100, 200, ... with precision 1e-8 within a few seconds (and this is in double precision !!!). 
For more information see this paper https://arxiv.org/abs/1601.04272 (draft)
#MathSoftDevelopment
