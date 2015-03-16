#include "BVP_PC_T.h"


using namespace std;
using namespace mpfr;
using namespace Eigen;


int T_BVP_PC_T::Cheb_Cof_N=21; // We will use Chebyshev approximations for 0F1 vith 89 coefficients
mpfr::mpreal T_BVP_PC_T::Tan_threshold=1000; // The threshold value that defines the type of approximation used

T_BVP_PC_T::T_BVP_PC_T(const std::string& str_l, const std::string& str_u0, const std::string& str_u1, const int& _n, mpfr::mpreal (*MF)(const mpfr::mpreal&), 
					   mpfr::mpreal (*dMF)(const mpfr::mpreal&), mpfr::mpreal (*MhF)(const mpfr::mpreal&)){
  //Initialization of the class
  lambda=str_l;
  u0=str_u0;
  u1=str_u1;
  x0="0";
  x1="1";
  if (_n<=0) {n=1; throw BadInitData();}// Invalid number of discretization
  else { n=_n;};
  C.resize(2*n,1);// Vector of coefficients for Multiple Shooting Procedure
  h=(x1-x0)/n;
  F=MF;
  dF=dMF;
  bpU= new mpreal*[n+1];
  for (int i=0; i<=n; i++){
      bpU[i]= new mpreal[3];
  }
  bpUD= new IntervalData[n+1];
  //===================== 
  // The definition of a (non)-uniform mesh;
  hF=MhF;
  h_array=new mpreal[n]; // Because we have exactly n subintervals 
  mpreal x, y0, y1;
  x=x0;
  bpUD[0].x=x0;
  y0=hF(x);
  for (int i=0; i<n; i++){
      y1=hF(x+h);
      h_array[i]=y1-y0;
      bpUD[i+1].x=bpUD[i].x+h_array[i]; 
	  y0=y1;
	  x=x+h;
  }
  mpfr_free_cache();
  //===================== 
// Definition of some constants which we will use for computation of A abd B-functions 
  _1div3.value=static_cast<mpreal>(1)/3;// 1/3
  _1div3.gamma=gamma(_1div3.value);

  _2div3.value=static_cast<mpreal>(2)/3;// 2/3 
  _2div3.gamma=gamma(_2div3.value);

  _4div3.value=static_cast<mpreal>(4)/3;// 4/3
  _4div3.gamma=gamma(_4div3.value);

  _5div3.value=static_cast<mpreal>(5)/3;// 5/3
  _5div3.gamma=gamma(_5div3.value);

  //===================== 
  //===================== 
  LoadChebPadeCoeffs(); // Loading the Cebyshev-Pade coefficients (see ChebPadeCoefficients.cpp)
  //===================== 
  //===================== 
};
  mpreal T_BVP_PC_T::ShootFunc(const mpfr::mpreal& du0){
    mpreal hNi, ShSqN_i, ChSqN_i;
	int i;
    bpU[0][0]=u0;
	bpU[0][1]=du0;
	for (i=0; i<n; i++) {
		bpU[i][2]=sqrt(F(bpU[i][0]));
		hNi=h_array[i]*bpU[i][2];
	    bpU[i][1]=bpU[i][1]/bpU[i][2];
		ShSqN_i=sinh(hNi);
		ChSqN_i=cosh(hNi);
        // u^{(0)}(t)=C1*sinh(t*Ni)+C2*cosh(t*Ni);
        // du^{(0)}(t)=Ni*(C2*sinh(t*Ni)+C1*cosh(t*Ni)); t\in [0, h]
		bpU[i+1][0]=bpU[i][1]*ShSqN_i+bpU[i][0]*ChSqN_i;
		bpU[i+1][1]=bpU[i][2]*(bpU[i][0]*ShSqN_i+bpU[i][1]*ChSqN_i);
	}
    return bpU[n][0];
  };
  mpfr::mpreal T_BVP_PC_T::AB_ShootFunc(const mpfr::mpreal& du0){ //The shooting function for solving the basic problem with the AB-algorithm;
/**	mpreal a,b,z;
	for (int i=0; i<10; i++){
       cin>>a;
	   cin>>b;
	   cin>>z;
		cout<<dB(a, b, z)<<endl;
	}

**/	//cout<<F01(_4div3, -2550)<<endl;
	//  cout<<CH_1div3<<endl;
/**	  mpreal a;
	  for (int i=0; i<16; i++){
		  cout<<"a="<<endl;
		  cin>>a;
	      cout<<"F01_1div3(a)="<<F01_1div3(a)<<endl;
	  }
**/

  Vector2mp U, C;
  Matrix2x2mp M;
  mpreal a,b,z;
  
  bpUD[0].U(0,0)=u0;
  bpUD[0].U(1,0)=du0;
  for (int i=0; i<n; i++){
     bpUD[i].a=dF(bpUD[i].U(0,0))*bpUD[i].U(1,0);
     bpUD[i].b=F(bpUD[i].U(0,0));
     GetMatchMatrix(bpUD[i].a, bpUD[i].b, bpUD[i].x, M);
     bpUD[i].C=M.colPivHouseholderQr().solve(bpUD[i].U);
     cout<<"i="<<i<<endl;
     cout<<"M=\n"<<M<<endl;
     cout<<"bpUD[i].U=\n"<<bpUD[i].U<<endl;
     cout<<"bpUD[i].C=\n"<<bpUD[i].C<<endl;
     cout<<"bpUD[i].a="<<bpUD[i].a<<endl;
     cout<<"bpUD[i].b="<<bpUD[i].b<<endl;
     cout<<"bpUD[i].z="<<bpUD[i].x<<endl;
     cout<<endl;
	 GetMatchMatrix(bpUD[i].a,bpUD[i].b,bpUD[i+1].x,M);
	 bpUD[i+1].U=M*bpUD[i].C;
     cout<<endl;
  }

return 0;
  }

//================================================================================
//This function implements a bisection procedure for solving the basic problem
//the only argument it takes is a tangent of the unknown solution u^{(0)}(x) 
//at the begining of the interval 
//================================================================================
  int T_BVP_PC_T::BisectionProcedure(const mpfr::mpreal& du0_left, const mpfr::mpreal& du0_right){    
	mpreal old_target, target, left_target, right_target, eps, du, left_du, right_du;
	int i=0;
	  left_du=du0_left;
	  right_du=du0_right;
      right_target=ShootFunc(right_du);
      left_target=ShootFunc(left_du);
	  if (right_target<left_target) { 
	     du=left_du;
		 left_du=right_du;
         right_du=du;
         target=right_target;
		 right_target=left_target;
		 left_target=target;
	  }
	  if ((right_target>=u1) && (left_target<=u1)) {// We have to ensure that the interval [du0_left, du0_right] contains the tangent of unknown solution
		  du=(right_du+left_du)/2;
		  target=ShootFunc(du);
		  old_target=u1;
		  while (old_target!=target){
			  i++;
			  if (target>=u1) {
				  right_du=du; right_target=target;
			  }  else {
				  left_du=du; left_target=target;	
			  };
		      du=(right_du+left_du)/2;
              old_target=target;
		      target=ShootFunc(du);
			  cout<<"Target="<<target<<endl;
		  }
	  } else {
		  cout<<"Invalid interval"<<endl;
		  cout<<"T0="<<left_target<<"; T1="<<right_target<<endl;
		  return 0;
	  }
      return i;
	}
  int T_BVP_PC_T::AdjustmentOfFireProcedure(const mpfr::mpreal& du0_left, const mpfr::mpreal& du0_right){ //This function implements a bisection procedure for solving the basic problem
                                                   //but in opposite to the BisectionProcedure its purpose is to obtain a 0-approximation  C0
												   //for the MultipleShootingProcedure presented below
	mpreal target, left_target, right_target, du, left_du, right_du;
	char decision;
	int i=0;
	  left_du=du0_left;
	  right_du=du0_right;
      right_target=ShootFunc(right_du);
      left_target=ShootFunc(left_du);
	  if (right_target<left_target) { 
	     du=left_du;
		 left_du=right_du;
         right_du=du;
         target=right_target;
		 right_target=left_target;
		 left_target=target;
	  }
	  if ((right_target>=u1) && (left_target<=u1)) {// We have to ensure that the interval [du0_left, du0_right] contains the tangent of unknown solution
		  du=(right_du+left_du)/2;
		  target=ShootFunc(du);
		  decision=13;
		  while ((int) decision ==13){
			  i++;
			  if (target>=u1) {
				  right_du=du; right_target=target;
			  }  else {
				  left_du=du; left_target=target;	
			  };
		      du=(right_du+left_du)/2;
		      target=ShootFunc(du);
			  cout<<"Target="<<target<<endl;
			  cout<<"Press Enter to make one more shot"<<endl;
			  cout<<"or press any other key to proceed to the Multiple Shooting Procedure."<<endl;
              decision=getch();
		  }
	  } else {
		  cout<<"Invalid interval"<<endl;
		  cout<<"T0="<<left_target<<"; T1="<<right_target<<endl;
		  return 0;
	  }
      
      C(0)=bpU[0][1];
	  for (i=1; i<n; i++){
      C(2*i-1)=bpU[i][0];
      C(2*i)=bpU[i][1];
	  }
      C(2*n-1)=bpU[n][1];
      return i;
  }

//================================================================================
  // This function implements a multiple shooting procedure for solving the basic problem
  // Vector C0 is an initial approximation for the Newton method used to solve a nonlinear 
//================================================================================
  int T_BVP_PC_T::MultipleShootingProcedure( VectorXmp& C0){
	  if (C0.rows()!=2*n) {throw MSPBadInitData(); return -1;} // incompatiple initial data

      VectorXmp E(2*n);//Vector of the equations evaluated in C_i 
      VectorXmp C1(2*n);//Vector for storing interim values of the Newton's method  
	  SMatrixXmp J(2*n, 2*n); // The Jacobi matrix
	  BiCGSTAB<SMatrixXmp > solver; // Our BiCGStab solver;
	  solver.setMaxIterations(100000);
	  //solver.setTolerance(10e-40);
	  mpreal SqN_i, dN_i,SqN_i1, dN_i1, ShSqN_i, ChSqN_i, eps; // Auxiliary variables
	  eps=machine_epsilon(mpreal::get_default_prec());
	  int IterCounter=0;
	  int i;
	  C1(0)=1;

	  while (C1.norm()>eps*10e5){
         IterCounter++;
	     if (n==1){ // if there i only one subiterval
//===========================
            SqN_i=sqrt(F(u0));
            ShSqN_i=sinh(h_array[0]*SqN_i);
		    ChSqN_i=cosh(h_array[0]*SqN_i);
//===========================
  // Vector
//===========================
            E(0)=u0*ChSqN_i+C0(0)*ShSqN_i-u1; 
            E(1)=SqN_i*(C0(0)*ChSqN_i+u0*ShSqN_i)-C0(1); 
//===========================
  // Matrix
//===========================
            J.coeffRef(0,0)=ShSqN_i;
            J.coeffRef(1,0)=SqN_i*ChSqN_i;
            J.coeffRef(0,1)=0;
            J.coeffRef(1,1)=-1;
//===========================
//===========================
	     } else { // there is more than one subintervals
//===========================
//===========================
            SqN_i=sqrt(F(u0));
            SqN_i1=sqrt(F(C0(1)));
		    dN_i1=dF(C0(1));
            ShSqN_i=sinh(h_array[0]*SqN_i);
		    ChSqN_i=cosh(h_array[0]*SqN_i);
//===========================
  // Vector (The first block)
//===========================
            E(0)=u0*ChSqN_i+C0(0)*ShSqN_i-C0(1); 
            E(1)=SqN_i*(C0(0)*ChSqN_i+u0*ShSqN_i)-sqrt(F(C0(1)))*C0(2); 
//===========================
  // Matrix (the first block)
//===========================
            J.coeffRef(0,0)=ShSqN_i;
            J.coeffRef(1,0)=SqN_i*ChSqN_i;
            J.coeffRef(1,1)=-C0(2)*dN_i1/(2*SqN_i1);
            J.coeffRef(0,1)=-1;
            J.coeffRef(1,2)=-SqN_i1;
            
			for (i=1; i<n-1; i++){ // Here we fill the all other blocks (except for the last one)

                SqN_i=SqN_i1;
                SqN_i1=sqrt(F(C0(2*i+1)));
                dN_i=dN_i1; 
		        dN_i1=dF(C0(2*i+1));
                ShSqN_i=sinh(h_array[i]*SqN_i);
		        ChSqN_i=cosh(h_array[i]*SqN_i);
			    // The vector
                E(2*i)=C0(2*i-1)*ChSqN_i+C0(2*i)*ShSqN_i-C0(2*i+1);
                E(2*i+1)=SqN_i*(C0(2*i)*ChSqN_i+C0(2*i-1)*ShSqN_i)-SqN_i1*C0(2*i+2);
				// The matrix
				J.coeffRef(2*i, 2*i-1)=h_array[i]*dN_i*(C0(2*i)*ChSqN_i+C0(2*i-1)*ShSqN_i)/(2*SqN_i)+ChSqN_i;
                J.coeffRef(2*i, 2*i)=ShSqN_i;
                J.coeffRef(2*i, 2*i+1)=-1;  
                J.coeffRef(2*i+1, 2*i-1)=dN_i*(C0(2*i)*ChSqN_i+C0(2*i-1)*ShSqN_i)/(2*SqN_i)+
					            h_array[i]*dN_i*(C0(2*i-1)*ChSqN_i+C0(2*i)*ShSqN_i)/2+SqN_i*ShSqN_i;
				J.coeffRef(2*i+1, 2*i)=SqN_i*ChSqN_i;
                J.coeffRef(2*i+1, 2*i+1)=-C0(2*i+2)*dN_i1/(2*SqN_i1);
                J.coeffRef(2*i+1, 2*i+2)=-SqN_i1;
			} // End For
//===========================
			i=2*(n-1);
            SqN_i=SqN_i1;
            dN_i=dN_i1; 
            ShSqN_i=sinh(h_array[n-1]*SqN_i);
		    ChSqN_i=cosh(h_array[n-1]*SqN_i);
//===========================
  // Vector (The last block)
//===========================
            E(i)=C0(i-1)*ChSqN_i+C0(i)*ShSqN_i-u1;
            E(i+1)=SqN_i*(C0(i)*ChSqN_i+C0(i-1)*ShSqN_i)-C0(i+1);
//===========================
  // Matrix (the last block)
//===========================
			J.coeffRef(i, i-1)=h_array[n-1]*dN_i*(C0(i)*ChSqN_i+C0(i-1)*ShSqN_i)/(2*SqN_i)+ChSqN_i;
            J.coeffRef(i, i)=ShSqN_i;
            J.coeffRef(i+1, i-1)=dN_i*(C0(i)*ChSqN_i+C0(i-1)*ShSqN_i)/(2*SqN_i)+
				            h_array[n-1]*dN_i*(C0(i-1)*ChSqN_i+C0(i)*ShSqN_i)/2+SqN_i*ShSqN_i;
			J.coeffRef(i+1, i)=SqN_i*ChSqN_i;
            J.coeffRef(i+1, i+1)=-1;
	     }// End of If

/**         cout<<"E=\n"<<E<<endl;
		 cout<<endl;
         cout<<"J=\n"<<J<<endl;
		 cout<<endl;
		 cin>>i;
**/
//============================
   // Here we try to solve the problem
//============================
         solver.compute(J);// Initialization of the solver
         C1=solver.solve(E);
		 cout<<"Error=\n"<<solver.error()<<endl;
		 cout<<"Iterations=\n"<<solver.iterations()<<endl;

		 if (solver.info()!=Success) {// the process is devergent 
		     return -1;
	     }
/**		 cout<<"C0(old)=\n"<<C0<<endl;
		 cout<<endl;
**/
		 C0=C0-C1;
		 cout<<"C1.norm=\n"<<C1.norm()<<endl;
		 cout<<endl;
/**		 cout<<"C1=\n"<<C1<<endl;
		 cout<<endl;
		 cout<<"C0=\n"<<C0<<endl;
		 cout<<endl;
		 cout<<endl;
		 cout<<"Iteration :"<<IterCounter<<endl;
**/
//============================
	  } // End of while
	      //Prepareing the data for further usage 
          bpU[0][0]=u0;
          bpU[0][1]=C0(0);
		  bpU[0][2]=sqrt(F(bpU[0][0]));
	  for (i=1; i<n; i++){
          bpU[i][0]=C0(2*i-1);
          bpU[i][1]=C0(2*i);
		  bpU[i][2]=sqrt(F(bpU[i][0]));
	  }
          bpU[n][0]=u1;
          bpU[n][1]=C0(2*n-1);
		  bpU[n][2]=sqrt(F(bpU[n][0]));
	  return IterCounter;
  } 
  int T_BVP_PC_T::BasicProblemSolve(){ // The procedure that solves the basic problem using Single and Multiple Shooting Procedures
	  int res;
      mpreal left_du0, right_du0;
	  res=0;
	  while (res==0){
		  cout<<"Please enter the value of left point of the interval (du0, du1)"<<endl;
		  cout<<"du0=";
		  cin>>left_du0;
		  cout<<"Please enter the value of right point of the interval (du0, du1)"<<endl;
		  cout<<"du1=";
		  cin>>right_du0;
		  res=AdjustmentOfFireProcedure(left_du0, right_du0);
	  }
	  cout<<"C=\n"<<C<<endl;
	  res=MultipleShootingProcedure(C);
	  cout<<res<<endl;
	  return res;
  }
  int T_BVP_PC_T::BasicProblemSolveSimple(){ // The procedure that solves the basic problem using only Single Shooting Procedure
	  int res;
      mpreal left_du0, right_du0;
	  res=0;
	  while (res==0){
		  cout<<"Please enter the value of left point of the interval (du0, du1)"<<endl;
		  cout<<"du0=";
		  cin>>left_du0;
		  cout<<"Please enter the value of right point of the interval (du0, du1)"<<endl;
		  cout<<"du1=";
		  cin>>right_du0;
		  res=BisectionProcedure(left_du0, right_du0);
	  }
	  return res;
  }

  mpfr::mpreal T_BVP_PC_T::F01_1div3(const mpfr::mpreal& z){// This is the hypergeometric function 0F1(;1/3;z);
      int i;
	  mpreal sum_num, sum_denom, tz;
	  if ((z<=-600) || (z>=0)){ // This case can be successfuly handled by the function F01
          return F01(_1div3, z);
	  } else {// And this case we treat with the Chebyshev - Pade approximation
		  int k;
		  k=-z.toLong()/50;    
		  sum_num=0;
		  sum_denom=0;
		  for (i=Cheb_Cof_N-1; i>0; i--){
		      sum_num=(sum_num+MCH_1div3[k](i,0))*(z/25+2*k+1);
		      sum_denom=(sum_denom+MCH_1div3[k](i,1))*(z/25+2*k+1);
		  } 
		  sum_num=sum_num+MCH_1div3[k](0,0);
		  sum_denom=sum_denom+MCH_1div3[k](0,1);
	  }
		  return sum_num/sum_denom;
  }
                                                     
  mpfr::mpreal T_BVP_PC_T::F01_2div3(const mpfr::mpreal& z){// This is the hypergeometric function 0F1(;2/3;z);
      int i;
	  mpreal sum_num, sum_denom, tz;
	  if ((z<=-600) || (z>=0)){ // This case can be successfuly handled by the function F01
          return F01(_2div3, z);
	  } else {// And this case we treat with the Chebyshev - Pade approximation
		  int k;
		  k=-z.toLong()/50;    
		  sum_num=0;
		  sum_denom=0;
		  for (i=Cheb_Cof_N-1; i>0; i--){
		      sum_num=(sum_num+MCH_2div3[k](i,0))*(z/25+2*k+1);
		      sum_denom=(sum_denom+MCH_2div3[k](i,1))*(z/25+2*k+1);
		  } 
		  sum_num=sum_num+MCH_2div3[k](0,0);
		  sum_denom=sum_denom+MCH_2div3[k](0,1);
	  }
		  return sum_num/sum_denom;
  }
  mpfr::mpreal T_BVP_PC_T::F01_4div3(const mpfr::mpreal& z){// This is the hypergeometric function 0F1(;4/3;z);
      int i;
	  mpreal sum_num, sum_denom, tz;
	  if ((z<=-600) || (z>=0)){ // This case can be successfuly handled by the function F01
          return F01(_4div3, z);
	  } else {// And this case we treat with the Chebyshev - Pade approximation
		  int k;
		  k=-z.toLong()/50;    
		  sum_num=0;
		  sum_denom=0;
		  for (i=Cheb_Cof_N-1; i>0; i--){
		      sum_num=(sum_num+MCH_4div3[k](i,0))*(z/25+2*k+1);
		      sum_denom=(sum_denom+MCH_4div3[k](i,1))*(z/25+2*k+1);
		  } 
		  sum_num=sum_num+MCH_4div3[k](0,0);
		  sum_denom=sum_denom+MCH_4div3[k](0,1);
	  }
		  return sum_num/sum_denom;
  }
  mpfr::mpreal T_BVP_PC_T::F01_5div3(const mpfr::mpreal& z){// This is the hypergeometric function 0F1(;5/3;z);
      int i;
	  mpreal sum_num, sum_denom, tz;
	  if ((z<=-600) || (z>=0)){ // This case can be successfuly handled by the function F01
          return F01(_5div3, z);
	  } else {// And this case we treat with the Chebyshev - Pade approximation
		  int k;
		  k=-z.toLong()/50;    
		  sum_num=0;
		  sum_denom=0;
		  for (i=Cheb_Cof_N-1; i>0; i--){
		      sum_num=(sum_num+MCH_5div3[k](i,0))*(z/25+2*k+1);
		      sum_denom=(sum_denom+MCH_5div3[k](i,1))*(z/25+2*k+1);
		  } 
		  sum_num=sum_num+MCH_5div3[k](0,0);
		  sum_denom=sum_denom+MCH_5div3[k](0,1);
	  }
		  return sum_num/sum_denom;
  }

  mpfr::mpreal T_BVP_PC_T::F01(const mp_pair& a1, const mpfr::mpreal& z){ //This is the hypergeometric function 0F1
	  // This function will be used for computing the A and B functions. 
	  //It is supposed that a1>0 and there can not be any troubles  with division by 0.
	  mpreal old_sum, new_sum, sum, r, sqrt_z, K1, K2, sum1, sum2, old_sum1, old_sum2, t1, t2, p1, p2, b1, b2, b3, b4, b5, b6;
	  int i=1;
/**	  // Simle Tailor series method (works bad for big arguments) 
	  if (abs(z)<=1000) {
	      sum="1";
	      old_sum="0";
	      r="1";
	      while (sum!=old_sum){
		      r=z*r/((a1.value+i-1)*i);
		      old_sum=sum;
		      sum=sum+r;
		      i++;
		      //cin>>eps;
	      }
	  } else 
**/	  
	  if (abs(z)<500)  {
          r=1/a1.value;
          old_sum=1+z*r;
	      r=1/(2*(a1.value+1));
	      new_sum=old_sum+sqr(z)/a1.value*r;
	      i=3; 
	      while ((abs(old_sum-new_sum)>0)){
		      r=1/(i*(a1.value+i-1));
		      sum=new_sum+(new_sum-old_sum)*r*z;
		      old_sum=new_sum;
		      new_sum=sum;
		      i++;
//			  cout<<"sum="<<sum<<endl;
//			  cout<<"i="<<i<<endl;
	      }
	  } else if (z>=500) {
          sqrt_z=sqrt(z);
		  K1=pow(sqrt_z, "0.5"-a1.value)*a1.gamma/(2*sqrt(const_pi()));
		  K2=K1;
		  K1=K1*exp(2*sqrt_z);
		  K2=K2*exp(-2*sqrt_z)*cos(const_pi()*("0.5"-a1.value));
//		  cout<<"Gamma="<<a1.gamma<<endl;
//		  cout<<"K1="<<K1<<endl;
//		  cout<<"K2="<<K2<<endl;
		  b1=(a1.value-"0.5");
		  b2=("1.5"-a1.value);

		  p1=b1*b2/(4*sqrt_z);
//		  cout<<"i="<<i<<" p1="<<p1<<endl;
		  t1=p1;
		  t2=-p1;
		  sum1=1+t1;
		  sum2=1+t2;
		  old_sum1=0;
		  old_sum2=0;
          i=1;
		  p1=(b1+i)*(b2+i)/(8*sqrt_z);
//		  cout<<"i="<<i<<" p1="<<p1<<endl;
		  while ((abs(p1)<1) && ((old_sum1!=sum1) || (old_sum2!=sum2))) {     
			  t1=t1*p1;
			  t2=-t2*p1;
		      old_sum1=sum1;
		      old_sum2=sum2;
		      sum1=sum1+t1;
		      sum2=sum2+t2;
			  i++;
              p1=(b1+i)*(b2+i)/(4*sqrt_z*(i+1));
//			  cout<<"i="<<i<<" p1="<<p1<<endl;
//			  cout<<" sum1-old_sum1="<<sum1-old_sum1<<endl;/
//			  cout<<" sum2-old_sum2="<<sum2-old_sum2<<endl;
//			  cout<<" sum1="<<sum1<<endl;
//			  cout<<" sum2="<<sum2<<endl;
		  }
		  return K1*sum1+K2*sum2;
	  } else if (z<=-500) {
          sqrt_z=sqrt(abs(z));
		  r=(2*a1.value-1)*const_pi()/4-2*sqrt_z;
		  K1=a1.gamma*pow(abs(z), (1-2*a1.value)/4)/sqrt(const_pi());
		  K2=K1;
		  K1=K1*cos(r);
		  K2=K2*sin(r)*(2*a1.value-1)*(2*a1.value-3)/(16*sqrt_z);
		  b6=a1.value/2;
		  b1="0.75"-b6;
		  b2="1.25"-b6;
		  b3=b6-"0.25";
		  b4=b6+"0.25";
		  b5="1.75"-b6;
		  b6=b6+"0.75";
/**		  cout<<"b1="<<b1<<endl;
		  cout<<"b2="<<b2<<endl;
		  cout<<"b3="<<b3<<endl;
		  cout<<"b4="<<b4<<endl;
		  cout<<"b5="<<b5<<endl;
		  cout<<"b6="<<b6<<endl;
**/		  i=0;
		  r=b2*b4/(4*z);
		  p1=r*b1*b3/"0.5";
		  p2=r*b5*b6/"1.5";

		  t1=1;
		  t2=1;
		  sum1=1;
		  sum2=1;
		  old_sum1=0;
		  old_sum2=0;
		  while ((abs(p1)<1) && (abs(p2)<1)  && ((old_sum1!=sum1) || (old_sum2!=sum2))) { 
			  t1=t1*p1;
			  t2=t2*p2;
			  old_sum1=sum1;
			  old_sum2=sum2;
			  sum1=sum1+t1;
			  sum2=sum2+t2;
			  i++;
		      r=(b2+i)*(b4+i)/(4*z*(i+1));
		      p1=2*r*(b1+i)*(b3+i)/(1+2*i);
		      p2=2*r*(b5+i)*(b6+i)/(3+2*i);
/**			  cout<<"i="<<i<<" r="<<r<<endl;
			  cout<<"i="<<i<<" p1="<<p1<<endl;
			  cout<<"i="<<i<<" p2="<<p2<<endl;
			  cout<<"i="<<i<<" sum1="<<sum1<<endl;
			  cout<<"i="<<i<<" sum2="<<sum2<<endl;
			  cin>>sum;
**/		  }    
		  return K1*sum1+K2*sum2;
	  }
//	  cout<<"i="<<i<<endl;
	  return sum;
/**
	  mpreal h_old, h_new, sum, old_sum, term, a; // Continued fractions (works not very good regarding the accuracy)
	  int i=2;
      h_new=1+z/(2*(a1+1));
      term=(1/h_new-1);
      sum=term;
      old_sum=0;
	  while (sum!=old_sum) {
          h_old=h_new;
		  a=z/((i+1)*(i+a1));
          h_new=1+a-a/h_old;
          term=term*a/(h_old*h_new);
          old_sum=sum;
          sum=sum+term;
		  i++;
	  }
      sum=1+z/(a1*(1+sum));
	  cout<<"i="<<i<<endl;
	  return sum;
**/
  }
  int T_BVP_PC_T::GetMatchMatrix(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z, Matrix2x2mp& M) { // The procedure that returns the "match matrix" 
      mpreal w, v, y;
	  if (abs(a)>=Tan_threshold) {
		  y=(a*z+b)/cbrt(sqr(a));
          w=y*sqr(y)/9;
		  v=(sqr(a*z+b)/(2*a));
		  M(0,0)=F01_2div3(w);
		  M(0,1)=y*F01_4div3(w);
		  M(1,0)=v*F01_5div3(w);
		  M(1,1)=cbrt(a)*F01_1div3(w);
	  } else if ((abs(a)<Tan_threshold) && (b>0)) {
		  v=sqrt(b);
          w=v*z;
		  M(0,0)=sinh(w);
		  M(0,1)=cosh(w);
		  M(1,0)=v*cosh(w);
		  M(1,1)=v*sinh(w);
	  } else if ((abs(a)<Tan_threshold) && (b<0)) {
		  v=sqrt(abs(b));
          w=v*z;
		  M(0,0)=sin(w);
		  M(0,1)=cos(w);
		  M(1,0)=v*cos(w);
		  M(1,1)=-v*sin(w);
	  } else {
		  M(0,0)=1;
		  M(0,1)=z;
		  M(1,0)=0;
		  M(1,1)=1;
	  }
	  return 0;
  }

  mpfr::mpreal T_BVP_PC_T::A(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z){ // This is the A-function 
	  if (abs(a)>=Tan_threshold) {
		  return F01_2div3(sqr(a*z+b)*(a*z+b)/(9*sqr(a)));
	  } else if ((abs(a)<Tan_threshold) && (b>0)) {
		  return sinh(sqrt(b)*z);
	  } else if ((abs(a)<Tan_threshold) && (b<0)) {
          return sin(sqrt(abs(b))*z);
	  } else {
	      return 1;
	  }
  }
  mpfr::mpreal T_BVP_PC_T::dA(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z){ // This is the derivative of A-function 
	  if (abs(a)>=Tan_threshold) {
		  return (sqr(a*z+b)/(2*a))*F01_5div3(sqr(a*z+b)*(a*z+b)/(9*sqr(a)));
	  } else if ((abs(a)<Tan_threshold) && (b>0)) {
		  return sqrt(b)*cosh(sqrt(b)*z);
	  } else if ((abs(a)<Tan_threshold) && (b<0)) {
          return sqrt(abs(b))*cos(sqrt(abs(b))*z);
	  } else {
	      return 0;
	  }
  }
  mpfr::mpreal T_BVP_PC_T::B(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z){ // This is the B-function 
	  if (abs(a)>=Tan_threshold) {
		  return (a*z+b)*F01_4div3(sqr(a*z+b)*(a*z+b)/(9*sqr(a)))/cbrt(sqr(a));
	  } else if ((abs(a)<Tan_threshold) && (b>0)) {
		  return cosh(sqrt(b)*z);
	  } else if ((abs(a)<Tan_threshold) && (b<0)) {
          return cos(sqrt(abs(b))*z);
	  } else {
	      return z;
	  }
  }
  mpfr::mpreal T_BVP_PC_T::dB(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& z){ // This is the derivative of B-function 
	  if (abs(a)>=Tan_threshold) {
		  return cbrt(a)*F01_1div3(sqr(a*z+b)*(a*z+b)/(9*sqr(a)));
	  } else if ((abs(a)<Tan_threshold) && (b>0)) {
		  return sqrt(b)*sinh(sqrt(b)*z);
	  } else if ((abs(a)<Tan_threshold) && (b<0)) {
          return -sqrt(abs(b))*sin(sqrt(abs(b))*z);
	  } else {
	      return 1;
	  }
  }
  int T_BVP_PC_T::SaveToFile(const std::string& filename){
	  ofstream outfile(&filename[0]);
	  mpreal x;
	  mpreal prec=trunc(mpreal::get_default_prec()*0.3);
	  outfile.precision(prec.toULong());
	  x=0;
	  for (int i=0 ; i<n; i++){
	    outfile<< x<< ' '<< bpU[i][0]<<endl;
		x=x+h_array[i];
	  }
	    outfile<< x<< ' '<< bpU[n][0]<<endl;
    return 0;
  } 
  int T_BVP_PC_T::SaveToFile(const std::string& filename, const int& substeps){
	  ofstream outfile(&filename[0]);
	  mpreal sub_h, x, u;
	  mpreal prec=trunc(mpreal::get_default_prec()*0.3);
	  outfile.precision(prec.toULong());
	  x=0;
	  for (int i=0 ; i<n; i++){
          sub_h=h_array[i]/substeps;
		  for (int j=0; j<substeps; j++){
			  u=bpU[i][1]*sinh(bpU[i][2]*sub_h*j)+bpU[i][0]*cosh(bpU[i][2]*sub_h*j);
	          outfile<<x+j*sub_h<<' '<<u<<endl;
		  }
		  x=x+h_array[i];
	  }
	      outfile<<x<<' '<<bpU[n][0]<<endl;
  return 0;
  } 

T_BVP_PC_T::~T_BVP_PC_T(){
	for (int i=0; i<=n; i++){
	  delete [] bpU[i];
	}
	delete [] bpU;
	delete [] bpUD;
	delete [] h_array;
}

int T_BVP_PC_T::LoadChebPadeCoeffs(){
  string number;
  ifstream myfile;
  int i,k;
	  // Fill arrays for Chebyshev coefficients
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  // 0F1(1/3, z);
  //===================================================================================================================================================== 
  //=====================================================================================================================================================
  myfile.open("_1div3.txt");
  for (k=0; k<Cheb_Intervals_N; k++){
	  MCH_1div3[k].resize(T_BVP_PC_T::Cheb_Cof_N, 2);// [-50*(k+1); -50*k]
	  for (i=0; i<T_BVP_PC_T::Cheb_Cof_N; i++) {
  	     myfile>>number; MCH_1div3[k](i,0)=number;
	     myfile>>number; MCH_1div3[k](i,1)=number;
	  }
  }
  myfile.close();
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  // 0F1(2/3, z);
  //===================================================================================================================================================== 
  //=====================================================================================================================================================
  myfile.open("_2div3.txt");
  for (k=0; k<Cheb_Intervals_N; k++){
	  MCH_2div3[k].resize(T_BVP_PC_T::Cheb_Cof_N, 2);// [-50*(k+1); -50*k]
	  for (i=0; i<T_BVP_PC_T::Cheb_Cof_N; i++) {
  	     myfile>>number; MCH_2div3[k](i,0)=number;
	     myfile>>number; MCH_2div3[k](i,1)=number;
	  }
  }
  myfile.close();
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  // 0F1(4/3, z);
  //===================================================================================================================================================== 
  //=====================================================================================================================================================
  myfile.open("_4div3.txt");
  for (k=0; k<Cheb_Intervals_N; k++){
	  MCH_4div3[k].resize(T_BVP_PC_T::Cheb_Cof_N, 2);// [-50*(k+1); -50*k]
	  for (i=0; i<T_BVP_PC_T::Cheb_Cof_N; i++) {
  	     myfile>>number; MCH_4div3[k](i,0)=number;
	     myfile>>number; MCH_4div3[k](i,1)=number;
	  }
  }
  myfile.close();
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  // 0F1(5/3, z);
  //===================================================================================================================================================== 
  //=====================================================================================================================================================
  myfile.open("_5div3.txt");
  for (k=0; k<Cheb_Intervals_N; k++){
	  MCH_5div3[k].resize(T_BVP_PC_T::Cheb_Cof_N, 2);// [-50*(k+1); -50*k]
	  for (i=0; i<T_BVP_PC_T::Cheb_Cof_N; i++) {
  	     myfile>>number; MCH_5div3[k](i,0)=number;
	     myfile>>number; MCH_5div3[k](i,1)=number;
	  }
  }
  myfile.close();
/**  MCH_1div3[0].resize(Cheb_Cof_N, 2);// [-50; 0]
  MCH_1div3[0]<<"-1.2442345400630397890686023328846678986421047516922584778074169121381899837068922",	"0.94984173413594851388283178179068234147935800297567330723074235423561891073616614",
"-6.3584401017025397011651957816925495408302845725751467021434021661442164820604398",	"-0.4324598839328199141801155592983248970811623554904274100246561637358706823666503",
"17.421996809419671587168064911749652877910797914213949629258016058697689241820122",	"0.098997017400202497905583135314089117918573713987750144904306870026269214605638015",
"27.995765831064321545493226270292993296338055338462234899183220088625952367387664",	"-0.015172906760067312295173260780794361038301550164339237249910387344122560441952144",
"-23.874693602763736413326515212301860141098342948905416285414228849430181938915688",	"0.0017489815519978873131877049288270235167614199750038772068008183350096789681454497",
"-42.62882117624970280430612163519230711539802746827485046818155897679920848669177",	"-0.00016144392434391981834975536397800302606235069018346589553346796760729025272242448",
"-7.0989188460389056227356569241394317716463408241789183606042751072134629888602018",	"1.2404376554234381548947322871939532572259247055542563085181971839215485567131256e-05",
"16.197760064629934012477749735122487543013779031519235543644410953330695338435046",	"-8.1387389704980235611558640252574429667374122167690806692678113665642019700329500e-07",
"13.312746843589541794132380758067408105981744344031077336301883265455896940605818",	"4.6404753198069519938826064154254252037360195451017734486261274350771918450629762e-08",
"5.2835745428213738265675608870824588009530545346425685549147253983806845618126367",	"-2.3268432229255951435294982455013656907336148930418549892513405714089055185810683e-09",
"1.3307479545316751345601927281690524300137294206220969267162961437071502629337552",	"1.0339806584136596228964988376066408405481565600709654737771361101685509716110585e-10",
"0.23280321344892742580356522729802026941912043375866653854395418912284755303667916",	"-4.0888056871050245896049889557140474239713725671092645925453644744877971673485451e-12",
"0.02953370306152901095624775325798051596556141687059341059455036182922478226783828",	"1.4400173670059727227357643418619549720999749188460244776906716826429354770006177e-13",
"0.0027781134911887279937971646700388500688519715978802919747258460995251157690675413",	"-4.5044984081057248332728185992370208202124124593210607765980503359997690076677336e-15",
"0.00019561213405138521388400016987108996528677301444677592983088819101512243182676051",	"1.2429350378878349046226613379289378294173228123787796767389941507828081228550874e-16",
"1.0299708794461324778373489266300376210087441403523713048306860990016585834871367e-05",	"-2.9882571935600508457620145414336780575343668546886892641392418673986921240504354e-18",
"4.0064779198762391257215988892120064423685553294485119407582384550104636211013202e-07",	"6.1340194437248714690903228650647469827219297961114254008214681812144352660180028e-20",
"1.1203880561957288795253380062762118662115620880843023438437765123655638251800880e-08",	"-1.0399634621291356197348874905037725755912969030745386611256157623497481439580648e-21",
"2.1366691811620361428286163850642315415331163026646187561539627228598066186174469e-10",	"1.3755446790896169575827547161751229432493233555440567816451804806635514947480439e-23",
"2.4962150491317089206687565282629791832910134491080771751587923575266054055200286e-12",	"-1.2712481701687678620540059267704599575598320224212872392528070750975053687892641e-25",
"1.3527453302944236729599121197021028983787876709979485319500673955287026360437482e-14",	"6.2038777291079373128436003022055657334079179099132446626470447857987421439157361e-28";


  MCH_1div3[1].resize(Cheb_Cof_N, 2);// [-100; -50]
 MCH_1div3[1]<< "0.6236249824698829824548648375074806493487145258480760578207591948942107728498",	"0.94619045192727972037742292030098003300580688917749346468031941898389384398135467",
"-5.94616129889963066936971930338729878151804243995144154619902672085048532751",	"-0.44639100610841922569246689798582707151795492162261736608192608163292767672927192",
"-0.1760548722509456982987621096270754003295969850026305099696596373962005904713681",	"0.1060846630415737422728941517764909650437054297856070092464874456041739330586589",
"8.1476465782017264411944775304499116936353226108790419887223807579105408821412",	"-0.016914027645849517290870866452474737446635031178538523760356551671937599797207286",
"-0.25901852769610864324180793597559246715643162107992755879708527385927917057710397",	"0.0020327061633293037139055313172296207245089999393679486895250821761196936817339959",
"-3.0842010288990881885805856114434375022847892740358234599875749995111856332716291",	"-0.00019609731300489220760396900295742813518026130207739714834898605726286910467281657",
"-0.1754585276436071249684841806801451351088780988663399165521364269284167988021",	"1.5788182483975188322868134082005906749901922132710993818904242675022845856600247e-05",
"0.48614350681109715562427985366771960726297042736866680938657851100509506275549",	"-1.0886261074055126971911331248527556755337392048327776750525860874071434317852019e-06",
"0.087754327037998820858896842611317722966184374309567203562002413928952900175637",	"6.5438276925122081123670073121363880828184066339052831550376426223730896761960914e-08",
"-0.02788483288372415800522496877645130414940444514584453573102093808411833150539468",	"-3.4714630549254252109358574576343921504133014778050516266745080173776378474253252e-09",
"-0.010555723752838238153451949846462580832988688930592152354534127929931095668201987",	"1.6384370702669475725387127030777486511878549830339999057747726636020174851424146e-10",
"-0.00059999459958106261465138757357499327560925524505460357558759261610464512292977",	"-6.9115727410871141090008810792937892703598691318216276277495958500286168887137503e-12",
"0.0002861388380673765918618042952089493793862375560066752487341751679328997922114022",	"2.6093456829701620932516671043055293278448466535241703494376519936380127401105856e-13",
"7.739154891691701430210762983282096244342917970962526928128731764117902460468430e-05",	"-8.7980396347608892356284301195714976984802080333353845253580064671390435651467303e-15",
"9.981893662651962998060433391420264412134362850912448548273894882771741487970570e-06",	"2.6331819265802707547632258741739856360737874339505236325752892766108251854844085e-16",
"8.076148045084769416404535359804926681575379465765200730132644258686611637911855e-07",	"-6.9159949024311179171623401486402329581303714626710659287997583940804747288895928e-18",
"4.379741640450494332952856077648142799385401981973122350421796112736149921408918e-08",	"1.5637988986685868082532168101205188395038596350814706266679249211868254718143946e-19",
"1.6038218328811744607114955003071409559230946133556696684949125030422009036382688e-09",	"-2.9488351119693144984364693591514441762249526839819582452879256719595589793926678e-21",
"3.8330647717932489950711540943510178433915858123442908304176933885551245951747686e-11",	"4.3880024235329465863550543722469302734844644086507880083643082211725566305724090e-23",
"5.4320225710585728249066588853741678897448512832885437579837680021219952058776491e-13",	"-4.6252971968842798854931976296185319830464687803363508327541579424560358705577303e-25",
"3.4821662252880642857289617642106502481488686843196046237217390753802325626720614e-15",	"2.6180370855949290693446292835875521083342243943434846077417320257439862040089694e-27";

  MCH_1div3[2].resize(Cheb_Cof_N, 2);// [-150; -100]
 MCH_1div3[2]<<"-1.713820231507864097734689261367012349597135479818090271707275894775501678156",  ".94226750017499859033584930659748302921878659637352055856454069440793960742428766",
"-1.960584397840385352062007041424350926796569311584978718624668655893251030",  "-.46068612031161267678431831698681877852219049178600966610184381508788439202869350",
"5.352721487006640926067736089337264810332648605110735736002644382795225926001653",  ".11367944367632918135150976681515832075731738244858537282732175391325281895287514",
".3521516425530051048719724983233374563274999417920134021432610729942004863",  "-.18859769585171486902602641437594180658996016324534899104094097916496615410180113e-1",
"-2.19578270866264424917612391380249468327138218395976252570321146771723073467",  ".23638991085646722976027903823789942072010812204716081152016139235978581501110466e-2",
".6843413619789201987260268910463480750626604940718651737137763234040109689248129e-1",  "-.23844721245349713624492037347740582672789815891021252839179938027475231707320244e-3",
".321532721733862886813371013001835643767913839979770854258265053160584547302744",  ".20129389163575148085475275439013274393614905432037821214475133573611484624925017e-4",
"-.5189584539464043938869215138056462295038569131679857664964794275192098335798e-2",  "-.14598014628334079234980915754492576263359597665358818328355355471529373062421672e-5",
"-.23099592900009509048093078215926722901823405235565463315651193206270952725010e-1",  ".92608309449873879237222846537695842455588931313196293221960828035638869836758866e-7",
"-.78191860627215103566920986867933736298505121919266786695191930765398595567e-3",  "-.52046703908650101678521363764046890043491126514116146252712675260272629448659348e-8",
".8664935920918376731762571422741720008900646203401657622559536946190788007484380e-3",  ".26135888684212240675189587497469122732133626699283394970803407232647578857416955e-9",
".8522351049637305106827078597567625447368919866272207810044865718857384682782036e-4",  "-.11787460007415756665368588167640713239779734547784784465330105076518240173371192e-10",
"-.134377374150245635345941226622053183364949644986860779121939823779492993030489e-4",  ".47842756669610587699203278816159374256665967054144234189881277284108463168799540e-12",
"-.274097494582719650007888308363928681963091151125283684088705978900052369912771e-5",  "-.17453445139025066253485701460573861262421140205654529152966797150973274605577239e-13",
"-.6872466228660431865483390440585054276270765979115468692488637296237703140860e-7",  ".56939171854461758291944574330814329730192433834936190428603361518997490998444179e-15",
".234530516980357311833090460544480238771217140332351443569586377917080048684271e-7",  "-.16444645713369296191752771971490457496136762890529458806794263907808876591875340e-16",
".30556132323037942710607109979274687249821254361654386770546540273938072815231497e-8",  ".41319368496953104783569567301060319656251863352670513921430700415138076709555772e-18",
".18273595870236234553846696441941905984572780936728478325906826974985778163955562e-9",  "-.87702900918975095104308086885175179519135692699616345286762070090099088493832553e-20",
".61582678999126740558028178265934081874774854278288265961916444490658850293094985e-11",  ".14929713311387612643095459148374026909798499120363589681946603607415883221225251e-21",
".11368841808134765666106658998938383326839694029994592216945610663825970242698492e-12",  "-.18387912068936532607819388262829478678222585696607739670934530608444450872509631e-23",
".90305564811511491136197256428207701182482876660431136203612404333253197658621085e-15",  ".12523798158097489619387538461091239183478406690029646995388657606996122744981537e-25";

  MCH_1div3[3].resize(Cheb_Cof_N, 2);// [-200; -150]
  MCH_1div3[3]<<"-.424041847517697903813668635758123506188672338750152187364539031098425654e-1",  ".93808936140191771132896757935319503799213415091568434779868824998060582403360290",
"4.141638571740455845462811863016374392181068800365740818294256172798397724",  "-.47522264839472363781485391642865615804050700162062526778022847142818617318534580",
"-1.9188951606598526688519666339617158941893986395919183632007903816733320034857",  ".12174500232456760736083966329010905023348728599013854780528717793761320874316408",
"-1.99383051887999495654661655711408040485547222424002609241360298350111596028",  "-.21014370999721925407302675391899937344803273783326682719544775020520831707145995e-1",
".9117497424671552673306094474895284630593603252569984059067160382831815535",  ".27468995979672153178968843190651578182191914971578800611821669574554829271237397e-2",
".232275717325600536799246027395456461252260260632143282377186056480658568986",  "-.28971141331546285227703459723188839311543701930620712765978723914013215895067680e-3",
"-.114883367860219192845033152365928793214941059855591619073020035549676868502",  ".25645302927481138686027388622734603917649279396896356027523364599945960884689942e-4",
"-.1192975519602080932576530405355407552899619308069476103716391499415938684335e-1",  "-.19564002914708084453629664442523097752638742651005235736008523601058702939305627e-5",
".62438682534627366255685417982418302042117799910225113251369957275580031995235e-2",  ".13102263424314365330307671684909447399084692814579690986911454057310155270235057e-6",
".422544900141837163773862416144677524595127578376739757049805474190602782233e-3",  "-.78047984013410854160713965936439906420094180744378224363694191146740211761883831e-8",
"-.176430421344333007510393347213818621023928213825630422614628174197368364606e-3",  ".41730250296076795221034625242527017115580236165332061820680447584590446814606013e-9",
"-.126123963890764113925764442213483132922190577739794701849175411409285966729e-4",  "-.20143627245532862898732163806958358216564100682458860625720402512917146134946335e-10",
".265785452184926143334085245547357312450848863760724907725228097987232749577179e-5",  ".88033959963824148574205996862027181993235857619590524095766679318282811282718535e-12",
".26061463842974426545278108185195607070089720345732145528138058610964094168619e-6",  "-.34825206359439194942511182838044061212758787894425770089634311626927961499324823e-13",
"-.16733055721092070796205513400057710523343723739022800300032595565817290435134e-7",  ".12423585446945582958227279324711644970222200282319110585476932241506575657082696e-14",
"-.287825182429893705112781390205305489555238967756030627649869282261494431806609e-8",  "-.39637618033069953098619010767947970674764673284960094999019334909574739149274373e-16",
"-.43618910199073309512462456691189944220539522913533387941580571953539602845452e-10",  ".11142739368705530196682011640874679163827322682432177499489760213617658201305885e-17",
".1109232281912331194214812649874432409824544870757243154921364250543914825377666e-10",  "-.26896558083532167711449446574205629279102242973189715273121618037959996849146733e-19",
".7968484954947283684494282930426778649336401915341674968739577170317683707664293e-12",  ".53233683198972472689528154055039029282228677465477411338882701359026789278943690e-21",
".22290766436195880849452230555154058547955485340068998551735283490118035924712500e-13",  "-.78712530467116168662133921249521092776179574941246997904673212469840061449380626e-23",
".23648591606909358215108329868881142574381212109992766785465193410881570712123376e-15",  ".68219948410037867471713847028592871305738811707040552804837935181378005251264290e-25";


  MCH_1div3[4].resize(Cheb_Cof_N, 2);// [-250; -200]
  MCH_1div3[4]<<".903560522090199942702203369021620217985991303820782088502556980212608093",  ".93371643815033353156240119799284557449590534613329144380860076786141013974001680",
"-3.854913860343216604583471980167018489056119045369092811603525058473122020",  "-.48969542650045015136294989483152324805355073533530354028399942582061220294531362",
".58161049852649697456906240308246555571546538541516604641971554221457356",  ".13015947367260190009935069562362858121887230414453372988942278155545662408422955",
"1.6973118090116866428904865956041229132884948117647200271273225296850823266336",  "-.23361463973615665574703145026932024758845523373885923526643798843443776547287866e-1",
"-.48057400435155680450743773031091209247985320119718019911722897890374451451",  ".31829413532852386073191853015266072164422418647459010612607382881011345119462742e-2",
"-.15929093381103064980280796857406873028493644108653628664773695448150141260",  "-.35082569061101863016981825055805629696597292542496980941269991938760727125547057e-3",
".574249392899904179112743910354648556506192998103539452657101566321441398989e-1",  ".32548450500290608101990089639313183893151865814200406184032234513437398113840013e-4",
".5024633409847056014272468813481003722206285635670234716021673199695556334606e-2",  "-.26107903895141480793204835149121264625349962345776672427488757787555725179420948e-5",
"-.26414564846295520284655120552788782136058947797750233675638429544600435890273e-2",  ".18450881221045998731063171350479203193167524232477831842402618989502787587650175e-6",
"-.5788052083461731243809529169864838761915408501617960803974461338235211245e-4",  "-.1164537704946005844564657393226992427494085975927997895539755772355130233533406e-7",
".61773130256588008494619709642518586675233431899745652803398077021793974808406e-4",  ".66279538256100780042289860883898153528418132864814746510908491677564757615898039e-9",
".2601464875243775990313530300721617711362596477374177045679056516876675260e-6",  "-.3423887140750148460558904369999421640230522012215906776455430022867783621463383e-10",
"-.836745506198997452454947284572529298912657494894377134772513124908555043067e-6",  ".1611389676684719021933943649695986408823575464150015779409636437237139154795169e-11",
"-.600454842930399112404756791768121539604640850811461501339476908914950245027e-8",  "-.6915307378085279477301474172459057289722968021693670741709747926167056266123392e-13",
".68087264329923925105344410945767440284230338680149223053101761388180138522748e-8",  ".27004862716993772765847478629623198003540196700685456504981881325787386288872636e-14",
".142302077305001581914478548104222551321119053605331083018924891448687341370404e-9",  "-.9535200847770684036134990561030538276455539204285410867960586822683039596122290e-16",
"-.30627687055189410081718471583090844230298995633998916453034183591577757590762e-10",  ".3009902268609525118547200151821551758145521824916276319924482601164425283664601e-17",
"-.124737104320633914349427319009060964320144071779417931314751231743310428966957e-11",  "-.8307214424369069505795092587415112135350224668560628964626222582472876016930448e-19",
".4996817595258048408211017562260996824885617299649348116203566011431458942758165e-13",  ".19373227020542389339710494961488711523726216338326900203170091648644318001725088e-20",
".3816738421117363818835054323523061452370613564880427279583530283893758372249971e-14",  "-.34946495458639057145873596997701245522885775735765067593558197980112528191984320e-22",
".62422932406068650106011221087839724077950464505462905342092079120518898508563219e-16",  ".42415980284050228102797408162155006372723863903516695622995454648872376860854614e-24";


  MCH_1div3[5].resize(Cheb_Cof_N, 2);// [-300; -250]
  MCH_1div3[5]<<"-.98937701520962912086087332717759609167008172482573556574626748654790694",  ".9537312030054585146380001793563517256746239727142767683543092873293134730610731",
"3.56166042957178282245014550416867047517888864575213431661993440097854951",  "-.4057144254712403768337894226836563773225434937000475359795853804707228306568374",
"-.25886508947956989361255784860713944019680855383350697681626477559360740",  ".912546019067419817957990722197108133555349133511574153159459596387819615182685e-1",
"-1.325639055213522123724235896494010087617125210341530170733156102032999835166",  "-.14171707573539749059795886523243667976059840932589787808412383013169432690372538e-1",
".2671184241576560784830746061398253964113989536562808148595813858980762096",  ".16990017187125172388908497631707120987364177568100597890621694471705451348285303e-2",
".11192807587747374173856524576178765413750779165034997241382865584552486150",  "-.1668099112725451281323319775666035185538436210393819599008610224361071568928931e-3",
"-.27669634303835314592547482658785343711514878522395533658967317184939093252e-1",  ".1392779919604964147887579576632118624690031737034279832871397817724349786064241e-4",
"-.3567976021044949659287119135723754735163489992805841581712435610226598947426e-2",  "-.10137428760729559331793262577964899453915037554905213301774131386979531953294584e-5",
".10727458825478586256496128270755188723402081395946016527419909864329302688110e-2",  ".65499335844118680411647112101428074572540073265897481618266945873986211306802103e-7",
".5582105285463393062413377570795942496217067325871263661348567089378705282e-4",  "-.380365765965054366348685865493758926485263677011897310657175775516993926045097e-8",
"-.20761520633017051197092756288694231911067073295358356883466755070860496477570e-4",  ".20046746989310806750745937469669610093823524453020854503061399765973629356555549e-9",
"-.57922407605274718879036807182796924319418354818246045151944323490882117239e-6",  "-.96431136220974408369819877481711792201403298218456739670721926271506045267790e-11",
".227248017971008331920339079378641269119229881791274450500248862677657728658e-6",  ".425402223360278427386014331234053055079185465250097226321585946549579784681377e-12",
".563617040383835081295331288112262685324533136808747190712892956777569091030e-8",  "-.1720590572881756459944888348865385569019830403390055521034474514703436643458994e-13",
"-.145101097181830202603881654391440888628581960415586732590986315605869242984839e-8",  ".6387569916086838955074641158135502606052053598819501684104319646027797887273013e-15",
"-.47116038683924077001411709308150268912322092854105630613418947228018908850164e-10",  "-.2155942077002840805005142231133500296357144715757252025826109321832591568088640e-16",
".49489731148681759567035002441566128020879534351730517588377562034297314162384e-11",  ".660819462618117577333645481759008991516218594051579362180179915332909293754993e-18",
".236068695233823166394222403374584889243486763634654333824660435670892847688942e-12",  "-.1776321309104332182390983758882217120518056561636306357514046602092482841495723e-19",
"-.556634865568507909795714605349079445564455434161460297649764379502942546166518e-14",  ".4212668958798534772028909936985019247223362347921274829829540026199391549959279e-21",
"-.4909318936740993156003071928826238953055643500191103846407147301191969756511945e-15",  "-.7568831284018792396258152786656631044236455192316507325482679760217744959124046e-23",
"-.76976756645383384939906634134038484042756407135351038815874325177670320760660291e-17",  ".11382411681061565454002591375848111543343694681144343867426510332605395953488958e-24";


  MCH_1div3[6].resize(Cheb_Cof_N, 2);// [-350; -300]
  MCH_1div3[6]<<".4503463342435430415244660464618869290295208468157663080294461101993058",  ".96577914585656196254342013592585842193342897548830979492523612131907926253618",
"-3.3806103069966359576379308813480486716242855385927168612482823018056444",  "-.343221776474385556408575072158654617856179903891469234802620083946862969481683",
".7017759049083568733393783941737551070001277872242084264226660399535668",  ".67671792587000736911512136461480736187534959236235094491346138762027883728117e-1",
".954541838909223155766903431818344223068485863314757515131010182127016919189",  "-.9380710784087658097762484082503052387247116472824721232654796919389831912564092e-2",
"-.2409502904092167320160192073420755666334019862970207994253789149806768244",  ".1020605214422607105990853595870990437075190278842101127889141494233707336504202e-2",
"-.6232062231381118266105067050933422974574396911536510116340034383891476033e-1",  "-.91909076955100725951829325919040764273633876481717229682054008156034610124162e-4",
".18087107759698002375525861345661559475998663674749533716223450579150267597e-1",  ".711304824447557147984229145463752577787205658556106886859200235930100976243164e-5",
".1514246835261438869371790116431246212181445075184895257538645572923143880739e-2",  "-.4835756324293811270891645323909681592398463047460385756444487987436062343152527e-6",
"-.5527266481268379550795880932722659280513606228630991949181129939337741603332e-3",  ".2942849149002218729559027188018267600596306488102891779369852119103985985473077e-7",
"-.1729795137474782994916437472997777499321218811128441632457313488037389851e-4",  "-.16195717189082172710338493801671749968585560670076482258914658535896011800166e-8",
".8664871579288053708088731304269250860736151016857447044949329849065979307784e-5",  ".815271499725802393313134383435349666691196407986851770876497539652730335254194e-10",
".13251907213579086685517359295876695247681944359223427633885372043950037769e-6",  "-.37649926513698864265727919510203478129088253099488728165705687710598886733832e-11",
"-.77739008042322191587966949240452675846647701887373000225674589557665010279e-7",  ".16082767699396074704859567485756532915699290217920632552434643225270144985355e-12",
"-.1175325453626862555896169648767889610596426690551070114754609433678796738758e-8",  "-.63224649724360636263001751179460406076461009919686866478994059609165266372285e-14",
".40983846557777730720476702757557909723500712495610110123810716638028072081357e-9",  ".2307932025834496652393824480470290793503479203615299457031527018055603314263154e-15",
".95359115495139813317206453863460782609915466130625139164905428858849777642073e-11",  "-.766132061922048075154067577490846831195763043168119226676737658813113190298839e-17",
"-.11675614242422678319924016907330200404207647028384433579332661902591295442304e-11",  ".235759339807971406363232008821795579386192368375261248989894105747162057794601e-18",
"-.430179104919688656807402552003522912516688433917512060365836447854623095023265e-13",  "-.629027435122489568065817767239038361196670010087520145688458828200929573598212e-20",
".1173874552360932251880167783684049464618730808730992792066448227106179144787787e-14",  ".1556899976976592450372175663389795338368313919728447340288828357103101505131533e-21",
".7678952026480271884985585825300716131771901425257075837076144450714999876366030e-16",  "-.27640575338299618944202366042430506203464015011643500029265911064746547972631781e-23",
".98831457828374892120914442099977214785789921805743205529928105468227697570549235e-18",  ".48742943310318258640311674207543987615891427132184343377784473435184871994822018e-25";


  MCH_1div3[7].resize(Cheb_Cof_N, 2);// [-400; -350]
  MCH_1div3[7]<<".6545290585527187143078214275017971096470789955249423370608101348738086",  ".97051538429717061192632969199408476826014564458300843245368900442488718319685",
"2.7704270072867604162468704077656939052514708598125571425531598529269170",  "-.31411618136968349923522398075729507114147158588237790147310329233731826046930",
"-1.4381829305824509395619437155894102677796553889950635673166418506350928",  ".5836326412036339524190504123329798781431089674102114532994195483438341313302e-1",
"-.50362077585233822373673220952241624437170349799473528462704475964543504815",  "-.769212569243189820958343341009364548237743586647943131199630425140069196670304e-2",
".25747908790129724878883847921994484392434420560533255455017634354569829403",  ".80357689254592170095510609693904648737982743452772913271958178765120137646510e-3",
".15344232259770012727453424979659616089582018007423485240093119205850768586e-1",  "-.69801815420357163432228865274583374249190085081361372371475734929371276412894e-4",
"-.13611067253183938007864455491492219149396434418663283422821397077271703378e-1",  ".52373842655142261580362215355987141869906533983236916965160345119056668552740e-5",
".158369932463116178913492854702970186586953454761805168106810409017886472871e-3",  "-.3460042301172794074362289772228132260331004407193555619222428694911485476666430e-6",
".3141126887808936730255484783277421231996320320804664435926100663614671042542e-3",  ".2052387515229341200895684488699740370935765787574401378391920666799576638526471e-7",
"-.10850150348792096366522282168949788027789688920490075967635277062884433687e-4",  "-.11018052689200858648614006796116143970113431062384690807491632356062333856895e-8",
"-.38405711797614712180516458339467957969366809398174415060910584290019802496201e-5",  ".542062121574080263821528421960286899553705004167447665818284022536829398032134e-10",
".149574880173896180947921109060958764253387785772907023434380269377674135173e-6",  "-.24450865988544869693819627504131795241635059176809362675007675614899379619694e-11",
".277693588375349556881002083515718110314977043753386200785615295234155261634e-7",  ".10215760707817897515991608408843430255656666784292810944100946423148844166701e-12",
"-.9061131407216594962811367350941702273112783167737713133533792680312986997417e-9",  "-.39194477726590049917887563160648213682660057444938209586494083722735473900734e-14",
"-.124214608769007853119927963856213632342128812362556738553996588591117681481234e-9",  ".139849373631193578453343324053472947624823167740122541382496306535101948512867e-15",
".24190197027190311845442498784605958825560681444622939688309683026646668140728e-11",  "-.451715421144588034777555043657611645457898532743410766428976739838127318489153e-17",
".329764136971026688989872544553784940628902722054603114299083416124706726155600e-12",  ".135666793903285810095537033955763840543153806244801321407220880890630181184683e-18",
"-.126762385539320964862744861411302051500614914481569788296808101750763231318292e-14",  "-.350194749479643320224858648958230415739016987129033710908133594546090277918145e-20",
"-.4220348367017998328666957357765759551170527716702366090003422023601901434916554e-15",  ".8448088088342268206989124615496562174342685190464583189060031399471724815256408e-22",
"-.3952812401798806473574817131531491926740697667821840483714520930058427829218424e-17",  "-.14368846686248938845733894677079241440078275458796829678667389047177734370539059e-23",
".90033724201229265002227495351043883172722156022780004792308741028170581188711078e-19",  ".24720420194853137933913517700405320254693175276012904306774606567143532732561653e-25";


  MCH_1div3[8].resize(Cheb_Cof_N, 2);// [-450; -400]
  MCH_1div3[8]<<"-1.9333064310519979882747884309145736038851361945991552808525299001149322",  ".9736541943375477368482739611272408108856741573633813345370844765327547048182",
"-1.1973814771596029938213877582236529083294458008376087833994416220487393",  "-.29817438094786825047586208362226996495093655693355757289087197922984640108896",
"1.848337253542206885245379570723674317380203996283891377746620352219056",  ".5221115306071971646994731683133696369796291117150012006298746153890962140025e-1",
"-.37753672202930144556949261483513408204130339746354524628974993582712280052e-1",  "-.648394363046211493519243700308987084140633242794287937804771894787375178437520e-2",
"-.215366896918764231829894710715420022509251000862734769230123005395506248",  ".63753843684971658312155028598313229548851825295091091791879531944143182043093e-3",
".217104960738341676632010744277568362944195107275005337826841383118444865e-1",  "-.52110287027524784861193871137493088419335493788036996693852741876780478490259e-4",
".7954391137706146237813300806294748970485128280327609376741207682270511284e-2",  ".36759530585382650185067942271637242262699306292273934695420463929027659268515e-5",
"-.1037024857931008939340864527171968658683509191535933109252144135732114227217e-2",  "-.2282369834513028710599349354451405133348752619426604291315160161139753002495315e-6",
"-.1289771275169743316257759250395930101889183079980853386554898095290142746690e-3",  ".1271100328200798936503970562056517580301052551971887626956223350069321465828104e-7",
".1964788073341500825990785377526748380499381687422127029522080756303032573e-4",  "-.6404073527496397052472507865794725960107513945858343312252881845901341055367e-9",
".1099747863271908422936894606476388873414903732431078083614730096138762633696e-5",  ".295293925254375685332451371263847320927932259355308172676600647456226800565296e-10",
"-.19039707546138602307780021533850815782155862920898132501737202126521657483e-6",  "-.12478748287003617744521236317175836368407748775966957584845625627203878946225e-11",
"-.57605709198886113845360101559055489482518339980692054980036818823442328218e-8",  ".4874833386831239375374729379307214619091673140417044542830067032740359402219e-13",
".1045915255122569591659619347756468111380421331833405776941863911920985288657e-8",  "-.17485041869435348434564758841592752913409425070279198547871179329485941038373e-14",
".2198946893968412733733772937395312939100923938937136315105574367226922912753e-10",  ".581280110018641378214327658962638331772865100543858965724182144869307725444730e-16",
"-.33435856189704620534718682154487912290989450544707122423499454212495296901543e-11",  "-.175112919271070906502378829840839458499625001761293935319271982328831323105054e-17",
"-.6587299680464360791672992821651827887004867832866980966637908752450213390216e-13",  ".48712125235532765498938211257660895980261023315552154529003168033654261529388e-19",
".58664981070124161878623208740667735332971116752256620528426164674263760878226e-14",  "-.1170576245479890814217977388623256256557451843734833152617049136470537294871724e-20",
".1310988125767156877345803646514098047976422771097127253874588065273490483906127e-15",  ".2582906738571053983066243564124011176861208107420429343505158999446829805006881e-22",
"-.4390405817915905825728142052717366891985160835555920881499781729914023622470996e-17",  "-.4098562125423943414544606556270743625583068711295938153742772639495735293804979e-24",
"-.11698285625685075997182021333564479951371483710364935442764016177706973146600730e-18",  ".62119189513234599653207670159557756320788255846605557721628182828571863506128444e-26";


  MCH_1div3[9].resize(Cheb_Cof_N, 2);// [-500; -450]
 MCH_1div3[9]<< "2.454692142734856309620391768430839967635854742566831005324877856740008",  ".9798952095605783807812196839525738405604232369992402801046580201130776477736",
"-1.019279626562771358169610711260589598605258153310721921705934688481341",  "-.2564738496234741611511912302353322268580944207880651117679180503639053997565",
"-1.419629026471020151661840802618169433503039981148764142716577979027959",  ".3991271486305546834334481008845889896211670805234254285150678046162720767699e-1",
".4418734305669512577852384494405605628747360613117180756027576062391231639",  "-.443379719768067381143139249420639642827496113660283183065533413938950694722750e-2",
".103872010051265266029772034109923631745805547099516090802155487104734631",  ".39424762430991105560921129833226362262921090964903597318330989277773895937942e-3",
"-.345523661055511927887228138228240086270036640379397337154008197800060145e-1",  "-.29247725666545319772605371561786756024450862125036050914092524264677182685878e-4",
"-.2141655292801945331795552062082658178644492674012109336198774855186867490e-2",  ".18836789795256007380288782998911080198037615686223175093073925119647879011173e-5",
".99770965799883884887038179394609518057565945422054836089307554810562745511e-3",  "-.1069879244798482427082907291084003532532572130191533833537724838336409864520107e-6",
".105401227467303394595648334071995189981554559216002961885458430063580661216e-4",  ".547269880945519870500108571676434113304337806550190821067965501583771323573801e-8",
"-.1401985228582434583386876645120183420413973777493268067484508510007887415e-4",  "-.2534239881132446793022659967867375986475793572787195959403348662008190772245e-9",
".101817211224976147159419664382424278842043770625833737456431536561315033562e-6",  ".107766291711765735130107767413017534720573836248892636409313513258187061394245e-10",
".11069629953423173252556498041294819214884232780799846945707803174502066684e-6",  "-.4197411051535820417769400232223373499454872949932499053582332159842511317556e-12",
"-.12303016099000649151552456614718439165754993708259032774932618316393145669e-8",  ".15167657292774116400542320970103209735671748386175031821802521716519386911264e-13",
"-.529534613568892947207838631038234895152509417568961319101267429610086750291e-9",  "-.50204029410515338393601284581426072840489180610203149193954796856601690996481e-15",
".3936640104485789422082019367093212842519637704896900040532249141615535602977e-11",  ".154824191843808789855380714505877608219605064438569387561560372091153786501763e-16",
".15592155067270021463751859214093152931437156621551684619249540768269441965205e-11",  "-.430180048428599271692912988725922429222089017926892884529727007908005173519077e-18",
".68302242153489543351417893985574165203620678403636595283349811173451179478e-15",  ".111477541932917595812391053039277877996010752146324101707058748467637346450032e-19",
"-.26441161709090227760723674722520263118737601680445858139184247394052562378288e-14",  "-.246201170357590710366043969894663008998495317077932171013646466144489918117802e-21",
"-.262994641770225356915830762312386818699402476533993627669107485541575315698966e-16",  ".5111578830807597022145987448626856636950814017937048882983578389391774016701518e-23",
".1992995584368288366393616224797111274174814670125625805957437735624424161559654e-17",  "-.7373434650527015998491362996649767087597078510929256699272337862756566264198978e-25",
".38293474981998224914295089939115982245800976447823900134094026072778275492538016e-19",  ".10865690539450133364057234948097376730306884797986299639333061287913565488727699e-26";


 MCH_1div3[10].resize(Cheb_Cof_N, 2);// [-550; -500]
 MCH_1div3[10]<<"-1.28205901778203694342159597315987076352276558976699328148689328175273",  ".984547820373397390361072183660706863824908422931929728002487289596743117781",
"2.640039793902011491722366008802648940078370511745403198425012162746203",  "-.2175282830351163640488936888745696319781705567108507291128919129303411300697",
".221130229969509253403137468557525435195056124097343354292137626869711",  ".307120005109184261188981749351637015959016413485690090179097479930307662422e-1",
"-.5462577630728396057229933354011585094429249165961644773887596671059644806",  "-.30973632466452258916267519754312882090479540327575773609446912602824488933364e-2",
".26208749464009317391423064786956962877379557068107959255632732288151450e-1",  ".2555848580330236692767263211759550677298141846914709574165635032440296109180e-3",
".279989865533021975001421989172781758677742262672499883373842536995915285e-1",  "-.1762235257104403995069402606278490059905909904993832539405563152489136924716e-4",
"-.213201769121597151602903594313569551333468962526253325418427283433609434e-2",  ".10697049411502106360833522599164852732047863930851894295959127040307042893845e-5",
"-.58920474880977066961364194575511469035490310717920756317262727201600176802e-3",  "-.572997475253111144697943004679613377395777612864106880654437588810463533563076e-7",
".48988921717351957186465636166417473775807724915717080496770271925245552676e-4",  ".279893379195529022468375290740018502596575422655518774122655842155042677975208e-8",
".645122689320409922259467947662971473106903212214341982881122930488042774e-5",  "-.1236509770631725568456319827494434505696206636352630247213770609526710102975e-9",
"-.516518385359081701153054711519822937483327174658215996772615282422425078097e-6",  ".50871439430889624680117450754378907169001573376236713266356874298784269165432e-11",
"-.4249282704352303262200900836310710359784154513294310424159211526381473352e-7",  "-.1909295008596257743124249044777677083003050591937334234282483315738248977198e-12",
".28788899151830292755229948929933210597115898976720109479889078356680150870e-8",  ".677693844069030155258174653555537687004815504379988678569186681292535771545e-14",
".183158034788958155043016018025576831086512908055757149489270723830412121048e-9",  "-.21800392142569107652562373086345548192420739784748758926446107716371859765550e-15",
"-.8672341212188576507952057539326051675237681125016653906907559977432614278008e-11",  ".67415918170895305155115499734618930098084065366642247679868345771864247513822e-17",
"-.52067460685023481034578014300949486780611481234470960324693302982689646873173e-12",  "-.183002367698210179703147140519094383599595053202985637885931984953764658369742e-18",
".12412262796480799993952374390414552409545304935458472534405220773470585807931e-13",  ".49182102140271355225903143312700174763176387572073521911354537012405269363806e-20",
".88780502670512850234208447825556103209003910076545869099422972670677599680985e-15",  "-.105646600095408098851813467883148497411244194840626851698531570904634498261934e-21",
"-.260890622685312676609290718516103610576709369387082318674088841090069144396368e-17",  ".2430770621684888588851322433138063445842252775441269869931379110887310925418023e-23",
"-.6801055531269607928155277505568585281781219909535414695463862036256651748783615e-18",  "-.3297972725229517524467690298843404058682928519167334734857386984319785938675454e-25",
"-.87831608090279978498890371817130558623369775774317551513634356928007582580807290e-20",  ".63230273982511169128253226154587587413695368062879180512050819692269033480701713e-27";


 MCH_1div3[11].resize(Cheb_Cof_N, 2);// [-600; -550]
 MCH_1div3[11]<<"-1.149806110572458703570811649496313711423254772977141563878059232674412",  ".984231640511804823134791939574450458688653366637248903605839489692608168094",
"-2.0824960234080566886050637310262411185660860377127210513177947418995376",  "-.221151493582644627364473570365360397970824187796141333441943521086727403038",
"1.097297353693172385312049619880532512664321937026551119856995550552508",  ".313390320940717462227839791981705767568537371697235690992763196407052659660e-1",
".22839460590868902792300966053886039444220186267953415768714295365670616516",  "-.31793395084827355087604821365766262799819142093028037269388477501334511843475e-2",
"-.115326249349322715323388566240251177960692468343938761776175233788830611",  ".2626637987139007091927193122353494599261974200323724596827260153581651205504e-3",
"-.37623007891045028407156704625420807790947719265192704996217716052784356e-2",  "-.1814987158969098594854386301169029287119550746981200377616162179994507178024e-4",
".3814434901335845542125567808630363229000412314071121602352057185931271727e-2",  ".1099946168548103825493960493521635089724816633475155660125976468575426496676e-5",
"-.5941243358547024748212456871217718847488647721049178983006885264790591837e-4",  "-.58866631304462831986380692843711620376132638021887136127695694847455383535394e-7",
"-.562820664862322223802608088518821742926712955587439889901193677782405940801e-4",  ".28606203884393579014285202266313445468894466319374559450513326842199901729679e-8",
".183695747047942184706451802892045151560959371104220206410518652416725458e-5",  "-.125849997538162875609795359931074035818233064184493846072360867102348906846e-9",
".441493156348398510696102421772554882834135484861310521579638844744442308795e-6",  ".51262083873231204400595849417523642794378485231962977642606990866494615460546e-11",
"-.16572863692373609581045260358618311832928303151309982826599361708062840040e-7",  "-.1909022773583461016588787059468550963004681761324271976692176571925924686637e-12",
"-.20373207051001887591827823611847010196526734617853090251480404258036592204e-8",  ".666291456727427686286163026512880596127528095181266894818555192991351345534e-14",
".721142700441974845137821642216418352523993659674807357072431475774594228942e-10",  "-.21189577448184050255911684544043485228545569854426453376688407017773930051804e-15",
".5821941070669162126356394434197687932809227525315159220783991502089625470010e-11",  ".63756695380734531072776470834006260695120642869199930709471974362918695479471e-17",
"-.16409469428368117867099630474109013083929629387668756046263808129505859561466e-12",  "-.170681452972168381393352641793085133736033283635188612992364683647867498257275e-18",
"-.10243859626854582802001441427465866308300088462209612754302508616286375762753e-13",  ".43826537527781459043324120503575066160151009869190114098690233887697873534681e-20",
".18357702493564279213657953535126746061181999756224050447302372248930237099695e-15",  "-.93075222053269785759724345244680746274434059380714913918703198946646163826307e-22",
".1005457030029538698417812104975567512342015828502164314668029904429273133969633e-16",  ".1973760076009326793940742587300457990726914048411597348850477018912079466902976e-23",
"-.744901148371692911249125337621250405755659010565075738183830844581233823377675e-19",  "-.2695485517629939373700493309688633574261372971832885247304105323848083076636958e-25",
"-.39829117085809957167851621469924835257748291461467446869202567195908725192975290e-20",  ".43862487842719714071154184725835263951398844073328298272114738961653863483349318e-27";
**/
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 
  // 0F1(2/3, z);
  //===================================================================================================================================================== 
  //===================================================================================================================================================== 

/**MCH_2div3[0].resize(Cheb_Cof_N, 2);// [-50; 0]     
MCH_2div3[0]<<"-.52946990163005309344103835186021354534273679273947373533311470431911010058546",  ".95089363090690083247171969571846659154232497627061763668114617738350944734988067",
"-.6439036649736542067628661008830877201874406195235162560464017356317577748944",  "-.42826672393890121746533268973398066118945560087643831029759974442805263989011171",
"6.667671806816613037383220601360398526245424528315124595457838075743413070775779",  ".96950379200058638160408999085371180164783474979306914706409585510150451761226425e-1",
"4.2598846953529923972697394702795762258237786515209758054084193444937866950198",  "-.14689922853184197415392687117205386423066623848354794549628900355984123898766211e-1",
"-10.962943099748796350381380109047375819265991827209344974832611462045049702453801",  ".16734645150851683664718531893535535751606382432522629466174484079239039697149659e-2",
"-10.725342929056431095570395504354915178219215413337030584322543476602916933699953",  "-.15260999067016962903426789702682498202202047506043994427510233390706113706164200e-3",
"1.0757892016876913814373222629657021435556640383747188751206859662232688671536",  ".11579961855556431064736867610610373877063125287378686854635516894934689985249701e-4",
"5.8877795459356455269516918255111494159911075255091918168460625359157348072543769",  "-.75005296219483704911348553277699770090262167084799520466129646882171042638516075e-6",
"3.832768195793321405050374644204426741215359366252487014631195661100173386055998",  ".42200876603261134559354437393579422314381407111137464418174159500279136940082524e-7",
"1.3649638724978161391720217510800404226764941521806070325597822049035458019818652",  "-.20871841473650032577769762095584031104134556603320774068077150312806803104926889e-8",
".31944133132579230319754939964582985991642413011898239745479552257682691827096758",  ".91440558504591661333942810467047950995899557716614996548218277435800679519967453e-10",
".52778495901563585887871906308384353609550432346105393825684581877292971892123632e-1",  "-.35632008309802945906491649419021080228167234234568740366798276588991665303198648e-11",
".63841133620670291526251948618594284808072808335031347209863676536521246893380490e-2",  ".12359427109458621775454136454686934517293937019944844633204961784004690225740908e-12",
".57624266162118115033918889861848989205441671834906603849202861817484753627819106e-3",  "-.38055397760789200144710952622676877975006296801951401314629693524687489687095660e-14",
".39112569284731934512321463313291062163909566319782590183413890961463437933130237e-4",  ".10329773394784727204437088838852473077037844037816013357293620078348066289860008e-15",
".19921980611939281703570543057621871401441966212243897089323243895233644195672177e-5",  "-.24414314060539097172694446780239740400202895602115433732299624033064402437339235e-17",
".75174013347835237196435622761289267004839154263575605048397132096273333034447257e-7",  ".49231473008247392956891746780380087496650060423278418987187788785417468379079721e-19",
".20439171924298288459954872401268889979071634773580096813319787628853174984302331e-8",  "-.81930468941104068850451785483072534571660113465004806913328164723559256908256469e-21",
".37971272296706133755959597121668882013096598849149775119856540241392915501665888e-10",  ".10628144345535803044266086837360131702500398108766925782013520732685161925121302e-22",
".43284831264962878651395065239426019416511447705848313445985903940764524550741500e-12",  "-.96240065008996928601557599621635257074601816321421234176404589534114770432286619e-25",
".22920544622822601341098258842667330534511290350476309178201415469973307652375888e-14",  ".45969966460139007745936013885193479904985811168050507483823378324373009154099270e-27";
     
MCH_2div3[1].resize(Cheb_Cof_N, 2);// [-100; -50]     
MCH_2div3[1]<<"-.11220319050746701169248401042221001970947146151651046875484427569707461625",  ".94736360108525673442630131500198477574861466043821587413583959668403357025131546",
"-1.37247682754828274528061308772956522619753277375731783834277162302204135017",  "-.44192169989034751661360192618533215805421469453557515779999399576176916135302204",
".961547445038746368475565877960351308812230204955572678687886380558005931792",  ".10380766441600117542470210305550961199498562558165849166904665746867718005984646",
"1.7377494372351302197942530215115097316558968616931184090390874866574669272078",  "-.16353753085783006222899459056411820665658365985736415404634268053928000183677522e-1",
"-.6538170332138823419625342550323766444303767787222253692373083029967990575675786",  ".19412250617553325540502032247561904802456709417761392258460903580915369053468214e-2",
"-.6863159723552202057425216576286790062672896232208651206001913490042293145538270",  "-.18489739899922856296382401008319472814914632219885203924368251259261248246186171e-3",
".895287586045440542136373022981932758398008280612897981315404974597375845437e-1",  ".14691452054510333167774418674963839777594432145343594938556421587897300533119418e-4",
".12346727737970737637456548462367103434150712220099631539489222782965388781817035",  "-.99928084025752904367146132211281797822736525060419119142005434855558616502919047e-6",
".8575356171540882007464806255057698941963148130929101178887121826236259190054e-2",  ".59225046410511523170837868123518559840095637909735993172340602705338901457754194e-7",
"-.907596509212895374900805615324019830929094072257754800257359909368487035147506e-2",  "-.30961586627974071883564550343821287398458077804950756163939555782604664202658257e-8",
"-.22497057913312475645537938857943915148045862886540505327004717377109378053348566e-2",  ".14392344059921233056499211718300954943725489469082104976706584662777129370544647e-9",
"-.32190849400474206504739601809792873660578821776906898736841047556418195548e-6",  "-.59758911974601636702267675967823849410560943244406148917457391957159387253773552e-11",
".85831140705668594767586985416507284624169042927163623554852822359677059982054e-4",  ".22191668210502488829212197417931210919864624331666474056837430384077558942755110e-12",
".18654397486252527448236695590498367370408720022079509181607745220082420932079904e-4",  "-.73545595831176105950412209475569164650666856788950349666114542049159985155758281e-14",
".21879542658678648699253769495981373615283880906747534121120187793338638329727128e-5",  ".21617724789360595297684612362175380454026332660744988176461001818635093243663763e-15",
".16617009564877229107105042750276319620178879407528918098446032708817642926372277e-6",  "-.55711801173274205803755260018315755746229445165880713170588572717454744821065609e-17",
".8580589828798439722123285369487551776073223284515864697211135878833969878953551e-8",  ".12347863722325654577338512491280928025296016753613705927846585708243334619514877e-18",
".30164700541897468302571366381581777748859395205403994320479502784275738582549507e-9",  "-.22796708669954782111840153435692335119430587817116944099306402741659689870314216e-20",
".69582524135244964167362480355968411291022765140069052612272868870199799565967256e-11",  ".33167390128104858947361949987437697484754645809507161114934262069094060426484612e-22",
".95543681874107765734707157924371729625482438067642092122969044199621552042258850e-13",  "-.34128142530922449148628511983338753286189255809919635783806420155639249517425491e-24",
".59518278932292333900494099647340548690852932192059080594526910567773783297888783e-15",  ".18820735394619449454891354601551308615587007884379981084644345831806599184274864e-26";
     
MCH_2div3[2].resize(Cheb_Cof_N, 2);// [-150; -100]     
MCH_2div3[2]<<"-.4794035462159639352968259116125382094217352893635810236498080598493697413",  ".94356924828621622598501858364529013796472546667388583428769000471925005049157005",
".1110860249291494580814799092100336420247656900414092244279601156012100308",  "-.45594992333070556235664534864963551002855524399654471707616053318684982050295575",
"1.1922553646985640096834900206102145367555737152905061690985267067925362946",  ".11115957757159871973506318181520602003566162548312838550786693716773562019398800",
"-.3474577929939249579272318004272643047534754263239151499580435255559261264961",  "-.18213123105415217136568079848347554449820452550651757762860822485022635434377776e-1",
"-.43603649005130560986116420438852095199681945667895356095309455877488130363899",  ".22536100259301002261720940999358483926763293326797735506553354691071327784448151e-2",
".10257451830889539502174003969675136599117408097544811363806222610803879975065",  "-.22431056953158288697505251716798642910488665717538291311545096706819914833909199e-3",
".625014486696914713211605224224586655563665994670146621137769535890432483487e-1",  ".18676142274303468384111967207023684498317474103474713536908661036500442418153942e-4",
"-.9181067701592274607925695849082926410267318691009231937868816155164562702714326e-2",  "-.13351329440237537011014750405787883766004086493448413502534197222273416059040360e-5",
"-.47801682682177756827928402443947974507588105386053167578230394877099060038436e-2",  ".83447077396072046902053839651096588896141746540501104264386659798861161612329408e-7",
".2185460138542353731899956455387607710248594366039364590448105782817260418648e-3",  "-.46176344230849610970358993254339714842820032315929919442194454222545076208919223e-8",
".20252615941031392786192954531479908652706220737891237467456853136307109401092915e-3",  ".22815929310417761078157861901273865312208666894245918523647353951487395407095994e-9",
".928903747915552065891860718215000227887502613885686324136775768817625188715e-5",  "-.10117548310851986441532618644058969110027458886896466623946521571289460721437519e-10",
"-.3904424064785804502430607854947873119201891186389168760653275124225563455163e-5",  ".40343095054742231397897112838783775980685479741611966501607448738848701778582419e-12",
"-.537194714205017660146621944153824576544466086142616350287373671839583933149005e-6",  "-.14445490533340456501580917419756650176774906212464103240833398271863370166548501e-13",
".1205284551407516265493023123888155030178231412420648603843735988883799063117886e-8",  ".46206832619106120944142550868091235584966321484212907444271959263945088434954379e-15",
".60430637831956438013107667403656293201853201319039344643248076835777861633029351e-8",  "-.13068929061420966097705439047931143776704768738263118258841630426433406562389481e-16",
".656574984185920913153518630414567056077020047799808920678817836991035437326560e-9",  ".32112558530702421537888084457970274642965334755571354858473177829609497444027071e-18",
".3614022353888185430750356916266074703346747773860813715176506576237230809133414e-10",  "-.66542867125807388629987978702652068707630952867751558988768668422507565283749559e-20",
".11511271726655334717905489971565410627390725682574905545053872402611152565433083e-11",  ".11035339797521940296656621931486899445395964162184138187491871960313334644344444e-21",
".20331434290573704186324482502555095889842567549526795124477553910941644465742726e-13",  "-.13204204952613935895902852180576609846426672094882430292841323567548127177440753e-23",
".15560717403799937610296911192420255182401719694286773816745931362613644991753591e-15",  ".87026024480936332544219641666824060061119137609558882279828808871116012334361398e-26";
     
MCH_2div3[3].resize(Cheb_Cof_N, 2);// [-200; -150]     
MCH_2div3[3]<<".2280200877058301834330999827894721429593508349360353484103371177913973190",  ".93952583472872708054516284986628383499519495147513014233792957575738519321134513",
".6582356326852413449612585596298957339830602471209059743027423342972335959",  "-.47023290040096559604531107127066245185440147236324979229289833489623115028762706",
"-.728091687926031967227552094300023434014304248586933929009940902552133165",  ".11897237071907634109686911717292175354505120119392769534338658716661948457596575",
"-.21086558192842637508524583467319508987615329810158427833703606703391122733134",  "-.20272454055127698115772624456941694018268924157428589155860220518270436257999383e-1",
".245531132952595860366542065753959728871535090013171301915041676439749588947",  ".26147503729419905120464861514018093017719195319405920218775404997636296500692114e-2",
".102917636766155473867827203952635390763291108106159320444943485806407082429e-1",  "-.27198071573670197483230427031661020910533040185920097412394835019495576237142790e-3",
"-.268198100256459095455068067718263032818828139253336054741954462934283240632e-1",  ".23732082356462690078418183855614334207211379439575622986221582971597010959234970e-4",
".19833571434400457131765153165769426214952018114182130383830751918632717434045e-3",  "-.17835777822863713047321363688142089036035749812817824707860432021382149363136307e-5",
".13826182290137534080363119414803170904269068904035970650211453035207785751517e-2",  ".11760208648952206326163854998539541508038102185965238915842101404333606335830224e-6",
"-.5611856627227507924565560655960193399191013213559472680191879817704820691e-5",  "-.68923301789253412511961232657861356826351120757710019272179315107246892024735155e-8",
"-.396046271435866572219600136594607913460613483541963567780118476035475035285799e-4",  ".36229353744941546592417267725249833350824741586993586663501848519335676601868998e-9",
"-.84190538906124833337502902830140766230879572003462752826076201999959254112e-6",  "-.17178435384403134095217479986044779966398426871358240684150195586052944308373358e-10",
".642655775822915389729976847043882246248018654877615318683690748236792013338e-6",  ".73674353466811125318316068828229819616921775067562085309946850914744913529922578e-12",
".36729380320643299517421918960275372120627265939151830257240418491394426412188e-7",  "-.28569609685203501120848047211294254708649768642690096038111268044584650468444778e-13",
"-.492510036872860102848369396828230045089722404008304999422847935132059640629871e-8",  ".99781368035892077866764022393954627711435931859482541792375980949430838007955741e-15",
"-.544189292074584080317575080087954525112783005637244525659043808351891100737683e-9",  "-.31120499283179091638731342736424364907517914016171660273392782984241041715280966e-16",
".3304393289883925592195652509239974132189287720736461543989263489154596805965e-12",  ".85362648550269616884286505453441172131172675285462151399390653132303326485104988e-18",
".2555065045851583777068847409492770335568855933268002220064593266273407588014853e-11",  "-.20058498538187481401518650284022994608284002075000183349361036155078536031846684e-19",
".15696768035301096112777322638028877011419433472555409160716895647345232720046882e-12",  ".38525618831060623083529449635373827409965645229045726815883058572055131074444208e-21",
".40757458752820245372203560557380422918112310709076204089885277710038759737529921e-14",  "-.55030846679595310523111596783923492257347144352016419147828064931679278505176928e-23",
".41057897000230733195813397851499334765255651754481557802688919403173685348225459e-16",  ".45653355693943370772920124812502152228771394787439551525698048172485801785556947e-25";
     
MCH_2div3[4].resize(Cheb_Cof_N, 2);// [-250; -200]     
MCH_2div3[4]<<"-.49603497636688536527575951061245977193116887136736646083577506002013704e-1",  ".93526099655660957397529106447910616429471296393576700646875097396780882016260525",
"-.7285775111710783970199626720316051929468783440640387916188814417931896088",  "-.48460849707052143480145386699373617451337972340542509205568067016388342895953617",
".425044241306469541551477738330992590734436845599907811928488203045797470",  ".12718816441524810403091145480760159402543459980010763519582118911301007536974255",
".22912835975988811327584917597261827162603514368896444024960943904239525539405",  "-.22529827353259547799005182806300471268887696653680777868246826545257913930614230e-1",
"-.14209825003465472832166288496285671038844933242354091814793265995535417721",  ".30279378079238115629396136272750719033466772134823828921269324666341504435719798e-2",
"-.11667162243535033217302932499560992612574934195426272309637316294897623267e-1",  "-.32902863676017521414318705235585813769851991519488777589236033399478264665035696e-3",
".131505168271201508005782084299770127061009027472507589279749646744431852793e-1",  ".30077577419306472094214422651877250498220904597768399937689024990030716919870877e-4",
"-.2023005506838459275512374389126384529085313990094408281773027613393141418132e-3",  "-.23756462299856262771906403963206653806461201586429584291037954811827322125234674e-5",
"-.53086152997288137342997831322568015227084004869563265039429905054511906088435e-3",  ".16520582684396854513855708067390361555035433911396224411620947542477346780554134e-6",
".22155809660072382176857595605690197401651203063832710478409391012249756489e-4",  "-.10252580993914221676501734003134463692435512544165380092319345285510291066727719e-7",
".116653415271620862240544233940355297749086436705675211582772906946848811659810e-4",  ".57327790954422111221980499350237103791862978601629244447022474167100150503203492e-9",
"-.47154598601336308621847280327477123931146378580221762992568101327541479507e-6",  "-.2906728289424583660666124126209988853828856252756992411396049740386338378524998e-10",
"-.157141217565569545409711169790102003089090663402873742483747653517956539329e-6",  ".1341271080580258737406208606827706762383211384402412275814895494936927503405239e-11",
".3677323873671018486474115258880156765477872219515728627005832728348418362628e-8",  "-.5636739222474406839977410307206031377630119094813704335640710237410713847428250e-13",
".133196789830373254859958925729649920401548353375455629194387476656602516034511e-8",  ".21523516104983538049707913899711124613565890426773389770179049634437060560312124e-14",
".1803774368160534513861742814546430913294562933759953364111219224519817127875e-11",  "-.7418557586677738838222400407682261093791646617603693940756960565267523308805213e-16",
"-.64553263213955327086407851786662552866329594797787550270951114557406760738442e-11",  ".2280546741632314123035793987518387700019627713690780499326699754684003507019092e-17",
"-.178879718104537695232966516243514902678154504208249511643508438324634963053729e-12",  "-.6114284070201761550895680809692387698795227044425503784273058996183793344835865e-19",
".1212041728498175667938273988790618442700848844580504221276499739999485839338900e-13",  ".13780856227467305891820353855486820939523939750113765569079275290408287546896125e-20",
".7262177983479121370526737373832037247593960643293972052914336539669115133988427e-15",  "-.23942618854561556142152603884217805374497819517576233860846652824397765926222309e-22",
".10934173726070890802341207110066986663972701275236041083141854924408168145906854e-16",  ".27264767532167058385021105338445426138805407814249582136189781097773070051730063e-24";
     
MCH_2div3[5].resize(Cheb_Cof_N, 2);// [-300; -250]     
MCH_2div3[5]<<".38464716919116923557147841273112371263743194560914494508303989160727069e-1",  ".9541443142937504727965804004718398379969450328727430554237060505224368847185936",
".669754654686700650856017222553479190643424792312197536623793546901558626",  "-.4036872510073365995054867937270155706870284822399900772726143662184425570975606",
"-.30952561557565393166007795640603815821031376747417169288845202766417073",  ".904464538710109585966016204403903458622523719322448439054450163327442853975677e-1",
"-.1868529058145369199193423572842934204703189126769416887087691547955425665112",  "-.14003062274379615403405858279093986338206684765014566385176781969640266307486693e-1",
".8678170094988426829718008257208046548494411514242705894191800844707753662e-1",  ".16750873414710269259917357080693576969334767023623752588620411639668829114629896e-2",
".10543561399984693009863020986878270636052231162290491847120474627507025721e-1",  "-.1642338634809199692224882703767535243638766280317599667591426606849849006222653e-3",
"-.67279095324367302426928201506789635305618450910558299016598928520685882274e-2",  ".1370657153942666748964442619190131794055217536974009850825083741557412464449283e-4",
"-.1184950478672704736174120097986847505486823838632779550870041930117870831547e-3",  "-.99817611835267096649441963560360886956746060422536841999109156401306689162845606e-6",
".22619396274232331338568147801796888891337224370193033554848548290611155391629e-3",  ".64603315521981933360809957938822089963465686472045997338818167911695343292604214e-7",
"-.2842301102917248538751684437413521879945332888028037322690997821118202560e-5",  "-.376269723075902286338013226945538680011771451644434503802075303588762576591033e-8",
"-.40697616680658894405207361014949657825294111307439439964400329618775742615954e-5",  ".19919403633543258957719533959749713216322498733160200833589243435806549159384604e-9",
".6600080391606529812648658281799260121336805466373743315296559750592572088e-7",  "-.96399429731251541251890269159927730409053516489950247831998954200595475100202e-11",
".436594824914741590900603797918639707479568324017275096796157550267307707719e-7",  ".428698968321144611957844058791857123988086628937598274299198185042646501725325e-12",
"-.278554867600783091955705558941512104335034660764442206448649018269291446884e-9",  "-.175140991352315550040687842117540315980905436291824427729803158561300786289968e-13",
"-.286018398863755721493238323619829341282164724022976253545467976240057882121449e-9",  ".6585855380373573989269301784403635681495753576951251756734539501015744375222134e-15",
"-.29748006212534413903331705640829775532992938221456769467437026754856715907567e-11",  "-.2256712795739349675304409051533079455998315210765450784105171328763476566947335e-16",
".104646885185307644353159877138529175507078248472687882795454547164283058075414e-11",  ".705262053739001253146161219366497564673735416763149637314976178446141299247918e-18",
".316149410989628705476173338564716010502667352303869162652748620108530358212334e-13",  "-.1936535047644083935282057038965989242822712491736616044351094355268846957513440e-19",
"-.1437562306251879459107477261991023520242646115733971497686590166370705695555487e-14",  ".4731647939775823421282249957291761022522738338334707819178156781653453701290262e-21",
"-.8480495967986884430213941756946143552068888015531584974872673336293809413075057e-16",  "-.8735846243944889098804690364793713839723646481141130905720030298890035831867281e-23",
"-.11583561377138855356333334259752103737717009238191854690443387553902007367868360e-17",  ".13878491565286315331912374571134808814922633663313909499768410111060418420105849e-24";
     
MCH_2div3[6].resize(Cheb_Cof_N, 2);// [-350; -300]     
MCH_2div3[6]<<"-.15016502768847433754240827485228419932137529120574435194070442791502327",  ".96389670998981057119245138405395380961527004369631939688321615571336898402247",
"-.54133587617562877996377228997650381797841331528769926688362903338995361",  "-.352494950120898933932550219628050313639046216688550383728378012720127280481116",
".33603288676805831604504021859146989110048465756487105997532777162069372",  ".71351059495034053346162324163280205241206017757244142784512766745174263615079e-1",
".106856560024053993624394433991062806701336043911751096428152897813712277007",  "-.10154220074867392310241376616251498324796388940667216681596610687038765720267725e-1",
"-.6684957059406780777314540074558124081364482619776590456511456346887405703e-1",  ".1133734822734979851168469977486166736648171891373075778797879546813204171228844e-2",
"-.3116027157495336588263373456816704613713015131253441363834140321244358582e-2",  "-.104777782684590722972209151603445902955842546177079774541223350532712323252673e-3",
".40045830931767711297411936237242644986564075898658088004242462206644399936e-2",  ".831929462237383191883472084878482705227976998769010408194965177108886527180817e-5",
"-.790314992535213826474189493816047873156337142878132566911084836871790653401e-4",  "-.5801894786635765598799761057572332708143196898228072234250341703867738402915521e-6",
"-.10643323531402719585445713938281379596467250038668744820924774408800681531735e-3",  ".3620068744414212944847062233116175632608436070001575871129715218899516151843955e-7",
".4138301871696132273662291236751546568261068040254753791114200098047483256e-5",  "-.20419934671162806042352171775967991973549441314342989770926721988221160966912e-8",
".15269817848431730281628391085283059783404020781506338301338480791384409533795e-5",  ".1052557226791010403220378284319818025895287359064553484202476932429414776380117e-9",
"-.59087601488285684872864642325598084755579557253760207135521137488814747643e-7",  "-.49746266225606450706085814332729667105846071751243586990606853463404182824536e-11",
"-.131853351215654403077593809376395894723938545274056432834405676003953198234e-7",  ".21710778274969059422207169042875629005034492900518152136369622574275693244732e-12",
".3564556291041937933935156103736586276653287368027265771385881560582579433784e-9",  "-.87146204263561374976567489267890985812025760332300269456680964773171079920866e-14",
".707913258537246384126513206839538123371006762514804445898533051075890986495039e-10",  ".323801312434319950503899641151288403404109925333001995218547111793483890408530e-15",
"-.63603628094139224540932365547984368711798304312423937913806573494102259575070e-12",  "-.1094071323160625710799239189250540331824150821622046954235963567111680930946098e-16",
"-.21991348551102494533952249913522915780891096150341357067931058957508116905880e-12",  ".340449816644788391754989738952092957937050112838042016904595628142507404487219e-18",
"-.230671011281481290255947382336352432951533082245655803295405233735918509636544e-14",  "-.921522100238562879658441230926053756254768908100741383160196293187136088537745e-20",
".2919185373689560124599682768919386223794127929024019627843738488798589677245908e-15",  ".2274985289966405691520739813973466895845024390577221905982270595827955315176711e-21",
".8543370135420772140135740368837679776237165151439362684918925801047775805676410e-17",  "-.4097923883411163626689153804682447304823483070527100198517447453716644513133839e-23",
".50276777621437385801089050517542734453598148479108558850344429150177296942908615e-19",  ".69226267341224267053309195937191312783706053488737571646732502983979481939431754e-25";
     
MCH_2div3[7].resize(Cheb_Cof_N, 2);// [-400; -350]     
MCH_2div3[7]<<".32473138252892629268778993735231522296503171021420715853004716239384984",  ".96756400859217219051001449632961194685648293671112986585637179190476563270547",
".294184556117293987203600939275188488822283339154282408903251701050625190",  "-.33317724102193900141282792351590703496322591307718138517829621358694110157638",
"-.37985640351509102885165365359389611138075594959986673009248905238904659",  ".6417134071868036650013322322381275079144644579863865736653108232922732275745e-1",
"-.8386310575477176086237392534051663495323666036823497622394564381842911877e-2",  "-.870846833066152349152066289473838703353578230769473887557869215119433864748342e-2",
".523682360949136719501512255037268272740552830037410553208225849215009108e-1",  ".92898093954608047887134902649437921287609779646144408519849831209226929703096e-3",
"-.432632625386152210714515562023367461233852821743898484880994561234587119e-2",  "-.82056141863054586922342469268358733082587659927001922602973148203114372596578e-4",
"-.227545086279443160991750576172474302952913289032285511939520322237916731033e-2",  ".62285344332826779994005191936218245302667675433347816006965900909404829110026e-5",
".2630710904337114487040467895684675915692967058682834158576231126918265645814e-3",  "-.4150572962129246240333846754330944989438008404882522240725996536479515354438833e-6",
".43756618433263441276623712473389975311326518945719837857678538199805287956161e-4",  ".2473334924683594754290166596559977846876402922483927027893142782246333555266342e-7",
"-.5819401017517230622747606010671387852576359964879689857162504179402583557e-5",  "-.13310847079528553781716508956208292400541931478711355113169580658277001262067e-8",
"-.4506937840013896818632241961193489094231045640021521883741455870927500696334e-6",  ".654006387638422187562065487719102705101571206556227591505367838589632996576568e-10",
".63950547601029167716230711242310553847868903337975915521984065176015442347e-7",  "-.29421689890983799971248622660122737027075995980396916693473178064504360469921e-11",
".28856892524873047610392121305368135195716398553330755625335917903364707956e-8",  ".12207220453878836657919675163505604944450903873400120480757436945760551451440e-12",
"-.3892117932221497524172112496580980673309527146563786469190692017624161751559e-9",  "-.46500078498139098166891846240688885020280402727166919925317192742867698911468e-14",
"-.12906752003915196904614425185308313283091829716159460888212730228634356198271e-10",  ".163711558343254263514693316564024072559770167342794615207903985795618402558760e-15",
".133614272499075504265868523516080987918246290218533100629489286642210929712953e-11",  "-.523028088245775963506346502524224245989484353877627389977716601397440245094411e-17",
".40623176489930649945658039929399169950921659938808167989035961429898120113160e-13",  ".153570415613022131152911115873615554502714165898851563506124249860665647202227e-18",
"-.238779221067588143855677454093646873322970102528772916138277162784276244594372e-14",  "-.391370697331181548075192456744421760268937911115665352513657800844506190984683e-20",
"-.761312693375474686433259442980509555473419867504071999969406471370098043263494e-16",  ".9060735616179016243325600500184992881386844711749364454675264257251314237627002e-22",
".1636494154912954648083178334343095511591843859127777573757186053547177439956193e-17",  "-.15300644923708478361977258448752182277182331594983158524163826332586884227329847e-23",
".56436196421991921873014384642593618144961209919140694462418604015134216504894144e-19",  ".23899951322642848343137420121236098229401960395066550888647967093116638778657860e-25";
     
MCH_2div3[8].resize(Cheb_Cof_N, 2);// [-450; -400]     
MCH_2div3[8]<<"-.4456529843482229284266704545296850298742822385760418512538282565077424",  ".9737805093165536041799197446786768541751322497063393452754037268544010592412",
".655767186488698311609604697356793960319362830539758575327030087447692e-1",  "-.29869367118375783962272629137287556899508504485023238940621100153471820133490",
".3243231356112082374934610700450643971999723824163021529026411939854348",  ".5196908958546994515346780586652205008497996474710376812093141397500668635634e-1",
"-.72854995766301044892114794672791939426962939513386717893449269915949551707e-1",  "-.640494513563881716342687493472586323313659882135514019412700885906171165332086e-2",
"-.299065732790579064544433548311558605758880122719613326384867392186468340e-1",  ".62358514349564018941324854474680702233794737258286106430235494754067426195236e-3",
".751499264522436601433963451868718889723866598415098322778506380047395157e-2",  "-.50427686530744850767325369828180188050468638617711086301991252896097999670565e-4",
".8340438756159120254047555726715569278225295364078687788928061767586397211e-3",  ".35143064257197515649221690653925140791593104092617950103860791469149551927314e-5",
"-.261539886735274236876360453792971908588246188857386375171530558497054067184e-3",  "-.2154244801352377106162399781931688148141304293577199473396083605964034178438199e-6",
"-.858696417767358726093578415846361212368941608557276690853257419731441279337e-5",  ".1182930343843289642827125077521662769308459055628660802391352642323951297108816e-7",
".4277138789810225151246946847870978637016059985199867529556021723471020812e-5",  "-.5873046124798696193173297080272134935475187949779768305720174292285579969906e-9",
".235278350307142296136482381959699895501230890819245062292951196448372383448e-7",  ".266474455664990713636616748032374513082737164697349423872650043259941381506472e-10",
"-.38503257017137566224224050783621009668791595382681532174893496816987943868e-7",  "-.11076010848089951991602762488420081933814023532104888551274319827999323708736e-11",
".1117885754629094796118956401680229090961911448192053995936811077895518775e-9",  ".4247515703541338042851705413655096683273753697058751923860004154747210246345e-13",
".2068461115300519192746242053741488717370032080042197736564565338612172708306e-9",  "-.149548717917622561528561574307561566408854332950420950779291973947080915062502e-14",
"-.393817052460322194603705160796377377802133383347366134068711389590490055766e-12",  ".486474946884945291616121594969166341614390713997211019459971639189441399044172e-16",
"-.67588779200372086133304798163475478350067000718729903273526940190983013580673e-12",  "-.1435589666742109714970889269069207867494426069058004949470705333889783490261506e-17",
"-.2863646294813007599833948547696421984109642361093836952667186814327207404688e-14",  ".388672538977309685773400120327751901388292143638189370939648137137642039973557e-19",
".126041859860000294949932535829290996209210082187106364350621620406610422338416e-14",  "-.913508888683108591086019328557008607657326160947661963136801365233090476063917e-21",
".159703087800865152976100182217320230137938541768223840952028673971443391725332e-16",  ".1939026356392254869072152325582209866334884874551232615343762571453250768395683e-22",
"-.10374851764112775587909218003325372904085455913940937839030801120233377097733702e-17",  "-.3015929924669163871497091649045383172661234933645249673122319563561738690618626e-24",
"-.22142365862516093025000990456702496258980535120344397806134579184957734256711601e-19",  ".42335070623313972528687409449529820629949554014647391767602024022353083366442897e-26";
     
MCH_2div3[9].resize(Cheb_Cof_N, 2);// [-500; -450]     
MCH_2div3[9]<<".354982204123036537122934812677852053723567313956514116852786433232330",  ".9805265863994530443705619447683620194490195925218859084045526455654836846054",
"-.4034509228788919854276764868639582861612654521703659778448809570992082",  "-.2516894272311098522276549526811260607547022362313096440391168365875354086560",
"-.1448502827525967582473354552281827942470294170256536116580336666718888",  ".3866630615123782713303224290089822593748606333883625077995444522609970773444e-1",
".10900732985169172111864360357003576060448044504323258308886817979033786345",  "-.423885008031102280286423818581938504966608048258369081082268685242811559023650e-2",
".46384163370603736346681297334070290391344261406985410678129823706761525e-2",  ".37257425801730720366171508687207946687631264094631675466915963393373098799022e-3",
"-.690272781235070999123318169450751392368306856540746402052923933074882386e-2",  "-.27314759573968751934646911837532282199550875381637828630255054949794841194123e-4",
".163886353423154086435314347408202493123436668081109484800621060283691244e-3",  ".17402317846839498850745753158022390183406518408389704932977575591253262254302e-5",
".175098782877266436503029118076102819610764792700282698668847228460055773402e-3",  "-.977464674555897213177762231621482483062448118453396728443538261592118543895379e-7",
"-.76235019423062663603307616344100257799591896176178343448004114203118542932e-5",  ".494951932650545542770687548682221059762740190677942247184507699765345098491452e-8",
"-.2266896209810686448035490240768547983801890280607647249116126866877513923e-5",  "-.2267909160594711226561561764312304796545037733998133328855043407268771027700e-9",
".1061056166283598651730447629487379207581279217065629322808840903219035461384e-6",  ".95559631321241920597028479393143878173274485473967205171342535591534170591809e-11",
".17204633823026743280683225451193973543199981280017741986666062973400724328e-7",  "-.3685336455713952264113816361027444666320849169727315407517010717631060040707e-12",
"-.6980449560054822187972743430278642931820222650385386064809139288470689025e-9",  ".13217549252726074111104003769338550424471289833786657721988011411342785747338e-13",
"-.825567856513354971465172002208780448664310527612021120876326397601021776924e-10",  "-.43354243867953863255822780883911380731885711392934591208230959677647998304211e-15",
".2297111520977349449166961531933389780342641231076703726453306583790113627990e-11",  ".133130287117785540785015424886915102945873397546945943241588335430957220248858e-16",
".25299629629591558539190950116471098632125081070914444194467922758419234584862e-12",  "-.366844480602631582385022955515278979206510313265299723825594877089972005101215e-18",
"-.3021343621510786159652390242716697380894204105073398255028332631001442414226e-14",  ".95330809678715668958230896385455316285794321486908425168974646703923488086808e-20",
"-.45684607008935678928618943541094355284842720889060231898793652429802450725842e-15",  "-.208681570398316608628481131470589333667930088570372875708801509933145265601715e-21",
"-.161126534096390509921353500221302682302334243348770349220206754247152817550967e-17",  ".4422855175914672353519978235417122379571141434773597421231945611933580071204731e-23",
".3694918028161041430826360167286631004940075261258942338176964416607328586499529e-18",  "-.6276553655624476791243430913192656022876455444703863897107293018377521831460179e-25",
".60844687086514695036670441586874945375682747487731529161907731449552918391906203e-20",  ".99795555266964714922723084059644761449489663377667461330968968990999233412241183e-27";
     
MCH_2div3[10].resize(Cheb_Cof_N, 2);// [-550; -500]     
MCH_2div3[10]<<"-.4168652041586062919524750321668242850429959281786470814578557619184e-2",  ".983648630089241924119942496798238121735341662680477830716240719989775846021",
".487442769142362059359116319266157738546782461840357684854666080502413",  "-.2244284949355940316968129908980841959705258146415164791627931043980734241501",
"-.100931805492284974799644403592174532480835708677094696610611651807192",  ".324884246336766182854823136831674480878174332243377450711514564574603407481e-1",
"-.8250184982230985094189788342885014816314734915593021883404468841119905568e-1",  "-.33657260674666506972829357517219226719801717470884340887388077365331136841698e-2",
".165649199948592268249852561993468638654486454174019102794690919503682088e-1",  ".2847065378769910880847948409408254064273655122699481510037021911281840835755e-3",
".344241523813255673423696256757254639890664288251530811197251146067538240e-2",  "-.2014808073076201599201465588594316673371546141723222598003289060550811505935e-4",
"-.718481715907895377043451067794628310817750840441719283171259074886585758e-3",  ".1253434403163902851578416028305989176211264812886005102747119648169409743778e-5",
"-.572062293693858905656218079096573545571706070850503064028626487631965278395e-4",  "-.68886442208125033648095134557924243638509805444686562642842091868981032825623e-7",
".131639071301870954129836231933334306106242917102433647104549266425196894428e-4",  ".344629477838938048184404244138898679763045158444145250916997225521594131154354e-8",
".482745152434361235907004525670876119492753315317087040488232182619026086e-6",  "-.156126803163629015412171182628726227928178858750191569619760200088781900026e-9",
"-.1241143651467287811621168442920499262834764086573229740703971678758378979813e-6",  ".65690936493207961170522752948301390072428999495759545255494405788271743667375e-11",
"-.2607874586772353800662731784085828359495011465018450530006298410861609700e-8",  "-.2526257381457451925360144110964429257307468609731223572673538749217456799503e-12",
".66138138658760172213676016982908026263375225436737339877379692555471533570e-9",  ".914517406250778276693687975749372317301679762186235370408423328434948850489e-14",
".1139966128934266324156008291035187289172341859455941647841952302589502343012e-10",  "-.30115174083641242413685320305377781779526635190429966029124800290415316950895e-15",
"-.20362343825235872451660740417414724098320731428830490691881520915970150573989e-11",  ".94498129765512569520444496214069512357885190280668307293333899697463486123260e-17",
"-.39999846138527343273439857500190635801842153458618535606272372844046161664764e-13",  "-.262492064521134185945751878010791723470872122021262297299011385547054665512805e-18",
".33732799403481039851646202765035994008983400663528885055120842969099954885071e-14",  ".70887423024957891247243265942909191323693178308504873887700756962726665063959e-20",
".87093875478153549694763787509122524836959179339262057149135354634612847217312e-16",  "-.156187933781313781531177552481092439489451032300189223100695072858191359854295e-21",
"-.2151456029550503130161181938918532794805560991693126874957936432696568014582580e-17",  ".3538931204438727794543371384626067802060719729158892172342895071588465048768788e-23",
"-.8116248795514915573407539892989991247224387918485071211552737769483062982349150e-19",  "-.4980400130650788152936829693163604700516646086331249909860858882348546256256884e-25",
"-.50234873858481807473282860349634422032701448733636497306592352890607062875559357e-21",  ".89808065524264318845369805942014228891801057600071251581721916683228026912667096e-27";
     
MCH_2div3[11].resize(Cheb_Cof_N, 2);// [-600; -550]     
MCH_2div3[11]<<"-.372306894365819424436548109749783636320345381776621747337292870706750",  ".983707140333888523380568677620301149196943229437923749946463152558286209169",
"-.164604711009804652119215265660918909807171065171093462873617159146758",  "-.227740168436437373186557807915701318709632286797778528179944812256666023690",
".244513365858847600315614204496075261160736181557750786666962577919078",  ".323812363078431109143052494674975731626608150449776924773836058100907018505e-1",
"-.277564268590090857112953637129574813073690738706103562288692648209736048e-2",  "-.32994008272471247948055077776006881419185671341692706699463863554627494792256e-2",
"-.20893206598524122118319644538269115157281829706163970006018748493921516e-1",  ".2717036138881774411063093585530142272505005469874215887914795317869435289085e-3",
".14971940785000065470497067309273227767456029742326274144335515800846260e-2",  "-.1871944181419752306493907154528251673047253928633892204892068221297370128308e-4",
".589576210581981093841187453211713932497720148883187992305915202565827449e-3",  ".1125981724916803476544389600340759670008127276742053668632066553266512596218e-5",
"-.53651548042657318259105281702821209203078141174964437795807387180143021580e-4",  "-.59838033947669755535417436340287905736207766976564077521085732181048678573034e-7",
"-.75750023018315362927475834548491232163885917626258314974503903358095740411e-5",  ".287617232833565686009481065382045823132015173839682326422955708320799829019140e-8",
".759482753215432373720285829784038003950885415312038833721826827050853563e-6",  "-.125280851041910053657327247772439735110993095646233652731134786892003630370e-9",
".528241692150434355174975937581235918419576758970173837320934475640229552466e-7",  ".50306642389774952145949195075157947614195071304433475301428764779436663692725e-11",
"-.5476621413212174323481891808571968920425488729713788725793135368143753129e-8",  "-.1850606288232628197151156449251756275041304190175990291771818396261849220555e-12",
"-.22659555773726547618839340455974526500656870513346185171011327927917191925e-9",  ".634276572763148708382736692889454320978355732455427004895353521792493728205e-14",
".222043154152734746354299562707410945448706314986811963775042700337964405213e-10",  "-.19891543375387207691500891017695366883028053602482787248918258469098937421083e-15",
".6563742253472031464599289503510732621555772812638292519407273100471606862449e-12",  ".58452958627420054844918428339185140335166946096318292763921494790757112794328e-17",
"-.51776249119620049935086834268020647631203710138754908117769389141250917807084e-13",  "-.154229586336128073039547365743047512131962708616882015829369916256138839564281e-18",
"-.13158021188748428723932487150711744674914528942059073604581699877719612301416e-14",  ".38313829132055808309432231227910314929863428861952439647329433368485389637412e-20",
".65426140098870321308557279512097712884789634095609171701756278044175031833061e-16",  "-.80407529092516287055641927886672532788253462673593380013461762161920364700646e-22",
".1648462413884958627399928441842366771169522769347382862135405865150344249472747e-17",  ".1617100360324850255609162724574359432247482781650187631401513857877631975608199e-23",
"-.3490968566968277847283379362196907762915967011120466469300358144557108864253988e-19",  "-.2208908683923628429507155539733296563018152492386734249074746109491796124002904e-25",
"-.93969553834291647000543210898217582929506426535292993216779271618787115347230965e-21",  ".32345723920435525386393337930946724429161826640949477510608201404614781512250853e-27";
**/

return 0;
};
