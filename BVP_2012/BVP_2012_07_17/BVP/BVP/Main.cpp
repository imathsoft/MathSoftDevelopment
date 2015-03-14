#include "BVP_PC_T.h"
#include "XCannon.h"
#include "TroeschProblem.h"
#include <time.h>
#include <omp.h>

using namespace std;
using namespace mpfr;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

int main(){
  mpreal::set_default_prec(128);
  mpreal prec=mpreal::get_default_prec()*0.3;
  cout.precision(prec.toULong());
  mpreal a,b,c;
  string st;

  XCannon<mpreal> xc(&Troesch<mpreal, 5>::Nonlin, &Troesch<mpreal, 5>::dNonlin, "0.01", &Troesch<mpreal, 5>::StepFunc, &Troesch<mpreal, 20>::CheckFunc, 1e-40);

  auto ic = xc.Shoot("0", "1", "0", "0.000000045825");
  double aa = 1.5;
  a = aa;
  a=a/7;
  cout << a.getPrecision() << endl;
  cout << a << endl;
  a = 1.5;
  cout << a.getPrecision() << endl;
  cout << a << endl;

  cout <<  ic;

  //cout << xc.DerivativeBisection("0", "1", "0", "1"); 
  cout << XI<mpreal>(10,-2,-30,14,"0.1",100) <<endl;
  
  /*
  while (getchar() != 32)
  {   
	  cout << "Derivative = ";
	  cin >> a;
      cout << xc.Shoot("0", "1", "0", a);;
  }
  */
  
  //cout << xc.Shoot("0", "1", "0", "0.00000001875");
  xc.SaveToFile("D:\\My Documents\\Maple\\foo.txt");

  GVTypes<mpreal>::FullGradientVector AA;
  AA << 110, 1, 0, 0, 0, 0;

  GVTypes<mpreal>::FullGradientVector BB;
  BB << -12, 0, 1, 0, 0, 0;

  GVTypes<mpreal>::FullGradientVector CC;
  CC << -310, 0, 0, 1, 0, 0;

  GVTypes<mpreal>::FullGradientVector DD;
  DD << 114, 0, 0, 0, 1, 0;

  GVTypes<mpreal>::FullGradientVector HH;
  HH << "0.1", 0, 0, 0, 0, 1;
   
  //cout << "start : " << currentDateTime() << endl;


  /*
  clock_t t;
  t =clock();
  for (int i = 0 ; i< 100; i++)
  {
	  AA[0,0] = mpfr::random(i);
	  BB[0,0] = mpfr::random(i);
	  CC[0,0] = mpfr::random(i);
	  DD[0,0] = mpfr::random(i);
	  cout << 
	  //X<double>(1,2,3,4,5,25) ;
	  X(AA,BB,CC,DD,HH,25) 
		  << endl;
  };
  t = clock() - t;
  printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  //cout << "start : " << currentDateTime() << endl;
  */
 
  //cout << X_Func(AA,BB,CC,DD,HH,1e-40) <<endl;
  //cout <<X3_Func(AA[0],BB[0],CC[0],DD[0],HH[0],1e-40) << endl;
  //cout << X3<mpreal>(10,-2,-30,14,"0.1",7) <<endl;
  //cout << X(AA,BB,CC,DD,HH,100) <<endl;
  //cin >> st;


  //T_BVP_PC_T Troesch("10","0","1",80, &Nonlin, &dNonlin, &StepFunc);
  /**
  a="0.05";
  c="-13.4";
  b=Troesch.ShootFunc(a);
  while ((b.toString()=="nan")||(b<=1)){
    cin>>st; 
    a=st;
    b=Troesch.ShootFunc(a);
    cout<<"b="<<b<<endl;
 }
 cout<<Troesch.BisectionProcedure(a)<<endl;

 Troesch.SaveToFile("proba.txt", 10);
 cout<<mpreal::get_default_prec()<<" "<<endl;
 cout<<a.toLong();
 **/

  VectorXmp C(160);
  C<<"0.009", "0.05", "0.05", "0.05", "0.1", "0.1", "0.1", "0.1","0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1",
	  "0.1", "0.1", "0.2", "0.2", "0.2", "0.2", "0.3", "0.3","0.3", "4";
  //cout<<"C=\n"<<C<<endl;
  //cout<<"Iterations:"<<Troesch.MultipleShootingProcedure(C)<<endl;;
  //Troesch.SaveToFile("proba.txt", 10);
  //cout<<"C=\n"<<C<<endl;

  getchar();

 // Troesch.BasicProblemSolve();
 // Troesch.SaveToFile("proba.txt", 10);
}