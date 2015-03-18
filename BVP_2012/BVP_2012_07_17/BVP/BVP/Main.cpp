#include "BVP_PC_T.h"
#include "XCannon.h"
#include "..\BVP\TroeschHybridCannon.h"
#include "..\BVP\BisectionComponent.h"
#include "..\BVP\Problems\TroeschProblem.h"
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
	mpreal::set_default_prec(64);
	char answer = 'y';
	while (answer == 'y')
	{
		cout << "lambda =";
		mpreal l;
		cin >> l;
		mpreal h;
		cout << "h = ";
		cin >> h;
		h = max(0.0001, min(0.1, h));
		cout <<"h= " << h<< endl;

		TroeschProblem<mpreal> tpf(l);
		TroeschHybridCannon<mpreal> thc(tpf.GetNonLin(), 
			tpf.GetDerivNonLin(), h, 
			tpf.GetStepFunc(), tpf.GetCheckFunc(), 
			(1e-3*h*h).toDouble());

		std::function<int(const InitCondition<mpreal>&)> evalFunc = 
					 [](const InitCondition<mpreal>& ic) { return sgn(ic.Value - ic.Argument); };
		BisectionComponent<mpreal> bc(thc);
		cout << currentDateTime() << endl;
		bc.DerivativeBisectionGen("0.0", "1.0", "0.0", "1.0", "0.0", "1.0", evalFunc);
		cout << currentDateTime() << endl;
		string fileName = "TroeschProblemSolution.txt";
		thc.SaveToFile(fileName.c_str());
		cout << "Result is saved to " << fileName << endl;
		cout << "Continue? y/n";
		cin >> answer;
	}
}