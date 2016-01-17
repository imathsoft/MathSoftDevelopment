#include "Cannon\TroeschHybridCannon.h"
#include "ShootingSimple\BisectionComponent.h"
#include "Problems\TroeschProblem.h"
#include "MultipleShooting\HybridMultipleShootingComponent.h"
#include "Utils\AuxUtils.h"
#include <time.h>
#include <boost\timer.hpp>
#include <omp.h>

using namespace std;
using namespace auxutils;

//typedef float_50_noet numType;
typedef double numType;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    localtime_s(&tstruct, &now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

int main(){
	char answer = 'y';
	while (answer == 'y')
	{
		cout << "lambda =";
		numType l;
		cin >> l;
		numType h;
		cout << "h = ";
		cin >> h;
		h = max((numType)0.000001, min((numType)0.1, h));
		cout <<"h= " << h<< endl;

		TroeschProblem<numType> tp(l);

		PointSimple<numType> ptLeft;
		ptLeft.Argument  = 0;
		ptLeft.Value  = 0;

		PointSimple<numType> ptRight;
		ptRight.Argument  = 1;
		ptRight.Value  = 1;

		HybridMultipleShootingComponent<numType> HMSComp(tp);

		//cout << currentDateTime() << endl;
		boost::timer t;
		std::vector<InitCondition<numType>> solution = HMSComp.Run(ptLeft, ptRight, h);
		//cout << currentDateTime() << endl;
		cout << t.elapsed() << endl;
		std::string fileName = "f:\\TroeschProblemSolution.txt";
		//auxutils::SaveToFile(solution, fileName.c_str());
		//cout << "Result is saved to " << fileName << endl;
		cout << "Continue? y/n";
		cin >> answer;
	}
}