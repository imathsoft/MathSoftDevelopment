#ifndef GUARD_UNITTESTAUX
#define GUARD_UNITTESTAUX

//#include <string>

//using namespace std;

#ifdef _DEBUG
	#ifdef WIN32
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/Win32/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/Win32/Release/mpfr.lib")
	#else
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/x64/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/x64/Release/mpfr.lib")
	#endif
#else
	#ifdef WIN32
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/Win32/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/Win32/Release/mpfr.lib")
	#else
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpir/dll/x64/Release/mpir.lib")
		#pragma comment(lib, "../../../../Math_C++/mpfr_mpir_x86_x64_msvc2010/mpfr/dll/x64/Release/mpfr.lib")
	#endif
#endif

namespace UnitTestAux
{
	static wchar_t* Message(const char* text)
	{
		size_t size = strlen(text) + 1;
		wchar_t* wa = new wchar_t[size];
		mbstowcs(wa,text,size);
		return wa;
	}
};

#endif