#include "AuxUtils.h"

namespace auxutils
{
    void WriteToStream(std::ofstream& stream, mpfr::mpreal value )
	{
		stream << value.toString() << " ";
	}

	void WriteToStream(std::ofstream& stream, float value )
	{
		mpfr::mpreal imv = value;
		WriteToStream(stream, imv);
	}

	void WriteToStream(std::ofstream& stream, double value )
	{
		mpfr::mpreal imv = value;
		WriteToStream(stream, imv);
	}

	void WriteToStream(std::ofstream& stream, float_50_noet value )
	{
		stream << value << " ";
	}
};