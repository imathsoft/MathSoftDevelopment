#include "AuxUtils.h"

namespace auxutils
{
	void WriteToStream(std::ofstream& stream, float value )
	{
		stream << ToString(value);
	}

	void WriteToStream(std::ofstream& stream, double value )
	{
		stream << ToString(value);
	}

	void WriteToStream(std::ofstream& stream, float_50_noet value )
	{
		stream << value << " ";
	}
};