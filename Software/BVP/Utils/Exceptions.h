#ifndef GUARD_EXCEPTIONS
#define GUARD_EXCEPTIONS

#include <stdexcept>

///Not implemented exception
class NotImplementedException : public std::logic_error
{
public:
explicit NotImplementedException()
		: std::logic_error("Function not yet implemented.")
		{	// construct from message string
		}
};

#endif