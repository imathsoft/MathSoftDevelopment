/***************************************************************************\
|* Function Parser for C++ v4.4.3                                          *|
|*-------------------------------------------------------------------------*|
|* Copyright: Juha Nieminen                                                *|
\***************************************************************************/

#ifndef ONCE_FPARSER_MPFR_H_
#define ONCE_FPARSER_MPFR_H_

#include "fparser.hh"
#include <mpreal.h>
#include <dlmalloc.h>

class FunctionParser_mpfr: public FunctionParserBase<mpfr::mpreal> {};

#endif
