#include "Params.h"

using namespace Params;

// I define this in the object module rather than in the .h file
// so that only this file needs to be recompiled when modifying
// the parameters.  (To avoid paying the price of function call
// overhead one would need to define and populate a variable.)
Params::ProblemType Params::get_problem() { return Params::GEM; }
