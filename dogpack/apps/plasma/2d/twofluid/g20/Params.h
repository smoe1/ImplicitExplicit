#ifndef Params_h
#define Params_h
namespace Params
{
  enum ProblemType{
    GEM=1,
    SHOCK=2,
    ODE=3,
    STEADY=4,
  };
  ProblemType get_problem();
};
#endif
