#ifndef AppSolver_h
#define AppSolver_h
#include "DogSolverCart4.h"

class AppSolverCart4 : public DogSolverCart4
{
 public:
  void initParams();
  
};

class AppSolver
{
  
  AppSolverCart4 appSolverCart4;
  
 public:
  
  int main(int argc, char* argv[])
  {
    // if there were multiple solvers we would need to choose the
    // appropriate solver to delegate to based on mesh_type in
    // parameters.ini
    return appSolverCart4.main(argc,argv);
  }
  
};

#endif // AppSolver_h
