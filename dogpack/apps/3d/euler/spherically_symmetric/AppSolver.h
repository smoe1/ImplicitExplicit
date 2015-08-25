#ifndef AppSolver_h
#define AppSolver_h
#include "DogSolverCart3.h"

class AppSolverCart3 : public DogSolverCart3
{
 public:
  void initParams();
  
};

class AppSolver
{
  
  AppSolverCart3 appSolverCart3;
  
 public:
  
  int main(int argc, char* argv[])
  {
    // if there were multiple solvers we would need to choose the
    // appropriate solver to delegate to based on mesh_type in
    // parameters.ini
    return appSolverCart3.main(argc,argv);
  }
  
};

#endif // AppSolver_h
