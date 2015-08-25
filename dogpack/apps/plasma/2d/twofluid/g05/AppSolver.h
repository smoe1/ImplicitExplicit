#ifndef AppSolver_h
#define AppSolver_h

#include "AppSolverCart2.h"

class AppSolver{
  AppSolverCart2 appSolverCart2;
 public:
  int main(int argc, char* argv[])
  {
    // if there were multiple solvers we would need to choose the
    // appropriate solver to delegate to based on mesh_type in
    // parameters.ini
    return appSolverCart2.main(argc,argv);
  }
};

#endif
