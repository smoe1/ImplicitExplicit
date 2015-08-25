#ifndef App_h
#define App_h

#include "AppStateCart2.h"
#include "PlasmaSolverCart2.h"

// Implementation file is App.cpp

class AppSolverCart2 : public PlasmaSolverCart2
{
 public:
  virtual void init();
  virtual void initParams();

  const AppStateCart2& get_state()
    {return (const AppStateCart2&)DogSolver::get_state();}
  AppStateCart2& fetch_state()
    {return (AppStateCart2&)DogSolver::fetch_state();}
};

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

#endif // App_h
