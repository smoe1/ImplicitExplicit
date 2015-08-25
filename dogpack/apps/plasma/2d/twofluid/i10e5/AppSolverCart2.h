#ifndef AppSolverCart2_h
#define AppSolverCart2_h
#include "AppStateCart2.h"
#include "PlasmaSolverCart2.h"

class AppSolverCart2 : public PlasmaSolverCart2
{
 public:
  virtual void init();
  virtual void initParams();

  // accessors
  //
  // (An alternative to such casting would be for AppSolverCart2::init
  // to populate pointers named e.g. AppSolverCart2::state and
  // AppSolverCart2::state_old) and then copy these pointers to
  // DogSolver::state and DogSolver::state_old.)
  //
  const AppStateCart2& get_state()
    {return (const AppStateCart2&)DogSolver::get_state();}
  AppStateCart2& fetch_state()
    {return (AppStateCart2&)DogSolver::fetch_state();}
  const AppStateCart2& get_state_old()
    {return (const AppStateCart2&)DogSolver::get_state_old();}
  AppStateCart2& fetch_state_old()
    {return (AppStateCart2&)DogSolver::fetch_state_old();}
 private:
  virtual void truncate_logfiles(double time)const;
};

#endif // AppSolverCart2_h
