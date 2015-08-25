#include "AppSolver.h"
#include "IniDocument.h" // for ini_doc

void AppSolverCart2::init()
{
  if(!&get_state())
  {
    set_state(new AppStateCart2);
    //set_state_old(new AppStateCart2);
    fetch_state().init();
    //fetch_state_old().init();
  }

  DogSolverCart2::init();
}

void AppSolverCart2::initParams()
{
    DogSolverCart2::initParams();
    InitApp(); // get application parameters
}
