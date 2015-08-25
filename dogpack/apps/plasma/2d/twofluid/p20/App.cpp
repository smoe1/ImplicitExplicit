#include "App.h"
#include "tensors.h"
#include "IniDocument.h" // for ini_doc

// ******* section: AppSolver *******
// ******* section: AppSolverCart2 *******

void AppSolverCart2::init()
{
  if(!&get_state())
  {
    set_state(new AppStateCart2);
    fetch_state().init();
  }

  DogSolverCart2::init();
}

void AppSolverCart2::initParams()
{
    DogSolverCart2::initParams();
    InitApp(); // get application parameters
}

