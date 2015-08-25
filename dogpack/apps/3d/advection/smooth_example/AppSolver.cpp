#include "AppSolver.h"
#include "IniDocument.h"
#include "InitialParams.h"

InitialParams initialParams;

void AppSolverCart3::initParams()
{
  DogSolverCart3::initParams();
  initialParams.init(ini_doc);
}
