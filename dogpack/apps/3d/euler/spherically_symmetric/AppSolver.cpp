#include "AppSolver.h"
#include "IniDocument.h"
#include "InitialParams.h"
#include "EulerParams.h"

EulerParams eulerParams;
InitialParams initialParams;

void AppSolverCart3::initParams()
{
  DogSolverCart3::initParams();
  initialParams.init(ini_doc);
  eulerParams.init(ini_doc);

  initialParams.set_derived_params(eulerParams.gamma);
}
