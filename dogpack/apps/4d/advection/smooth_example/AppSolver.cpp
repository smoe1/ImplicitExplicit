#include "AppSolver.h"
#include "IniDocument.h"
#include "InitialParams.h"

InitialParams initialParams;

void AppSolverCart4::initParams()
{
    DogSolverCart4::initParams();
    initialParams.init(ini_doc);
}
