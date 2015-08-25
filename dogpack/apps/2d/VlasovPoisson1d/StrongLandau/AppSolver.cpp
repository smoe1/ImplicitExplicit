#include "IniDocument.h"
#include "VlasovParams.h"
#include "SLParams.h"
#include "AppSolver.h"

VlasovParams vlasovParams;
SLParams     slParams;

void AppSolverCart2::initParams()
{
  DogSolverCart2::initParams();
  vlasovParams.init(ini_doc);
  slParams.init(ini_doc);  
}
