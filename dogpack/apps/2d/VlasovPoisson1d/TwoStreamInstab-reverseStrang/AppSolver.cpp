#include "IniDocument.h"
#include "VlasovParams.h"
#include "SLParams.h"
#include "AppSolver.h"
#include "VPDataCart2.h"

VlasovParams vlasovParams;
SLParams     slParams;
VPDataCart2  vpDataCart2;

void AppSolverCart2::initParams()
{
  DogSolverCart2::initParams();
  vlasovParams.init(ini_doc);
  slParams.init(ini_doc);  
  vpDataCart2.init();
}
