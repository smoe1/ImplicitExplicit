#include "IniDocument.h"
#include "AcousticParams.h"
#include "InitApp.h"

AcousticParams acousticParams;

//void AppSolverCart2::initParams()
void InitApp()
{
    acousticParams.init(ini_doc);
}
