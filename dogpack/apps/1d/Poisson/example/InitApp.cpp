#include "IniDocument.h"
#include "PoissonParams.h"

PoissonParams poissonParams;

void InitApp(IniDocument& ini_doc)
{
  poissonParams.init(ini_doc);
}
