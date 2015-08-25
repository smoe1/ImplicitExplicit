#include "IniDocument.h"
#include "QuadMomentParams.h"

QuadMomentParams quadMomentParams;

void InitApp(IniDocument& ini_doc)
{
  quadMomentParams.init(ini_doc);
}
