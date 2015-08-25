// Borrowed from apps/2d/mhd/OrszagTang (DS)
#include "IniDocument.h"
#include "SLParams.h"

SLParams slParams;

void InitApp(IniDocument& ini_doc)
{
    slParams.init(ini_doc);
}
