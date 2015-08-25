#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

void RunStartScript(int ndims)
{
  char command_str[1024];
  const char* get_outputdir();
  int numchars;
  int exit_status;

  // run startscript
  // to create output directory if it does not exist
  // and copy data files to output directory
  numchars = snprintf(command_str,1024,
		      "if test -f startscript_moments && test -x startscript_moments;\n"
		      "then ./startscript %s %d\n"
		      "else ${DOGPACK}/scripts/startscript_moments %s %d\n"
		      "fi", get_outputdir(),ndims,get_outputdir(),ndims);
  assert_lt(numchars,1023);
  assert_gt(numchars,0);
  exit_status = system(command_str);
}
