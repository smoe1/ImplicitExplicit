#include <stdlib.h>
#include "assert.h"

void RunMatlabCopyScript(const char* outputdir)
{
  char command_str[1024];
  int numchars;
  int exit_status;

  // script to call MATLAB routines if necessary
  numchars = snprintf(command_str,1024,
		      "${DOGPACK}/scripts/matlabcopy_script %s/matlab_output\n",
		      outputdir);
  assert(numchars<1023);
  assert(numchars>0);
  exit_status = system(command_str);
}
