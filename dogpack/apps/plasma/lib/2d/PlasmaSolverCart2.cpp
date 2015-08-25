#include "PlasmaSolverCart2.h"
#include <string>
#include "GEMparams.h"

void PlasmaSolverCart2::initOutputDirectory(int idx, const char* framedir)
{
  DogSolverCart2::initOutputDirectory(idx, framedir);
  std::string filename = std::string(framedir)+"/out_parameters.ini";
  gemParams.append_to_output_parameters_ini(filename.c_str());
}
