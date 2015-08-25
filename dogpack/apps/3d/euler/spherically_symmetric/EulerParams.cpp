#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include <cstring>
#include <iomanip>
#include "dog_ini.h"
#include "EulerParams.h"
#include "DogSolver.h"

void EulerParams::init(IniDocument& ini_doc)
{
  string section_label = "eulerparams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [eulerparams]
  //
  string    s_gamma = ini_sec["gamma"];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream    is_gamma (s_gamma);

  //
  // Store information in global parameter values
  //
  gamma = 5.0/3.0;

  is_gamma    >> gamma;

  // Output values to screen
  cout << scientific;
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   gamma     :  " << setprecision(4) << gamma << endl;
  cout << endl;

  // Output values to plotting help file
  write_plot_help();
}

// from DogSolver.h
const char* get_outputdir();

void EulerParams::write_plot_help()
{
  stringstream ss;
  string outputdir;
  ss << get_outputdir();
  ss >> outputdir;
  string filename = outputdir + "/qhelp_EulerParams.dat";
  FILE* file = fopen(filename.c_str(),"w");
  fprintf(file,"%16.8e : gamma\n", gamma);
  fclose(file);
}
