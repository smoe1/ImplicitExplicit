#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include <cstring>
#include "dog_ini.h"
#include "InitialParams.h"
#include "DogSolver.h"

void InitialParams::init(IniDocument& ini_doc)
{
  string section_label = "initialparams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [initialparams]
  //
  string    s_x0 = ini_sec["x0"];
  string    s_y0 = ini_sec["y0"];
  string    s_z0 = ini_sec["z0"];
  string    s_u  = ini_sec["u"];
  string    s_v  = ini_sec["v"];
  string    s_w  = ini_sec["w"];
  string s_width = ini_sec["width"];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream    is_x0 (s_x0);
  istringstream    is_y0 (s_y0);
  istringstream    is_z0 (s_z0);  
  istringstream is_width (s_width);
  istringstream     is_u (s_u);
  istringstream     is_v (s_v);
  istringstream     is_w (s_w);

  //
  // Store information in global parameter values
  //
  x0 = 0.5;
  y0 = 0.5;
  z0 = 0.5;
  width = 0.3;
  u = 0.0;
  v = 0.0;
  w = 0.0;

  is_x0    >> x0;
  is_y0    >> y0;
  is_z0    >> z0;
  is_width >> width;
  is_u     >> u;
  is_v     >> v;
  is_w     >> w;

  // Output values to screen
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   x0        :  " << x0 << endl;
  cout << "   y0        :  " << y0 << endl;
  cout << "   z0        :  " << z0 << endl;
  cout << "   width     :  " << width << endl;
  cout << "   u         :  " << u  << endl;
  cout << "   v         :  " << v  << endl;
  cout << "   w         :  " << w  << endl;
  cout << endl;

  // Output values to plotting help file
  write_plot_help();
}

// from DogSolver.h
const char* get_outputdir();

void InitialParams::write_plot_help()
{
  stringstream ss;
  string outputdir;
  ss << get_outputdir();
  ss >> outputdir;
  string filename = outputdir + "/qhelp_InitialParams.dat";
  FILE* file = fopen(filename.c_str(),"w");
  fprintf(file,"%16.8e : x0\n", x0);
  fprintf(file,"%16.8e : y0\n", y0);
  fprintf(file,"%16.8e : z0\n", z0);
  fprintf(file,"%16.8e : width\n",width);
  fprintf(file,"%16.8e : u\n", u);
  fprintf(file,"%16.8e : v\n", v);
  fprintf(file,"%16.8e : w\n", w);
  fclose(file);
}
