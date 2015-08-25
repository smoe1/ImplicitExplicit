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
  string    s_x0  = ini_sec["x0"];
  string    s_y0  = ini_sec["y0"];
  string    s_z0  = ini_sec["z0"];
  string    s_w0  = ini_sec["w0"];
  string    s_u1  = ini_sec["u1"];
  string    s_u2  = ini_sec["u2"];
  string    s_u3  = ini_sec["u3"];
  string    s_u4  = ini_sec["u4"];
  string s_width  = ini_sec["width"];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream    is_x0 (s_x0);
  istringstream    is_y0 (s_y0);
  istringstream    is_z0 (s_z0);  
  istringstream    is_w0 (s_w0);  
  istringstream is_width (s_width);
  istringstream     is_u1 (s_u1);
  istringstream     is_u2 (s_u2);
  istringstream     is_u3 (s_u3);
  istringstream     is_u4 (s_u4);

  //
  // Store information in global parameter values
  //
  x0 = 0.5;
  y0 = 0.5;
  z0 = 0.5;
  w0 = 0.5;
  width = 0.3;
  u1 = 0.0;
  u2 = 0.0;
  u3 = 0.0;
  u4 = 0.0;

  is_x0    >> x0;
  is_y0    >> y0;
  is_z0    >> z0;
  is_w0    >> w0;
  is_width >> width;
  is_u1    >> u1;
  is_u2    >> u2;
  is_u3    >> u3;
  is_u4    >> u4;

  // Output values to screen
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   x0        :  " << x0 << endl;
  cout << "   y0        :  " << y0 << endl;
  cout << "   z0        :  " << z0 << endl;
  cout << "   w0        :  " << w0 << endl;
  cout << "   width     :  " << width << endl;
  cout << "   u1        :  " << u1 << endl;
  cout << "   u2        :  " << u2 << endl;
  cout << "   u3        :  " << u3 << endl;
  cout << "   u4        :  " << u4 << endl;
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
    fprintf(file,"%16.8e : w0\n", w0);
    fprintf(file,"%16.8e : width\n",width);
    fprintf(file,"%16.8e : u1\n", u1);
    fprintf(file,"%16.8e : u2\n", u2);
    fprintf(file,"%16.8e : u3\n", u3);
    fprintf(file,"%16.8e : u4\n", u4);
    fclose(file);

}
