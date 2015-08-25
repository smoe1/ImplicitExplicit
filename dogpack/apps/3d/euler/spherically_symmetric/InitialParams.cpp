#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include <cstring>
#include <iomanip>
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
  string    s_rhol = ini_sec["rhol"];
  string    s_u1l  = ini_sec["u1l"];
  string    s_u2l  = ini_sec["u2l"];
  string    s_u3l  = ini_sec["u3l"];
  string    s_pl   = ini_sec["pl"];
  string    s_rhor = ini_sec["rhor"];
  string    s_u1r  = ini_sec["u1r"];
  string    s_u2r  = ini_sec["u2r"];
  string    s_u3r  = ini_sec["u3r"];
  string    s_pr   = ini_sec["pr"];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream    is_rhol (s_rhol);
  istringstream    is_u1l  (s_u1l);
  istringstream    is_u2l  (s_u2l);
  istringstream    is_u3l  (s_u3l);
  istringstream    is_pl   (s_pl);
  istringstream    is_rhor (s_rhor);
  istringstream    is_u1r  (s_u1r);
  istringstream    is_u2r  (s_u2r);
  istringstream    is_u3r  (s_u3r);
  istringstream    is_pr   (s_pr);

  //
  // Store information in global parameter values
  //
  rhol = 0.0;
  u1l  = 0.0;
  u2l  = 0.0;
  u3l  = 0.0;
  pl   = 0.0;
  rhor = 0.0;
  u1r  = 0.0;
  u2r  = 0.0;
  u3r  = 0.0;
  pr   = 0.0;

  is_rhol  >> rhol;
  is_u1l   >> u1l;
  is_u2l   >> u2l;
  is_u3l   >> u3l;
  is_pl    >> pl;
  is_rhor  >> rhor;
  is_u1r   >> u1r;
  is_u2r   >> u2r;
  is_u3r   >> u3r;
  is_pr    >> pr;

  // Output values to screen
  cout << scientific;
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   rhol, rhor:  " << setprecision(4) << rhol << ", " << setprecision(4) << rhor << endl;
  cout << "   u1l,  u1r :  " << setprecision(4) << u1l  << ", " << setprecision(4) << u1r  << endl;
  cout << "   u2l,  u2r :  " << setprecision(4) << u2l  << ", " << setprecision(4) << u2r  << endl;
  cout << "   u3l,  u3r :  " << setprecision(4) << u3l  << ", " << setprecision(4) << u3r  << endl;
  cout << "   pl,   pr  :  " << setprecision(4) << pl   << ", " << setprecision(4) << pr   << endl;
  cout << endl;

  // Output values to plotting help file
  write_plot_help();
}

void InitialParams::set_derived_params(double gamma)
{
  m1l = rhol*u1l;
  m2l = rhol*u2l;
  m3l = rhol*u3l;
  El  = pl/(gamma-1.0) + 0.5*rhol*(u1l*u1l+u2l*u2l+u3l*u3l);

  m1r = rhor*u1r;
  m2r = rhor*u2r;
  m3r = rhor*u3r;
  Er  = pr/(gamma-1.0) + 0.5*rhor*(u1r*u1r+u2r*u2r+u3r*u3r);
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
  fprintf(file,"%16.8e : rhol\n", rhol);
  fprintf(file,"%16.8e : u1l\n",  u1l);
  fprintf(file,"%16.8e : u2l\n",  u2l);
  fprintf(file,"%16.8e : u3l\n",  u3l);
  fprintf(file,"%16.8e : pl\n",   pl);
  fprintf(file,"%16.8e : rhor\n", rhor);
  fprintf(file,"%16.8e : u1r\n",  u1r);
  fprintf(file,"%16.8e : u2r\n",  u2r);
  fprintf(file,"%16.8e : u3r\n",  u3r);
  fprintf(file,"%16.8e : pr\n",   pr);
  fclose(file);
}
