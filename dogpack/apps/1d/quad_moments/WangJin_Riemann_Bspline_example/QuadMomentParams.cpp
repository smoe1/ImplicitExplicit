#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include<cstring>
#include<string>
#include "dog_str.h"
#include "dog_ini.h"
#include "QuadMomentParams.h"
#include "tensors.h"
#include <new>

void QuadMomentParams::init(IniDocument& ini_doc)
{
  string section_label = "QuadMomentParams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [mhdparams]
  //
  string s_use_collision       = ini_sec["use_collision"   ];
  string s_epsilon             = ini_sec["epsilon"         ];
  //string s_rho0                = ini_sec["rho0"            ];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_use_collision      (s_use_collision   );
  istringstream is_epsilon            (s_epsilon         );
  //istringstream is_rho0               (s_rho0            );
 
  //
  // Store information in global parameter values
  //
  use_collision = 1;
  epsilon = 0.001;
 
  
  is_use_collision >> use_collision;
  is_epsilon >> epsilon; 
  //is_rho0 >> rho0; 

  //
  // Output values to screen
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   use_collision       :  " << use_collision     << endl;
  cout << "   epsilon             :  " << epsilon           << endl;
  //cout << "   rho_0               :  " << rho0              << endl;
  cout << endl;
}


void QuadMomentParams::write_quadMomenthelp(const char* filename)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"%f\n", epsilon);
  fprintf(file,"%i\n", use_collision);
  //fprintf(file,"%f\n", rho0);
  fclose(file);
} 
