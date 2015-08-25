#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "CollisionParams.h"

void CollisionParams::init(IniDocument& ini_doc)
{
  string section_label = "collisionParams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [mhdparams]
  //
  string s_epsilon             = ini_sec["epsilon"         ];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_epsilon           (s_epsilon         );

  //
  // Store information in global parameter values
  //
  epsilon = 0.000000000001;

  is_epsilon >> epsilon;  

  //
  // Output values to screen
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   epsilon             :  " << epsilon           << endl;
  cout << endl;
}

void CollisionParams::write_collisionhelp(const char* filename)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"%24.16e\n", epsilon);
  fclose(file);
} 
