#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "TwoFluidParams.h"

void TwoFluidParams::init(IniDocument& ini_doc)
{
  string section_label = "twofluid";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [mhdparams]
  //
  string s_gamma = ini_sec["gamma"];
  string s_mass_ratio = ini_sec["mass_ratio"];
  string s_light_speed = ini_sec["light_speed"];
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_gamma(s_gamma);
  istringstream is_mass_ratio(s_mass_ratio);
  istringstream is_light_speed(s_light_speed);

  //
  // Store information in global parameter values
  //
  gamma = 5.0/3.0;
  mass_ratio = 1.0;
  light_speed = 1.0;
  is_gamma >> gamma;
  is_mass_ratio >> mass_ratio;
  is_light_speed >> light_speed;

  //
  // Output values to screen
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   gamma            :  " << gamma << endl;
  cout << "   mass_ratio       :  " << mass_ratio << endl;
  cout << "   light_speed      :  " << light_speed << endl;
  cout << endl;
}
