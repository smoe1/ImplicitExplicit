#include "dog_ini.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "MHDparams.h"
#include "debug.h"

MHDparams mhdParams;

static void not_set_error(const string& varname)
{
  eprintf("In parameters.ini, section [mhd]"
       ", you must set the parameter %s", varname.c_str());
}

static bool invalid_value(string varname, string val)
{
  eprintf("invalid value: %s = [%s]", varname.c_str(), val.c_str());
  return false;
}

// access defaults from a file that other applications
// (e.g. matlab) might also want to read
void MHDparams::init(IniDocument& ini_doc)
{
  string section_label = "mhd";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  vector<string> option_names_list;
  option_names_list.push_back("gamma");
  option_names_list.push_back("resistivity");
  option_names_list.push_back("viscosity"  );
  option_names_list.push_back("cc"         );

  // get defaults from e.g. $DOGPACK/config/ini_defaults/mhd.ini
  get_defaults      (ini_sec, section_label, option_names_list);
  verify_options_set(ini_sec, section_label, option_names_list);

  string s_gamma       = ini_sec["gamma"      ];
  string s_resistivity = ini_sec["resistivity"];
  string s_viscosity   = ini_sec["viscosity"  ];
  string s_cc          = ini_sec["cc"         ];
  //
  // convert strings to numbers
  //
  istringstream is_gamma       (s_gamma       );
  istringstream is_resistivity (s_resistivity );
  istringstream is_viscosity   (s_viscosity   );
  istringstream is_cc          (s_cc          );
  //
  (is_gamma      >> gamma      )||invalid_value("gamma"      , s_gamma      );
  (is_resistivity>> resistivity)||invalid_value("resistivity", s_resistivity);
  (is_viscosity  >> viscosity  )||invalid_value("viscosity"  , s_viscosity  );
  (is_cc         >> cc         )||invalid_value("cc"         , s_cc         );

  checkParameters();
  reportParameters();
}

void MHDparams::checkParameters()
{
  assert_gt(gamma,0.);
  assert_ge(resistivity,0.);
  assert_ge(viscosity,0.);
  assert_gt(cc,0.);
}

void MHDparams::reportParameters()
{
  // Output parameters to screen
  cout << "   === parameters from [mhd] ===" << endl;
  cout << "   gamma:        " << gamma       << endl;
  cout << "   resistivity:  " << resistivity << endl;
  cout << "   viscosity  :  " << viscosity   << endl;
  cout << "   cc         :  " << cc          << endl;
  cout << endl;
}

