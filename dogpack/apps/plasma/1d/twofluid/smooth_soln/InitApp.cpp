#include "IniDocument.h"
#include "TwoFluidParams.h"
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;

TwoFluidParams twofluidParams;

void InitApp(IniDocument& ini_doc)
{
  twofluidParams.init(ini_doc);

  ofstream massratio_file("mass_ratio.dat", ios::out );
  massratio_file << setprecision(16);
  massratio_file << setw(24) << scientific << twofluidParams.mass_ratio << endl;
}
