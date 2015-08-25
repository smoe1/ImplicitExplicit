#include "tensors.h"
#include "Params.h"
#include "IniDocument.h"
#include <fstream>
#include <iostream>
using namespace std;

AppParams appParams;

void nextline(ifstream& read_file, int n=1)
{
    char buffer[256];
    for(int i=0;i<n;i++) read_file.getline(buffer,256);
}

void getAppParams()
{
    ifstream read_file("param.data", ios::in);
    //
    // Configured Application Parameters
    //
                         read_file >> appParams.gamma;
    nextline(read_file); read_file >> appParams.mass_ratio;
    nextline(read_file); read_file >> appParams.debye;
    nextline(read_file); read_file >> appParams.cs_light;
    nextline(read_file); read_file >> appParams.larmor_radius;
    nextline(read_file); read_file >> appParams.ion_rlx_period;
    nextline(read_file); read_file >> appParams.elc_rlx_period;
    read_file.close();
    //
    // set derived parameters
    //
    // assume that the *sum* of the particle masses is 1.
    appParams.ion_mass = appParams.mass_ratio/(1.+appParams.mass_ratio);
    //appParams.ion_mass = 1.0; // Hakim?
    appParams.elc_mass = appParams.ion_mass/appParams.mass_ratio;
    cout << "appParams.gamma:          " << appParams.gamma           << endl;
    cout << "appParams.mass_ratio:     " << appParams.mass_ratio      << endl;
    cout << "appParams.debye:          " << appParams.debye           << endl;
    cout << "appParams.cs_light:       " << appParams.cs_light        << endl;
    cout << "appParams.larmor_radius:  " << appParams.larmor_radius   << endl;
    cout << "appParams.ion_rlx_period: " << appParams.ion_rlx_period  << endl;
    cout << "appParams.elc_rlx_period: " << appParams.elc_rlx_period  << endl;
    cout << endl;
    cout << "appParams.ion_mass:       " << appParams.ion_mass        << endl;  
    cout << "appParams.elc_mass:       " << appParams.elc_mass        << endl;  
    cout << endl;
}

void InitApp(IniDocument& ini_doc)
{
    getAppParams();
}
