#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "VlasovParams.h"

void VlasovParams::init(IniDocument& ini_doc)
{

    string section_label = "vlasovparams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [vlasovparams]
    //
    string s_rho0              = ini_sec["rho0"];
    string s_alpha             = ini_sec["alpha"];
    string s_temp              = ini_sec["temp"];
    string s_k                 = ini_sec["k"];

    //
    // Convert read-in strings to numerical values
    //
    istringstream is_rho0            (s_rho0         );
    istringstream is_alpha           (s_alpha        );
    istringstream is_temp            (s_temp         );
    istringstream is_k               (s_k            );

    //
    // Store information in global parameter values
    //

    // default values
    alpha  = 0.;
    temp   = 1.;
    rho0   = 1.0;
    k      = 0.5;

    is_alpha  >> alpha;  
    is_temp   >> temp;  
    is_rho0   >> rho0;  
    is_k      >> k;  

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   Note: many of these are not used for all VP applications " << endl;
    cout << "   alpha             :  " << alpha         << endl;
    cout << "   temp              :  " << temp          << endl;
    cout << "   rho_0             :  " << rho0          << endl;
    cout << "   k                 :  " << k             << endl;
    cout << endl;
}

const double& VlasovParams::get_integration_constant(void)
{ return integration_constant; }

void VlasovParams::set_integration_constant(double c)
{ integration_constant = c; }

const double& VlasovParams::get_avg_momentum(void)
{ return rho_u0; }

void VlasovParams::set_avg_momemtum(double c)
{ rho_u0 = c; }
