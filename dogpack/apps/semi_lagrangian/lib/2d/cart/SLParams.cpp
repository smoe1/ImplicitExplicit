#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "SLParams.h"

void SLParams::init(IniDocument& ini_doc)
{

    string section_label = "hybridparams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [mhdparams]
    //
    string s_rk_nystrom_flag   = ini_sec["rk_nystrom_flag"];
    string s_Rk_Cfl            = ini_sec["Rk_Cfl"];
    string s_hybrid_flag       = ini_sec["hybrid_flag"];
//  string s_k                 = ini_sec["k"];

    //
    // Convert read-in strings to numerical values
    //
//  istringstream is_rho0            (s_rho0         );
//  istringstream is_alpha           (s_alpha        );
//  istringstream is_temp            (s_temp         );
    istringstream is_Rk_Cfl          (s_Rk_Cfl       );
    istringstream is_hybrid_flag     (s_hybrid_flag            );
    istringstream is_rk_nystrom_flag (s_rk_nystrom_flag            );
//    istringstream is_k               (s_k            );
//  istringstream is_k               (s_k            );

    //
    // Store information in global parameter values
    //

    // default values .. (?) - DS
    Rk_Cfl          = .1;
    hybrid_flag     = 0;
    rk_nystrom_flag = 0;

    is_Rk_Cfl >> Rk_Cfl;  
    is_rk_nystrom_flag      >> rk_nystrom_flag;  
    is_hybrid_flag      >> hybrid_flag;  
    //is_k      >> k;  

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   Note: many of these are not used for all Semi-Lagrangian applications:" << endl;
    cout << "   rk_nystrom_flag   :  " << rk_nystrom_flag << endl;
    cout << "   hybrid_flag       :  " << hybrid_flag     << endl;
    cout << "   Rk_Cfl            :  " << Rk_Cfl          << endl;
    cout << endl;

}
