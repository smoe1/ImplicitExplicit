#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "PoissonParams.h"

void PoissonParams::init(IniDocument& ini_doc)
{
    string section_label = "poisson";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [poisson]
    //
    string s_lft_bc_type           = ini_sec["lft_bc_type"          ];
    string s_rgt_bc_type           = ini_sec["rgt_bc_type"          ];
    string s_lft_val               = ini_sec["lft_val"              ];
    string s_rgt_val               = ini_sec["rgt_val"              ];

    //
    // Convert read-in strings to numerical values
    //
    istringstream is_lft_bc_type           (s_lft_bc_type           );
    istringstream is_rgt_bc_type           (s_rgt_bc_type           );
    istringstream is_lft_val               (s_lft_val               );
    istringstream is_rgt_val               (s_rgt_val               );

    //
    // Store information in global parameter values
    //
    is_lft_bc_type    >> lft_bc_type;
    is_rgt_bc_type    >> rgt_bc_type;
    is_lft_val        >> lft_val;
    is_rgt_val        >> rgt_val;

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   lft_bc_type       :  " << lft_bc_type      << endl;
    cout << "   rgt_bc_type       :  " << rgt_bc_type      << endl;    
    cout << "   lft_val           :  " << lft_val          << endl;    
    cout << "   rgt_val           :  " << rgt_val          << endl;
    cout << endl;
}
