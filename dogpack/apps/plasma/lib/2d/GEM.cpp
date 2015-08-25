#include<assert.h>
#include<iostream>
#include<sstream>
#include <cmath>
#include "dog_ini.h"
#include "DogParamsCart2.h"
#include "GEMparams.h"
#include "GEM.h"
#include "debug.h"
using namespace std;

// GEM parameters

GEMparams gemParams;

void GEMparams::set_domain_scaling(double in)
{
    x_angular_frequency = M_PI/dogParamsCart2.get_xhigh();
    y_angular_frequency = M_PI/(2.*dogParamsCart2.get_yhigh());

    // assumes L_x = 2*L_y
    //
    domain_scaling = in;
    if(true || domain_scaling != 1.)
    {
      const double angular_period = 4.0*domain_scaling;
      const double angular_frequency = 1./angular_period;
      assert_almost_eq(angular_frequency,x_angular_frequency);
      assert_almost_eq(angular_frequency,y_angular_frequency);
    }
}

void GEMparams::set_B_0(double in, double psi_0_over_B_0)
{
    B_0 = in;
    // psi_0 = 0.1*B_0;
    psi_0 = psi_0_over_B_0 * B_0;
    // rescale psi_0 so that J3|_Xpoint remains unchanged
    psi_0_rescaled = psi_0*domain_scaling*domain_scaling;
}

void GEMparams::set_B_guide(double in)
{
    B_guide = in;
}

static void not_set_error(const string& varname)
{
  eprintf("In parameters.ini, section [plasma]"
       ", you must set the parameter %s", varname.c_str());
}

static bool invalid_value(string varname, string val)
{
  eprintf("invalid value: %s = [%s]", varname.c_str(), val.c_str());
  return false;
}

void GEMparams::append_to_output_parameters_ini(const char* filename) const
{
  FILE* file = fopen(filename,"a");
  fprintf(file,"[GEM]\n");
  fprintf(file,"neg_x_is_virtual = %d\n",!x_symmetry_not_enforced());
  fprintf(file,"neg_y_is_virtual = %d\n",!y_symmetry_not_enforced());
  fprintf(file,"enforcing_rotational_symmetry = %d\n",
    get_enforcing_rotational_symmetry());
  fclose(file);
}

void GEMparams::init(IniDocument& ini_doc)
{
  string section_label = "GEM";
  IniDocument::Section& ini_sec = ini_doc[section_label];
  vector<string> option_names_list;
  //
  option_names_list.push_back("temp_ratio"           );
  option_names_list.push_back("B_0"                  );
  option_names_list.push_back("psi_0_over_B_0"       );
  option_names_list.push_back("B_guide"              );
  option_names_list.push_back("n_0"                  );
  option_names_list.push_back("n_b"                  );
  option_names_list.push_back("yBCs"                 );
  option_names_list.push_back("xBCs"                 );
  option_names_list.push_back("domain_scaling"       );
  option_names_list.push_back("sheet_thickness"      );
  option_names_list.push_back("enforced_symmetry"    );
  option_names_list.push_back("enforce_symmetry"     );
  option_names_list.push_back("enforce_flip_symmetry");
  option_names_list.push_back("mesh_is_virtual"     );

  // get defaults from e.g. $DOGPACK/config/ini_defaults/GEM.ini
  get_defaults      (ini_sec, section_label, option_names_list);
  verify_options_set(ini_sec, section_label, option_names_list);

  string s_temp_ratio            = ini_sec["temp_ratio"           ];
  string s_B_0                   = ini_sec["B_0"                  ];
  string s_psi_0_over_B_0        = ini_sec["psi_0_over_B_0"       ];
  string s_B_guide               = ini_sec["B_guide"              ];
  string s_n_0                   = ini_sec["n_0"                  ];
  string s_n_b                   = ini_sec["n_b"                  ];
  string s_yBCs                  = ini_sec["yBCs"                 ];
  string s_xBCs                  = ini_sec["xBCs"                 ];
  string s_domain_scaling        = ini_sec["domain_scaling"       ];
  string s_sheet_thickness       = ini_sec["sheet_thickness"      ];
  string s_enforced_symmetry     = ini_sec["enforced_symmetry"    ];
  string s_enforce_symmetry      = ini_sec["enforce_symmetry"     ];
  string s_enforce_flip_symmetry = ini_sec["enforce_flip_symmetry"];
  string s_mesh_is_virtual       = ini_sec["mesh_is_virtual"      ];

  istringstream is_temp_ratio            (s_temp_ratio            );
  istringstream is_B_0                   (s_B_0                   );
  istringstream is_psi_0_over_B_0        (s_psi_0_over_B_0        );
  istringstream is_B_guide               (s_B_guide               );
  istringstream is_n_0                   (s_n_0                   );
  istringstream is_n_b                   (s_n_b                   );
  istringstream is_domain_scaling        (s_domain_scaling        );
  istringstream is_sheet_thickness       (s_sheet_thickness       );
  istringstream is_enforced_symmetry     (s_enforced_symmetry     );
  istringstream is_enforce_symmetry      (s_enforce_symmetry      );
  istringstream is_enforce_flip_symmetry (s_enforce_flip_symmetry );
  istringstream is_mesh_is_virtual       (s_mesh_is_virtual       );

  double sheet_thickness;
  double psi_0_over_B_0;
  //
  (is_temp_ratio    >> temp_ratio    )||invalid_value("temp_ratio"    , s_temp_ratio    );
  (is_domain_scaling>> domain_scaling)||invalid_value("domain_scaling", s_domain_scaling);
  (is_sheet_thickness>> sheet_thickness)||invalid_value("sheet_thickness", s_sheet_thickness);
  (is_B_0           >> B_0           )||invalid_value("B_0"           , s_B_0           );
  (is_psi_0_over_B_0>> psi_0_over_B_0)||invalid_value("psi_0_over_B_0", s_psi_0_over_B_0);
  (is_B_guide       >> B_guide       )||invalid_value("B_guide"       , s_B_guide       );
  (is_n_0           >> n_0           )||invalid_value("n_0"           , s_n_0           );
  (is_n_b           >> n_b           )||invalid_value("n_b"           , s_n_b           );
  (is_enforced_symmetry >> enforced_symmetry)
    ||invalid_value("enforced_symmetry", s_enforced_symmetry);
  (is_enforce_symmetry >> enforce_symmetry)
    ||invalid_value("enforce_symmetry", s_enforce_symmetry);
  int enforce_flip_symmetry;
  (is_enforce_flip_symmetry >> enforce_flip_symmetry)
    ||invalid_value("enforce_flip_symmetry", s_enforce_flip_symmetry);
  int int_mesh_is_virtual;
  (is_mesh_is_virtual >> int_mesh_is_virtual)
    ||invalid_value("mesh_is_virtual", s_mesh_is_virtual);
  mesh_is_virtual = !!int_mesh_is_virtual;
  enforcing_flip_symmetry = enforce_flip_symmetry;
  if(B_guide!=0.) enforcing_flip_symmetry = 0.;
  const int use_X_symmetry = enforce_symmetry;
  const int use_Y_symmetry = enforcing_flip_symmetry;
  // except for restarts from old data enforced_symmetry should always be negative
  if(enforced_symmetry < 0)
  {
    enforcing_rotational_symmetry = (!enforcing_flip_symmetry) && enforce_symmetry;
    enforced_symmetry = use_X_symmetry + 2*use_Y_symmetry + 4*(!mesh_is_virtual);
  }
  else
  {
    Wprintf("overriding values of enforce_symmetry,\n"
      " enforce_flip_symmetry, and mesh_is_virtual\n"
      " with value of deprecated parameter enforced_symmetry=%d",enforced_symmetry);
    enforce_symmetry = enforced_symmetry&1;
    enforcing_flip_symmetry = !!(enforced_symmetry&2);
    enforcing_rotational_symmetry = (!enforcing_flip_symmetry) && enforce_symmetry;
    mesh_is_virtual = !(enforced_symmetry&4);
  }

  if(s_yBCs=="open")
    yBCs = OPEN_BOUNDARY;
  else if(s_yBCs=="copy")
    yBCs = COPY_BOUNDARY;
  else if(s_yBCs=="conducting_wall")
    yBCs = CONDUCTING_WALL;
  else if(s_yBCs=="forced")
    yBCs = FORCED;
  else
    invalid_value_error(s_yBCs.c_str());

  if(s_xBCs=="open")
    xBCs = OPEN_BOUNDARY;
  else if(s_xBCs=="copy")
    xBCs = COPY_BOUNDARY;
  else if(s_xBCs=="periodic")
    xBCs = PERIODIC;
  else
    invalid_value_error(s_xBCs.c_str());
  
  //int int_yBCs;
  //int int_xBCs;
  //yBCs = (BoundaryConditionType)int_yBCs;
  //xBCs = (BoundaryConditionType)int_xBCs;

  bool resetDogParamsCart2(double domain_scaling, int enforced_symmetry);
  bool changed_dims = resetDogParamsCart2(domain_scaling, enforced_symmetry);
  if(changed_dims)
  {
    cout << "   Reset parameters from [grid] based on enforced_symmetry="
         << enforced_symmetry << " to be as follows:" << endl;
    dogParamsCart2.reportParameters();
  }

  set_domain_scaling(domain_scaling);
  half_sheet_thickness = sheet_thickness/2.;

  set_B_0(B_0, psi_0_over_B_0);
  set_B_guide(B_guide);

  // report parameters
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   temp_ratio        :  " << temp_ratio       << endl;
  cout << "   domain_scaling    :  " << domain_scaling   << endl;    
  cout << "   sheet_thickness   :  " << sheet_thickness  << endl;    
  cout << "   B_0               :  " << B_0              << endl;          
  cout << "   B_guide           :  " << B_guide          << endl;          
  cout << "   n_0               :  " << n_0              << endl;          
  cout << "   n_b               :  " << n_b              << endl;          
  cout << "   yBCs              :  " << yBCs             << endl;
  cout << "   xBCs              :  " << xBCs             << endl;
  cout << "   --- parameters derived from [" << section_label << "] ---" << endl;
  cout << "   x_angular_frequency :  " << x_angular_frequency << endl;
  cout << "   y_angular_frequency :  " << y_angular_frequency << endl;
  cout << "   psi_0             :  " << psi_0             << endl;
  cout << "   psi_0_rescaled    :  " << psi_0_rescaled    << endl;
  cout << "   enforced_symmetry :  " << enforced_symmetry << endl;
  cout << "   enforcing_rotational_symmetry :  " << enforcing_rotational_symmetry << endl;
  cout << "   enforcing_flip_symmetry :  " << enforcing_flip_symmetry << endl;
  cout << "   mesh_is_virtual   :  " << mesh_is_virtual << endl;
  cout << endl;

  // check parameters
  //
  assert(temp_ratio>0);
  assert(sheet_thickness>0);
  assert(B_0 > 0);
  assert(yBCs > 0);
  assert(xBCs > 0);
}

// GEM methods

// second component of electric field
// (rather than 0 to avoid massive oscillations)
double GEM::get_E2(double y,
  double ion_mass,
  double elc_mass,
  double number_density,
  double ion_pressure_portion,
  double elc_pressure_portion)
{
  double total_mass = ion_mass+elc_mass;
  double tanh_y_l = tanh(y/gemParams.half_sheet_thickness);
  double E2 = (ion_mass*elc_pressure_portion - 
               elc_mass*ion_pressure_portion)
               * gemParams.B_0*gemParams.B_0 / (total_mass*number_density)
               * tanh_y_l*(1-tanh_y_l*tanh_y_l);
  return E2;
}

double GEM::get_profile_density(double y)
{
    return gemParams.n_b + (1.0-pow(tanh(y/gemParams.half_sheet_thickness),2));
}

double GEM::get_number_density_from_profile_density(double profile_density)
{
    return(profile_density*gemParams.n_0);
}

double GEM::get_pressure_from_profile_density(double profile_density)
{
    return 0.5*gemParams.B_0*gemParams.B_0*profile_density;
}

double GEM::get_theta_x(double x)
{
    return x*gemParams.x_angular_frequency;
}
double GEM::get_theta_y(double y)
{
    return y*gemParams.y_angular_frequency;
}

void GEM::get_B(double theta_x,double theta_y,double y,
    double& B1, double& B2)
{
    //double amplitude = gemParams.angular_frequency*gemParams.psi_0_rescaled;
    double amplitude_1 = gemParams.y_angular_frequency*gemParams.psi_0_rescaled;
    double amplitude_2 = gemParams.x_angular_frequency*gemParams.psi_0_rescaled;
    double B1_perturb= -amplitude_1* cos(theta_x)*sin(theta_y);
    double B2_perturb=  amplitude_2* sin(theta_x)*cos(theta_y);
    double B1_equilib = gemParams.B_0 * tanh(y/gemParams.half_sheet_thickness);
    B1 = B1_perturb + B1_equilib;
    B2 = B2_perturb;
}

double GEM::get_J3(double theta_x, double theta_y, double y)
{
    double J3_equilib = -(gemParams.B_0/gemParams.half_sheet_thickness)
        * (1.0-pow(tanh(y/gemParams.half_sheet_thickness),2));
    double amplitude_scale = pow(gemParams.x_angular_frequency,2)
                           + pow(gemParams.y_angular_frequency,2);
    const double amplitude = gemParams.psi_0_rescaled*amplitude_scale;

    // assumes L_x = 2*L_y
    // (by design current at origin is invariant under domain_scaling)
    if(true || gemParams.get_domain_scaling() != 1.)
    {
      const double amp2 = gemParams.psi_0/8.0;
      assert_almost_eq(amp2,amplitude);
    }

    double J3_perturb = amplitude*cos(theta_x)*cos(theta_y);
    double J3 = J3_equilib + J3_perturb;
    //double J3 = J3_equilib; // Hakim
    return J3;
}

// use perturbation that does not reference Lx and Ly
void GEM::get_BJ_largeDomain(double x,double y,
    double& B1, double& B2, double &J3)
{
    const double w0 = gemParams.half_sheet_thickness;
    const double wx = w0*8.;
    const double wy = w0*4.;
    const double w0inv = 1./w0;
    const double wxinv = 1./wx;
    const double wyinv = 1./wy;
    const double wxi2 = wxinv*wxinv;
    const double wyi2 = wyinv*wyinv;
    const double xscaled = x*wxinv;
    const double yscaled = y*wyinv;
    const double xs2 = xscaled*xscaled;
    const double ys2 = yscaled*yscaled;
    const double window = .1*w0*exp(-xs2-ys2);
    const double B1_perturb= -2.*yscaled*window;
    const double B2_perturb=  2.*xscaled*window;
    const double B1_equilib = gemParams.B_0 * tanh(yscaled);
    B1 = B1_perturb + B1_equilib;
    B2 = B2_perturb;
    //double J3 = GEM::get_J3(theta_x, theta_y, y); // current = curl of B
    const double J3_equilib = -(gemParams.B_0*w0inv)*(1.0-pow(tanh(yscaled),2));
    const double J3_perturb = -2.*w0*window*((2*xs2-1)*wxi2+(2*ys2-1)*wyi2);
    J3 = J3_equilib + J3_perturb;
    //J3 = J3_equilib; // Hakim
}

