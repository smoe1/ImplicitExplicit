#include <stdlib.h>
#include <cmath>
#include <string>
#include <assert.h>
#include "debug.h"
#include "constants.h"
#include "dog_ini.h"
#include "dog_str.h"
#include "dog_io.h"
#include "PlasmaParams.h"
#include "DogParams.h"
#include <cstring>

PlasmaParams plasmaParams;

PlasmaParams::~PlasmaParams()
{
  free(spc_mass_mode);
  delete [] riemann_subsystem_firstindices;
}

static void not_set_error(const char* varname)
{
  eprintf("In parameters.ini, section [plasma]"
       ", you must set the parameter %s", varname);
}

static bool invalid_value(const char* varname, const char* val)
{
  eprintf("invalid value: %s = [%s]", varname, val);
  return false;
}

static SourceType::Enum get_source_type(string source_type_s)
{
  SourceType::Enum source_type;
  if(source_type_s=="simple")
    source_type = SourceType::SIMPLE;
  else if(source_type_s=="split")
    source_type = SourceType::SPLIT;
  else
    invalid_value_error(source_type_s.c_str());
  return source_type;
}

// access defaults from a file that other applications
// (e.g. matlab) might also want to read
void PlasmaParams::init(IniDocument& ini_doc)
{
  string section_label = "plasma";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  vector<string> option_names_list;
  option_names_list.push_back("gamma"                );
  option_names_list.push_back("mass_ratio"           );
  option_names_list.push_back("spc_mass_mode"        );
  option_names_list.push_back("total_mass"           );
  option_names_list.push_back("ion_mass"             );
  option_names_list.push_back("cs_light"             );
  option_names_list.push_back("model_name"           );
  option_names_list.push_back("source_type_Emom"     );
  option_names_list.push_back("source_type_rotP"     );
  option_names_list.push_back("source_type_isoP"     );
  option_names_list.push_back("ion_prandtl_number"   );
  option_names_list.push_back("elc_prandtl_number"   );
  option_names_list.push_back("iso_period_type"      );
  option_names_list.push_back("ion_iso_period"       );
  option_names_list.push_back("elc_iso_period"       );
  option_names_list.push_back("base_iso_period"      );
  option_names_list.push_back("ion_heat_conductivity");
  option_names_list.push_back("elc_heat_conductivity");
  option_names_list.push_back("ion_lat_iso_timescale");
  option_names_list.push_back("elc_lat_iso_timescale");
  option_names_list.push_back("slowing_period"       );
  option_names_list.push_back("trigger_mui3"         );
  option_names_list.push_back("slowing_rate_slope"   );
  option_names_list.push_back("allow_fast_sound"     );
  option_names_list.push_back("cc"                   );
  option_names_list.push_back("eps_Bcp"              );
  option_names_list.push_back("eps_Ecp"              );
  option_names_list.push_back("clean_E_field"        );
  option_names_list.push_back("riemann_subsystem_firstindices");

  // get defaults from e.g. $DOGPACK/config/ini_defaults/plasma.ini
  get_defaults      (ini_sec, section_label, option_names_list);
  verify_options_set(ini_sec, section_label, option_names_list);

  const char* s_gamma                 = ini_sec["gamma"                ].c_str();
  const char* s_mass_ratio            = ini_sec["mass_ratio"           ].c_str();
  string      spc_mass_mode_s         = ini_sec["spc_mass_mode"        ];
  const char* s_total_mass            = ini_sec["total_mass"           ].c_str();
  const char* s_ion_mass              = ini_sec["ion_mass"             ].c_str();
  string        model_name_s          = ini_sec["model_name"           ];
  const char* s_cs_light              = ini_sec["cs_light"             ].c_str();
  string        source_type_Emom_s    = ini_sec["source_type_Emom"     ];
  string        source_type_rotP_s    = ini_sec["source_type_rotP"     ];
  string        source_type_isoP_s    = ini_sec["source_type_isoP"     ];
  string        iso_period_type_s     = ini_sec["iso_period_type"      ];
  const char* s_ion_prandtl_number    = ini_sec["ion_prandtl_number"   ].c_str();
  const char* s_elc_prandtl_number    = ini_sec["elc_prandtl_number"   ].c_str();
  const char* s_ion_iso_period        = ini_sec["ion_iso_period"       ].c_str();
  const char* s_elc_iso_period        = ini_sec["elc_iso_period"       ].c_str();
  const char* s_base_iso_period       = ini_sec["base_iso_period"      ].c_str();
  const char* s_ion_heat_conductivity = ini_sec["ion_heat_conductivity"].c_str();
  const char* s_elc_heat_conductivity = ini_sec["elc_heat_conductivity"].c_str();
  const char* s_ion_lat_iso_timescale = ini_sec["ion_lat_iso_timescale"].c_str();
  const char* s_elc_lat_iso_timescale = ini_sec["elc_lat_iso_timescale"].c_str();
  const char* s_slowing_period        = ini_sec["slowing_period"       ].c_str();
  const char* s_trigger_mui3          = ini_sec["trigger_mui3"         ].c_str();
  const char* s_slowing_rate_slope    = ini_sec["slowing_rate_slope"   ].c_str();
  const char* s_allow_fast_sound      = ini_sec["allow_fast_sound"     ].c_str();
  const char* s_cc                    = ini_sec["cc"                   ].c_str();
  const char* s_eps_Bcp               = ini_sec["eps_Bcp"              ].c_str();
  const char* s_eps_Ecp               = ini_sec["eps_Ecp"              ].c_str();
  const char* s_clean_E_field         = ini_sec["clean_E_field"        ].c_str();
  const char* s_riemann_subsystem_firstindices = ini_sec["riemann_subsystem_firstindices"  ].c_str();
  //
  // populate class with parameter data.
  //
  double gamma     ;
  double cs_light  ;
  double cc        ;
  double eps_Bcp   ;
  double eps_Ecp   ;
  int clean_E_field;
  //
  sscanf(s_gamma                ,"%lf",&gamma                )|| invalid_value("gamma"                 , s_gamma                 );
  sscanf(s_mass_ratio           ,"%lf",&mass_ratio           )|| invalid_value("mass_ratio"            , s_mass_ratio            );
  sscanf(s_total_mass           ,"%lf",&total_mass           )|| invalid_value("total_mass"            , s_total_mass            );
  sscanf(s_ion_mass             ,"%lf",&ion_mass             )|| invalid_value("ion_mass"              , s_ion_mass              );
  sscanf(s_cs_light             ,"%lf",&cs_light             )|| invalid_value("cs_light"              , s_cs_light              );
  sscanf(s_ion_prandtl_number   ,"%lf",&ion_prandtl_number   )|| invalid_value("ion_prandtl_number"    , s_ion_prandtl_number    );
  sscanf(s_elc_prandtl_number   ,"%lf",&elc_prandtl_number   )|| invalid_value("elc_prandtl_number"    , s_elc_prandtl_number    );
  sscanf(s_ion_iso_period       ,"%lf",&ion_iso_period       )|| invalid_value("ion_iso_period"        , s_ion_iso_period        );
  sscanf(s_elc_iso_period       ,"%lf",&elc_iso_period       )|| invalid_value("elc_iso_period"        , s_elc_iso_period        );
  sscanf(s_base_iso_period      ,"%lf",&base_iso_period      )|| invalid_value("base_iso_period"       , s_base_iso_period       );
  sscanf(s_ion_heat_conductivity,"%lf",&ion_heat_conductivity)|| invalid_value("ion_heat_conductivity" , s_ion_heat_conductivity );
  sscanf(s_elc_heat_conductivity,"%lf",&elc_heat_conductivity)|| invalid_value("elc_heat_conductivity" , s_ion_heat_conductivity );
  sscanf(s_ion_lat_iso_timescale,"%lf",&ion_lat_iso_timescale)|| invalid_value("ion_lat_iso_timescale" , s_ion_lat_iso_timescale );
  sscanf(s_elc_lat_iso_timescale,"%lf",&elc_lat_iso_timescale)|| invalid_value("elc_lat_iso_timescale" , s_elc_lat_iso_timescale );
  sscanf(s_slowing_period       ,"%lf",&slowing_period       )|| invalid_value("slowing_period"        , s_slowing_period        );
  sscanf(s_trigger_mui3         ,"%lf",&trigger_mui3         )|| invalid_value("trigger_mui3"          , s_trigger_mui3          );
  sscanf(s_slowing_rate_slope   ,"%lf",&slowing_rate_slope   )|| invalid_value("slowing_rate_slope"    , s_slowing_rate_slope    );
  sscanf(s_allow_fast_sound     ,"%d" ,&allow_fast_sound     )|| invalid_value("allow_fast_sound"      , s_allow_fast_sound      );
  sscanf(s_cc                   ,"%lf",&cc                   )|| invalid_value("cc"                    , s_cc                    );
  sscanf(s_eps_Bcp              ,"%lf",&eps_Bcp              )|| invalid_value("eps_Bcp"               , s_eps_Bcp               );
  sscanf(s_eps_Ecp              ,"%lf",&eps_Ecp              )|| invalid_value("eps_Ecp"               , s_eps_Ecp               );
  sscanf(s_clean_E_field        ,"%d" ,&clean_E_field        )|| invalid_value("clean_E_field"         , s_clean_E_field         );
  //
  // parse lists of numbers
  //
  const int meqn = dogParams.get_meqn();
  num_riemann_subsystems = new_array_from_str(
    s_riemann_subsystem_firstindices,riemann_subsystem_firstindices,1,',');
  assert_ge(num_riemann_subsystems,1);
  assert_le(num_riemann_subsystems,MAX_NUM_SUBSYSTEMS);
  assert_eq(riemann_subsystem_firstindices[1],1);
  for(int i=1;i<num_riemann_subsystems;i++)
    assert_lt(riemann_subsystem_firstindices[i],riemann_subsystem_firstindices[i+1]);
  for(int i=1;i<=num_riemann_subsystems;i++)
  {
    if(riemann_subsystem_firstindices[i]>meqn)
    {
      num_riemann_subsystems = i-1;
      Wprintf("ignoring riemann subsystem %d with first index %d > meqn = %d;"
        "\n\ttruncating number of riemann subsystems to %d",
        i,riemann_subsystem_firstindices[i],meqn,
        num_riemann_subsystems);
      break;
    }
  }
  assert_le(riemann_subsystem_firstindices[num_riemann_subsystems],meqn);

  fiveMomentParams.gamma      = gamma;
  maxwellParams.set_cs_light(cs_light);
  maxwellParams.set_cp_speed(cc*cs_light);
  maxwellParams.set_eps_Bcp(eps_Bcp);
  maxwellParams.set_eps_Ecp(eps_Ecp);
  maxwellParams.set_clean_E_field(clean_E_field);
  //
  using namespace PlasmaModelType;
  if(model_name_s=="mhd")
    model = mhd;
  else if(model_name_s=="g05")
    model = g05;
  else if(model_name_s=="g10")
    model = g10;
  else if(model_name_s=="g20")
    model = g20;
  else if(model_name_s=="i10e5")
    model = i10e5;
  else if(model_name_s=="p05")
    model = p05;
  else if(model_name_s=="p10")
    model = p10;
  else if(model_name_s=="p20")
    model = p20;
  else
    invalid_value_error(model_name_s.c_str());

  // set source_type
  source_type_Emom = ::get_source_type(source_type_Emom_s);
  source_type_rotP = ::get_source_type(source_type_rotP_s);
  source_type_isoP = ::get_source_type(source_type_isoP_s);

  // set iso_period_type
  {
    using namespace IsoPeriodType;
    if(iso_period_type_s=="constant")
      iso_period_type = constant;
    else if(iso_period_type_s=="trace")
      iso_period_type = trace;
    else if(iso_period_type_s=="det")
      iso_period_type = det;
    else if(iso_period_type_s=="default")
    {
      if(base_iso_period >= 0.)
        iso_period_type = det;
      else
        iso_period_type = constant;
    }
    else
      invalid_value_error(iso_period_type_s.c_str());
  }

  // set derived parameters
  //
  switch(iso_period_type)
  {
    using namespace IsoPeriodType;
    case det:
    case trace:
      ion_base_iso_period = base_iso_period;
      elc_base_iso_period = base_iso_period;
      break;
    case constant:
      ion_base_iso_period = ion_iso_period;
      elc_base_iso_period = elc_iso_period;
      break;
    default:
      invalid_value_error(iso_period_type);
  }
  ion_base_iso_rate = 1./ion_base_iso_period;
  elc_base_iso_rate = 1./elc_base_iso_period;
  // assume that the *sum* of the particle masses is 1.
  spc_mass_mode = strdup(spc_mass_mode_s.c_str());
  if(spc_mass_mode_s=="total_mass")
  {
    elc_mass = total_mass/(1.+mass_ratio);
    ion_mass = mass_ratio*elc_mass;
  }
  else if(spc_mass_mode_s=="ion_mass")
  {
    assert(ion_mass > EPSILON);
    elc_mass = ion_mass/mass_ratio;
    total_mass = ion_mass + elc_mass;
  }
  else
  {
    invalid_value_error(spc_mass_mode);
  }
  ion_mass_inv = 1./ion_mass;
  elc_mass_inv = 1./elc_mass;
  ion_q_over_m = ion_mass_inv;
  elc_q_over_m = -elc_mass_inv;
  assert_ge(cc,1.);
  min_speed =-maxwellParams.get_cp_speed();
  max_speed = maxwellParams.get_cp_speed();

  checkParameters();
  reportParameters();
}

void PlasmaParams::checkParameters()
{
  assert_gt(fiveMomentParams.gamma, 0.);
  assert_gt(maxwellParams.get_cs_light(), 0.);
  assert_gt(maxwellParams.get_cp_speed(), 0.);
  assert_gt(mass_ratio, 0.);
  //
  assert(ion_mass  > 0);
  assert(elc_mass  > 0);
  assert(min_speed < 0);
  assert(max_speed > 0);

  assert(ion_heat_conductivity >= 0);
  assert(elc_heat_conductivity >= 0);

  bool test_equal(double a, double b);
  if(!(test_equal(ion_heat_conductivity*ion_heat_conductivity*ion_mass,
                  elc_heat_conductivity*elc_heat_conductivity*elc_mass)))
  {
    Wprintf("warning: ion heat conductivity should equal"
            " sqrt(1/mass_ratio) times electron heat conductivity:"
            "\n  ion_heat_conductivity = %f"
            "\n  elc_heat_conductivity = %f"
            "\n  mass_ratio = %f",
            ion_heat_conductivity,elc_heat_conductivity,mass_ratio)
  }
  if(source_type_isoP==SourceType::SIMPLE)
  {
    assert_ne(ion_base_iso_period,0.);
    assert_ne(elc_base_iso_period,0.);
  }
  if(iso_period_type==IsoPeriodType::det)
  {
    assert_lt(ion_iso_period,0.);
    // we might want to isotropize electrons instantaneously
    // (we essentially do in the five-moment model)
    assert_le(elc_iso_period,0.);
  }
  if(ion_iso_period < 0)
  {
    if(elc_iso_period>=0)
    {
      Wprintf("ion_iso_period < 0 yet elc_iso_period = %f", elc_iso_period);
    }
  }
  else
  {
    if(!test_equal(ion_iso_period, sqrt(mass_ratio) * elc_iso_period))
    {
      Wprintf("warning: ion isotropization period should equal"
           " sqrt(mass_ratio) times electron isotropization period:"
           "\n  ion_iso_period = %f"
           "\n  elc_iso_period = %f"
           "\n  mass_ratio = %f",
           ion_iso_period, elc_iso_period, mass_ratio);
    }
  }
  bool completely_splitting_source = 
   get_source_type_Emom()==SourceType::SPLIT &&
   get_source_type_rotP()==SourceType::SPLIT &&
   get_source_type_isoP()==SourceType::SPLIT;
  if(dogParams.get_source_term())
  {
    if(completely_splitting_source)
    {
      Wprintf("parameters.ini: you can set source_term = 0");
    }
  }
  else
  {
    if(!completely_splitting_source)
    {
      eprintf("parameters.ini: you need to set source_term = 1");
    }
  }
    
}

inline const char* get_SourceType_str(
  SourceType::Enum sourceType)
{
  switch(sourceType)
  {
    using namespace SourceType;
    case SIMPLE:
      return "simple";
    case SPLIT:
      return "split";
    default:
      invalid_value_error(sourceType);
  }
}

const char* PlasmaParams::get_iso_period_type_string() const
{
  switch(iso_period_type)
  {
    using namespace IsoPeriodType;
    case det:
      return "det";
    case trace:
      return "trace";
    case constant:
      return "constant";
    default:
      invalid_value_error(iso_period_type);
  }
}

void PlasmaParams::reportParameters()
{
  // Output parameters to screen
  printf("   === parameters from [plasma]\n");
  printf("   gamma                 : %g\n", fiveMomentParams.gamma);
  printf("   mass_ratio            : %g\n", mass_ratio);
  printf("   spc_mass_mode         : %s\n", spc_mass_mode);
  printf("   total_mass            : %g\n", total_mass);
  printf("   ion_mass              : %g\n", ion_mass);
  printf("   cs_light              : %g\n", maxwellParams.get_cs_light());
  printf("   source_type_Emom      : %s\n", get_SourceType_str(source_type_Emom));
  printf("   source_type_rotP      : %s\n", get_SourceType_str(source_type_rotP));
  printf("   source_type_isoP      : %s\n", get_SourceType_str(source_type_isoP));
  printf("   iso_period_type       : %s\n", get_iso_period_type_string());
  printf("   ion_prandtl_number    : %g\n", ion_prandtl_number);
  printf("   elc_prandtl_number    : %g\n", elc_prandtl_number);
  printf("   ion_base_iso_period   : %g\n", ion_base_iso_period);
  printf("   elc_base_iso_period   : %g\n", elc_base_iso_period);
  printf("   base_iso_period       : %g\n", base_iso_period);
  printf("   ion_heat_conductivity : %g\n", ion_heat_conductivity);
  printf("   elc_heat_conductivity : %g\n", elc_heat_conductivity);
  printf("   ion_lat_iso_timescale : %g\n", ion_lat_iso_timescale);
  printf("   elc_lat_iso_timescale : %g\n", elc_lat_iso_timescale);
  printf("   slowing_period        : %g\n", slowing_period);
  printf("   trigger_mui3          : %g\n", trigger_mui3);
  printf("   slowing_rate_slope    : %g\n", slowing_rate_slope);
  printf("   allow_fast_sound      : %d\n", allow_fast_sound);
  printf("   cp_speed              : %g\n", maxwellParams.get_cp_speed());
  printf("   eps_Bcp               : %g\n", maxwellParams.get_eps_Bcp());
  printf("   eps_Ecp               : %g\n", maxwellParams.get_eps_Ecp());
  printf("   clean_E_field         : %d\n", maxwellParams.get_clean_E_field());
  printf("   --- parameters derived from [plasma] ---\n");
  printf("   ion_mass              : %g\n", ion_mass);
  printf("   elc_mass              : %g\n", elc_mass);
  printf("   min_speed             : %g\n", min_speed);
  printf("   max_speed             : %g\n", max_speed);
  if(num_riemann_subsystems>1){
    printf(
        "    num_riemann_subsystems:  %d\n"
        "riemann_subsystem_firstindices:  ",
      num_riemann_subsystems);
    fprint_array(stdout,riemann_subsystem_firstindices,1,num_riemann_subsystems);
  }
  printf("\n");
}

int PlasmaParams::get_riemann_subsystem_firstindex(int i) const
{
  assert_ge(i,1);
  assert_le(i,num_riemann_subsystems);
  return riemann_subsystem_firstindices[i];
}

int PlasmaParams::get_riemann_subsystem_lastindex(int i) const
{
  assert_ge(i,1);
  assert_le(i,num_riemann_subsystems);
  if(i<num_riemann_subsystems)
    return riemann_subsystem_firstindices[i+1] - 1;
  else
    return dogParams.get_meqn();
}

int PlasmaParams::get_riemann_subsystem_size(int i) const
{
  return get_riemann_subsystem_lastindex(i) - get_riemann_subsystem_firstindex(i) + 1;
}

