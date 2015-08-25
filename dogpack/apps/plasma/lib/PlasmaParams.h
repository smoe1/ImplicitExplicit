#ifndef PLASMAPARAMS_H
#define PLASMAPARAMS_H
#include "FiveMomentParams.h"
#include "MaxwellParams.h"
#include "TwoFluidParams.h"
using namespace std;

#define MAX_NUM_SUBSYSTEMS 10
namespace SpeciesType{
  enum Enum{
    ION = 1,
    ELC = 2,
  };
}

namespace GasModelType{
  enum Enum{
      gas05    = 1,
      gas10    = 2,
      gas20    = 3,
  };
}
namespace PlasmaModelType{
  enum Enum{
      mhd    = 1,
      g05    = 2,
      g10    = 3,
      g20    = 4,
      i10e5  = 5,
      p05    = 6,
      p10    = 7,
      p20    = 8,
  };
}
namespace SourceType{
  enum Enum{
      SIMPLE = 1,
      SPLIT  = 2,
  };
}
namespace IsoPeriodType{
  enum Enum{
      constant = 1,
      trace  = 2,
      det    = 3,
  };
}

class IniDocument;
struct PlasmaParams{
 private:
   //FiveMomentParams fiveMomentParams;
   //MaxwellParams maxwellParams;
   //TwoFluidParams twoFluidParams;
   double mass_ratio           ;
   char*  spc_mass_mode        ;
   double total_mass           ;
   PlasmaModelType::Enum model ;
   SourceType::Enum source_type_Emom;
   SourceType::Enum source_type_rotP;
   SourceType::Enum source_type_isoP;
   IsoPeriodType::Enum iso_period_type;
   double ion_prandtl_number  ;
   double elc_prandtl_number  ;
   double ion_base_iso_period  ;
   double elc_base_iso_period  ;
   double ion_iso_period  ; // no external access
   double elc_iso_period  ; // no external access
   double base_iso_period      ;
   double ion_heat_conductivity;
   double elc_heat_conductivity;
   double ion_lat_iso_timescale;
   double elc_lat_iso_timescale;
   double slowing_period       ;
   double trigger_mui3         ;
   double slowing_rate_slope   ;
   int    allow_fast_sound     ;
   int num_riemann_subsystems; int* riemann_subsystem_firstindices;
   //
   // derived parameters
   //
   double ion_mass;
   double elc_mass;
   double ion_mass_inv;
   double elc_mass_inv;
   double ion_q_over_m;
   double elc_q_over_m;
   double min_speed;
   double max_speed;
   double ion_base_iso_rate;
   double elc_base_iso_rate;
 //
 // methods
 //
 private:
   void checkParameters();
   void reportParameters();
 public:

   PlasmaParams():spc_mass_mode(0),
     num_riemann_subsystems(0), riemann_subsystem_firstindices(0)
   {}
   ~PlasmaParams();
   void init(IniDocument& ini_doc);

   double get_mass_ratio            () const {return mass_ratio           ;}
   const char* get_spc_mass_mode    () const {return spc_mass_mode        ;}
   double get_total_mass            () const {return total_mass           ;}
   PlasmaModelType::Enum get_model  () const {return model                ;}
   SourceType::Enum get_source_type_Emom () const {return source_type_Emom;}
   SourceType::Enum get_source_type_rotP () const {return source_type_rotP;}
   SourceType::Enum get_source_type_isoP () const {return source_type_isoP;}
   IsoPeriodType::Enum get_iso_period_type () const {return iso_period_type;}
   const char* get_iso_period_type_string() const;
   double get_ion_prandtl_number    () const {return ion_prandtl_number   ;}
   double get_elc_prandtl_number    () const {return elc_prandtl_number   ;}
   double get_ion_base_iso_period   () const {return ion_base_iso_period  ;}
   double get_elc_base_iso_period   () const {return elc_base_iso_period  ;}
   double get_ion_base_iso_rate     () const {return ion_base_iso_rate    ;}
   double get_elc_base_iso_rate     () const {return elc_base_iso_rate    ;}
   double get_base_iso_period       () const {return base_iso_period      ;}
   double get_ion_heat_conductivity () const {return ion_heat_conductivity;}
   double get_elc_heat_conductivity () const {return elc_heat_conductivity;}
   double get_ion_lat_iso_timescale () const {return ion_lat_iso_timescale;}
   double get_elc_lat_iso_timescale () const {return elc_lat_iso_timescale;}
   double get_slowing_period        () const {return slowing_period       ;}
   double get_trigger_mui3          () const {return trigger_mui3         ;}
   double get_slowing_rate_slope    () const {return slowing_rate_slope   ;}
   int    get_allow_fast_sound      () const {return allow_fast_sound     ;}
   double get_ion_mass  () const {return ion_mass ;}
   double get_elc_mass  () const {return elc_mass ;}
   double get_ion_mass_inv  () const {return ion_mass_inv ;}
   double get_elc_mass_inv  () const {return elc_mass_inv ;}
   double get_ion_q_over_m() const {return ion_q_over_m;}
   double get_elc_q_over_m() const {return elc_q_over_m;}
   double get_min_speed () const {return min_speed;}
   double get_max_speed () const {return max_speed;}
   //
   // subsystem access
   //
   int get_num_riemann_subsystems()
     const { return num_riemann_subsystems; }
   int get_riemann_subsystem_firstindex(int i) const;
   int get_riemann_subsystem_lastindex(int i) const;
   int get_riemann_subsystem_size(int i) const;
};

extern PlasmaParams plasmaParams;

#endif
