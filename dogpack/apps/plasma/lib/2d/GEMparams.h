#ifndef _GEMPARAMS_H_
#define _GEMPARAMS_H_
#include "Boundaries.h"

class IniDocument;
class GEMparams{
    friend class GEM;
    // primary parameters
    //
    // GEM problem is on domain [-4*pi..4pi]x[-2*pi..2*pi].
    // multiply width and length of GEM problem by domain_scaling 
    // so that values and derivatives match in vicinity of origin.
    double domain_scaling; // period_rescale_factor
    double half_sheet_thickness;
    double B_0;
    double B_guide;
    double n_0; // typical number density of both species
    double n_b; // background number density per n_0.
    double temp_ratio;
    int enforce_symmetry;
    int enforcing_flip_symmetry;
    int enforcing_rotational_symmetry;
    int enforced_symmetry;
    bool mesh_is_virtual;
    BoundaryConditionType yBCs;
    BoundaryConditionType xBCs;
    //
    // derived parameters
    //
    double x_angular_frequency;
    double y_angular_frequency;
    double psi_0; // GEM's psi_0
    double psi_0_rescaled;

  private:
    void resetDogParams();
    void set_domain_scaling(double in);
    void set_B_0(double in, double psi_0_over_B_0);
    void set_B_guide(double in);
  public:
    void init(IniDocument& ini_doc);
    double get_domain_scaling()const{return domain_scaling;}
    double get_half_sheet_thickness()const{return half_sheet_thickness;}
    double get_B_0()const{return B_0;}
    double get_B_guide()const{return B_guide;}
    double get_n_0()const{return n_0;}
    double get_n_b()const{return n_b;}
    double get_temp_ratio(){return temp_ratio;}
    //int get_enforce_symmetry()const{return enforce_symmetry;}
    int get_enforcing_flip_symmetry()const{return enforcing_flip_symmetry;}
    int get_enforcing_rotational_symmetry()
      const{return enforcing_rotational_symmetry;}
    int get_enforced_symmetry(){return enforced_symmetry;}
    bool x_symmetry_not_enforced()const{return !(enforced_symmetry&X);}
    bool y_symmetry_not_enforced()const{return !(enforced_symmetry&Y);}
    bool get_mesh_is_virtual(){return mesh_is_virtual;}
    int get_num_quadrants()const{ switch(enforced_symmetry&XY){
      case X:case Y:return 2;case XY: return 1;} return 4;}
    BoundaryConditionType get_yBCs()const{return yBCs;};
    BoundaryConditionType get_xBCs()const{return xBCs;};
    void append_to_output_parameters_ini(const char* filename) const;
};
extern GEMparams gemParams;

#endif
