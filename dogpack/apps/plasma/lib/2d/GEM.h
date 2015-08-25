#ifndef _GEM_H_
#define _GEM_H_

class GEM{
  public:
    static double get_profile_density(double y);
    static double get_number_density_from_profile_density(double profile_density);
    static double get_pressure_from_profile_density(double profile_density);
    static double get_theta_x(double x);
    static double get_theta_y(double y);
    static void get_B(double theta_x,double theta_y,double y,
        double& B1, double& B2);
    static double get_J3(double theta_x, double theta_y, double y);
    static void get_BJ_largeDomain(double x,double y,
        double& B1, double& B2, double &J3);
    static double get_E2(double y,
             double ion_mass,
             double elc_mass,
             double number_density,
             double ion_pressure_portion,
             double elc_pressure_portion);
};

void increment_steps_that_have_been_accepted();
int get_steps_that_have_been_accepted();

#endif
