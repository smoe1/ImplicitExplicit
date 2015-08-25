#include <cmath>
#include "dogdefs.h"
#include "FiveMomentParams.h"
#include "DogParamsCart2.h"
#include "PlasmaParams.h"
#include "Components.h"
#include "gas05.h"
#include "GEM.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  if(1)
  {
    void QinitGEM(const dTensor2& xpts, dTensor2& qvals);
    QinitGEM(xpts, qvals);
  }
  else
  {
    void QinitShock(const dTensor2& xpts, dTensor2& qvals);
    QinitShock(xpts, qvals);
  }
}

void QinitShock(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);
    const int meqn = qvals.getsize(2);
    const double& gamma = fiveMomentParams.gamma;

    const int x_cycles = 1;
    const int y_cycles = 1;
    const double x_wave_num =
      2*pi*x_cycles/(dogParamsCart2.get_xhigh()-dogParamsCart2.get_xlow());
    const double y_wave_num =
      2*pi*y_cycles/(dogParamsCart2.get_yhigh()-dogParamsCart2.get_ylow());

    // Brio-Wu parameters
    //
    const double& ion_mass = plasmaParams.get_ion_mass();
    const double& elc_mass = plasmaParams.get_elc_mass();
    //
    // left values
    //
    double const& left_rho_i   = 1.0*ion_mass;
    double const& left_u1_i    = 0.0;
    double const& left_u2_i    = 0.0;
    double const& left_u3_i    = 0.0;
    double const& left_press_i = 0.5;
    //
    double const& left_B1      = 0.75;
    double const& left_B2      = 1.0;
    double const& left_B3      = 0.0;
    double const& left_E1      = 0.0;
    double const& left_E2      = 0.0;
    double const& left_E3      = 0.0;	    
    //
    // rght values
    //
    double const& rght_rho_i   = 0.125*ion_mass;
    double const& rght_u1_i    = 0.0;
    double const& rght_u2_i    = 0.0;
    double const& rght_u3_i    = 0.0;
    double const& rght_press_i = 0.05;
    //
    double const& rght_B1      = 0.75;
    double const& rght_B2      =-1.0;
    double const& rght_B3      = 0.0;
    double const& rght_E1      = 0.0;
    double const& rght_E2      = 0.0;
    double const& rght_E3      = 0.0;

    //
    // derived values
    //
    double const left_M1_i = left_rho_i*left_u1_i;
    double const left_M2_i = left_rho_i*left_u2_i;
    double const left_M3_i = left_rho_i*left_u3_i;
    //
    const double& left_P11_i = left_press_i;
    const double  left_P12_i = 0.0;
    const double  left_P13_i = 0.0;
    const double& left_P22_i = left_press_i;
    const double  left_P23_i = 0.0;
    const double& left_P33_i = left_press_i;

    double const rght_M1_i = rght_rho_i*rght_u1_i;
    double const rght_M2_i = rght_rho_i*rght_u2_i;
    double const rght_M3_i = rght_rho_i*rght_u3_i;
    //
    const double& rght_P11_i = rght_press_i;
    const double  rght_P12_i = 0.0;
    const double  rght_P13_i = 0.0;
    const double& rght_P22_i = rght_press_i;
    const double  rght_P23_i = 0.0;
    const double& rght_P33_i = rght_press_i;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
	double y = xpts.get(i,2);
                
        double cycle_progress = (x_wave_num*x + y_wave_num*y)/(2*pi);
        double intpart = floor(cycle_progress);
        double cycle_fraction = cycle_progress - intpart;
        //double cycle_fraction = modf(cycle_progress, &intpart);

        if(cycle_fraction < 0.5)
        {
            qvals.set(i,1,  left_rho_i );
            qvals.set(i,2,  left_M1_i );
            qvals.set(i,3,  left_M2_i );
            qvals.set(i,4,  left_M3_i );
            qvals.set(i,5,  left_M1_i*left_u1_i + left_P11_i); // energy11 (ion)
            qvals.set(i,6,  left_M1_i*left_u2_i + left_P12_i); // energy12 (ion)
            qvals.set(i,7,  left_M1_i*left_u3_i + left_P13_i); // energy13 (ion)
            qvals.set(i,8,  left_M2_i*left_u2_i + left_P22_i); // energy22 (ion)
            qvals.set(i,9,  left_M2_i*left_u3_i + left_P23_i); // energy23 (ion)
            qvals.set(i,10, left_M3_i*left_u3_i + left_P33_i); // energy33 (ion)

	    qvals.set(i,11, left_B1 );
	    qvals.set(i,12, left_B2 );
	    qvals.set(i,13, left_E3 );
	    qvals.set(i,14, 0.0 );
            if(meqn >=_entropy_i)
            {
              const double left_entropy_i
                = get_entropy_per_vol05(left_rho_i,left_press_i);
	      qvals.set(i,_entropy_i, left_entropy_i );
            }
        }
        else
        {
            qvals.set(i,1,  rght_rho_i );
            qvals.set(i,2,  rght_M1_i );
            qvals.set(i,3,  rght_M2_i );
            qvals.set(i,4,  rght_M3_i );
            qvals.set(i,5,  rght_M1_i*rght_u1_i + rght_P11_i); // energy11 (ion)
            qvals.set(i,6,  rght_M1_i*rght_u2_i + rght_P12_i); // energy12 (ion)
            qvals.set(i,7,  rght_M1_i*rght_u3_i + rght_P13_i); // energy13 (ion)
            qvals.set(i,8,  rght_M2_i*rght_u2_i + rght_P22_i); // energy22 (ion)
            qvals.set(i,9,  rght_M2_i*rght_u3_i + rght_P23_i); // energy23 (ion)
            qvals.set(i,10, rght_M3_i*rght_u3_i + rght_P33_i); // energy33 (ion)

	    qvals.set(i,11, rght_B1 );
	    qvals.set(i,12, rght_B2 );
	    qvals.set(i,13, rght_E3 );
	    qvals.set(i,14, 0.0 );
            if(meqn >=_entropy_i)
            {
              const double rght_entropy_i
                = get_entropy_per_vol05(rght_rho_i,rght_press_i);
	      qvals.set(i,_entropy_i, rght_entropy_i );
            }
        }
    }
}

void QinitGEM(const dTensor2& xpts, dTensor2& qvals)
{
    int numpts=xpts.getsize(1);
    const int meqn = qvals.getsize(2);

    double gamma = fiveMomentParams.gamma;
    double ion_mass = plasmaParams.get_ion_mass();
    double elc_mass = plasmaParams.get_elc_mass();
    double reduced_mass = ion_mass*elc_mass/(ion_mass+elc_mass);
    //
    // The pressure portions should add to 1.0.
    //
    // The ratio of the pressures equals the ratio
    // of the temperatures by quasineutrality.
    //
    double ion_pressure_portion = 0.5;
    
    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
	double y = xpts.get(i,2);
        //
        // gas-dynamic
        //
        double profile_density = GEM::get_profile_density(y);
        double number_density =
            GEM::get_number_density_from_profile_density(profile_density);
        double pressure =
            GEM::get_pressure_from_profile_density(profile_density);
        //
        // magnetic field
        //
        double theta_x = GEM::get_theta_x(x);
        double theta_y = GEM::get_theta_y(y);
        double B1, B2;
        GEM::get_B(theta_x,theta_y,y, B1,B2);
        double J3 = GEM::get_J3(theta_x, theta_y, y); // current = curl of B
        //
        // derived gas-dynamic quantities
        double rho_i = ion_mass * number_density;
        double M3_i = reduced_mass*J3; // exact
        //double M3_i = 0.0; // Hakim
	double press_i = ion_pressure_portion*pressure;

        assert(rho_i>EPSILON);
        qvals.set(i,_rho_i, rho_i);
        qvals.set(i,_M1_i , 0.0  );
        qvals.set(i,_M2_i , 0.0  );
        qvals.set(i,_M3_i , M3_i );
        qvals.set(i,_N11_i, press_i);
        qvals.set(i,_N12_i, 0.);
        qvals.set(i,_N13_i, 0.);
        qvals.set(i,_N22_i, press_i);
        qvals.set(i,_N23_i, 0.);
        qvals.set(i,_N33_i, M3_i*M3_i/rho_i + press_i );

        qvals.set(i,_B1   , B1 );
        qvals.set(i,_B2   , B2 );
        qvals.set(i,_E3   , 0. );

        if(_psi <= meqn)
	  qvals.set(i,_psi  , 0. );

        if(_entropy_i <= meqn)
        {
          qvals.set(i,_entropy_i, get_entropy_per_vol05(rho_i,press_i));
        }
    }
}
