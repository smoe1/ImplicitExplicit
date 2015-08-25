#include <cmath>
#include "dogdefs.h"
#include "FiveMomentParams.h"
#include "DogParamsCart2.h"
#include "PlasmaParams.h"
#include "GEMparams.h"
#include "GEM.h"
#include "gas05.h"
#include "Components.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts,
        dTensor2& qvals)
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
    const double& gamma = fiveMomentParams.get_gamma();

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
    double const& left_rho_e   = 1.0*elc_mass;
    double const& left_u1_e    = 0.0;
    double const& left_u2_e    = 0.0;
    double const& left_u3_e    = 0.0;
    double const& left_press_e = 0.5;
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
    double const& rght_rho_e   = 0.125*elc_mass;
    double const& rght_u1_e    = 0.0;
    double const& rght_u2_e    = 0.0;
    double const& rght_u3_e    = 0.0;
    double const& rght_press_e = 0.05;
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
    double const left_M1_e = left_rho_i*left_u1_e;
    double const left_M2_e = left_rho_i*left_u2_e;
    double const left_M3_e = left_rho_i*left_u3_e;
    //
    double const left_energy_i = left_press_i/(gamma-1.0e0) + 0.5e0*left_rho_i
        *(left_u1_i*left_u1_i + left_u2_i*left_u2_i + left_u3_i*left_u3_i);
    double const left_energy_e = left_press_e/(gamma-1.0e0) + 0.5e0*left_rho_e
        *(left_u1_e*left_u1_e + left_u2_e*left_u2_e + left_u3_e*left_u3_e);
    const double  left_entropy_i=get_entropy_per_vol05(left_rho_i,left_press_i);
    const double  left_entropy_e=get_entropy_per_vol05(left_rho_e,left_press_e);
    //
    double const rght_M1_i = rght_rho_i*rght_u1_i;
    double const rght_M2_i = rght_rho_i*rght_u2_i;
    double const rght_M3_i = rght_rho_i*rght_u3_i;
    //
    double const rght_M1_e = rght_rho_i*rght_u1_e;
    double const rght_M2_e = rght_rho_i*rght_u2_e;
    double const rght_M3_e = rght_rho_i*rght_u3_e;
    //
    double const rght_energy_i = rght_press_i/(gamma-1.0e0) + 0.5e0*rght_rho_i
        *(rght_u1_i*rght_u1_i + rght_u2_i*rght_u2_i + rght_u3_i*rght_u3_i);
    double const rght_energy_e = rght_press_e/(gamma-1.0e0) + 0.5e0*rght_rho_e
        *(rght_u1_e*rght_u1_e + rght_u2_e*rght_u2_e + rght_u3_e*rght_u3_e);
    const double  rght_entropy_i=get_entropy_per_vol05(rght_rho_i,rght_press_i);
    const double  rght_entropy_e=get_entropy_per_vol05(rght_rho_e,rght_press_e);

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
            qvals.set(i,_rho_i,  left_rho_i );
            qvals.set(i,_M1_i ,  left_M1_i );
            qvals.set(i,_M2_i ,  left_M2_i );
            qvals.set(i,_M3_i ,  left_M3_i );
            qvals.set(i,_N_i  ,  left_energy_i );

	    qvals.set(i,_rho_e,  left_rho_e );
            qvals.set(i,_M1_e ,  left_M1_e );
            qvals.set(i,_M2_e ,  left_M2_e );
            qvals.set(i,_M3_e ,  left_M3_e );
            qvals.set(i,_N_e  ,  left_energy_e );

	    qvals.set(i,_B1   , left_B1 );
	    qvals.set(i,_B2   , left_B2 );
	    qvals.set(i,_B3   , left_B3 );
	    qvals.set(i,_E1   , left_E1 );
	    qvals.set(i,_E2   , left_E2 );
	    qvals.set(i,_E3   , left_E3 );
	    qvals.set(i,_psi  , 0.0 );
	    qvals.set(i,_phi  , 0.0 );
            if(meqn >=_entropy_i)
            {
	      qvals.set(i,_entropy_i, left_entropy_i );
	      qvals.set(i,_entropy_e, left_entropy_e );
            }
        }
        else
        {
            qvals.set(i,_rho_i,  rght_rho_i );
            qvals.set(i,_M1_i ,  rght_M1_i );
            qvals.set(i,_M2_i ,  rght_M2_i );
            qvals.set(i,_M3_i ,  rght_M3_i );
            qvals.set(i,_N_i  ,  rght_energy_i );

	    qvals.set(i,_rho_e,  rght_rho_e );
            qvals.set(i,_M1_e ,  rght_M1_e );
            qvals.set(i,_M2_e ,  rght_M2_e );
            qvals.set(i,_M3_e ,  rght_M3_e );
            qvals.set(i,_N_e  , rght_energy_e );

	    qvals.set(i,_B1   , rght_B1 );
	    qvals.set(i,_B2   , rght_B2 );
	    qvals.set(i,_B3   , rght_B3 );
	    qvals.set(i,_E1   , rght_E1 );
	    qvals.set(i,_E2   , rght_E2 );
	    qvals.set(i,_E3   , rght_E3 );
	    qvals.set(i,_psi  , 0.0 );
	    qvals.set(i,_phi  , 0.0 );
            if(meqn >=_entropy_i)
            {
	      qvals.set(i,_entropy_i, rght_entropy_i );
	      qvals.set(i,_entropy_e, rght_entropy_e );
            }
        }
    }
}

void QinitGEM(const dTensor2& xpts, dTensor2& qvals)
{
    int numpts=xpts.getsize(1);
    const int meqn = qvals.getsize(2);

    double gamma = fiveMomentParams.get_gamma();
    double ion_mass = plasmaParams.get_ion_mass();
    double elc_mass = plasmaParams.get_elc_mass();
    double reduced_mass = ion_mass*elc_mass/(ion_mass+elc_mass);
    //
    // The pressure portions should add to 1.0.
    //
    // The ratio of the pressures equals the ratio
    // of the temperatures by quasineutrality.
    //
    double elc_pressure_portion = 1./(1.+gemParams.get_temp_ratio());
    double ion_pressure_portion = gemParams.get_temp_ratio()*elc_pressure_portion;
    
    double B_guide = gemParams.get_B_guide();

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
        //
        // densities
        //
        double rho_i = ion_mass * number_density;
        double rho_e = elc_mass * number_density;
        //
        // momenta
        //
        double M3_i = reduced_mass*J3; // exact
        //double M3_i = 0.0; // Hakim
        double M3_e = -reduced_mass*J3; // exact
        //double M3_e = -elc_mass_portion*J3; // Hakim
        //
        // energy
        //
	double press_i = ion_pressure_portion*pressure;
	double press_e = elc_pressure_portion*pressure;
        double energy_i = press_i/(gamma-1.0e0) + 0.5e0*(M3_i*M3_i)/rho_i;
        double energy_e = press_e/(gamma-1.0e0) + 0.5e0*(M3_e*M3_e)/rho_e;

        double E2 = GEM::get_E2(y, ion_mass, elc_mass,
          number_density, ion_pressure_portion, elc_pressure_portion);

        qvals.set(i,_rho_i, ion_mass * number_density);
        qvals.set(i,_M1_i , 0.0  );
        qvals.set(i,_M2_i , 0.0  );
        qvals.set(i,_M3_i , M3_i );
        qvals.set(i,_N_i  , energy_i );

        qvals.set(i,_rho_e, elc_mass * number_density);
        qvals.set(i,_M1_e , 0.0  );
        qvals.set(i,_M2_e , 0.0  );
        qvals.set(i,_M3_e , M3_e );
        qvals.set(i,_N_e  , energy_e );

        qvals.set(i,_B1   , B1 );
        qvals.set(i,_B2   , B2 );
        qvals.set(i,_B3   , B_guide );
        qvals.set(i,_E1   , 0. );
        qvals.set(i,_E2   , E2 );
        qvals.set(i,_E3   , 0. );

        if(meqn<_psi) continue;
	qvals.set(i,_psi  , 0. );
        if(meqn<_phi) continue;
	qvals.set(i,_phi  , 0. );

        if(meqn<_entropy_i) continue;
        qvals.set(i,_entropy_i, get_entropy_per_vol05(rho_i,press_i));
        qvals.set(i,_entropy_e, get_entropy_per_vol05(rho_e,press_e));
    }
}
