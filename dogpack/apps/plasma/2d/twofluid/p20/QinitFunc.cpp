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
  void QinitGEM(const dTensor2& xpts, dTensor2& qvals);
  QinitGEM(xpts, qvals);
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
        const double u3_i = M3_i/rho_i;
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
        qvals.set(i,_Q111_i, 0.);
        qvals.set(i,_Q112_i, 0.);
        qvals.set(i,_Q113_i, u3_i*press_i);
        qvals.set(i,_Q122_i, 0.);
        qvals.set(i,_Q123_i, 0.);
        qvals.set(i,_Q133_i, 0.);
        qvals.set(i,_Q222_i, 0.);
        qvals.set(i,_Q223_i, u3_i*press_i);
        qvals.set(i,_Q233_i, 0.);
        qvals.set(i,_Q333_i, 3*u3_i*press_i + M3_i*u3_i*u3_i);


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
