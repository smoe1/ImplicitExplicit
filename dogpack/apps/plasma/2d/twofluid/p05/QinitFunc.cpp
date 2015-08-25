#include "FiveMomentParams.h"
#include "tensors.h"
#include "DogParamsCart2.h"
#include "PlasmaParams.h"
#include "GEMparams.h"
#include "GEM.h"
#include "gas05.h"
#include "Components.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    void QinitGEM(const dTensor2& xpts, dTensor2& qvals);
    QinitGEM(xpts, qvals);
}

void QinitGEM(const dTensor2& xpts, dTensor2& qvals)
{
    int numpts=xpts.getsize(1);
    const int meqn = qvals.getsize(2);

    const double& gamma = fiveMomentParams.gamma;
    const double& ion_mass = plasmaParams.get_ion_mass();
    double reduced_mass = ion_mass/2.;
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
        //
        double rho_i = ion_mass * number_density;
        double M3_i = reduced_mass*J3; // exact
        //double M3_i = 0.0; // Hakim
	double press_i = ion_pressure_portion*pressure;
        double energy_i = press_i/(gamma-1.0e0) + 0.5e0*(M3_i*M3_i)/rho_i;

        qvals.set(i,_rho_i, rho_i );
        qvals.set(i,_M1_i , 0.   );
        qvals.set(i,_M2_i , 0.   );
        qvals.set(i,_M3_i , M3_i );
        qvals.set(i,_N_i  , energy_i );

        qvals.set(i,_B1   , B1 );
        qvals.set(i,_B2   , B2 );
        qvals.set(i,_E3   , 0. );

        if(_psi <= meqn)
	  qvals.set(i,_psi  , 0.0 );

        if(meqn >=_y_i)
        {
          qvals.set(i,_x_i, rho_i*x);
          qvals.set(i,_y_i, rho_i*y);
        }
    }
}
