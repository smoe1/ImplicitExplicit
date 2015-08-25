#include "tensors.h"
#include "constants.h"
#include <cmath>
#include "PlasmaParams.h"
#include "Maxwell.h"
#include "MaxwellParams.h"
#include "Components.h"
#include "DogState.h"
#include "DogSolver.h"
using namespace std;

// This is a user-supplied routine that sets the source term
void SourceTermFunc( const dTensor2& xpts, const dTensor2& qvals,
     const dTensor2& auxvals, dTensor2& source)
{
    // Parameters
    const double ion_mass = plasmaParams.get_ion_mass();
    const double one_over_epsilon = maxwellParams.get_cs_light_squared();
    const double cp_speed_squared = maxwellParams.get_cp_speed_squared();
    const double slowing_period = plasmaParams.get_slowing_period();
    const double slowing_rate_slope = plasmaParams.get_slowing_rate_slope();
    const double trigger_mui3 = plasmaParams.get_trigger_mui3();
    
    const double time = DogSolver::get_time_hack();

    int numpts=xpts.getsize(1);
    // Loop over each point
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
	const double y = xpts.get(i,2);

        // Variables
        const double& rho_i = qvals.get(i,_rho_i);
        const double& M1_i  = qvals.get(i,_M1_i );
        const double& M2_i  = qvals.get(i,_M2_i );
        const double& M3_i  = qvals.get(i,_M3_i );
        
        const double& B1    = qvals.get(i,_B1);
        const double& B2    = qvals.get(i,_B2);
        const double& E3    = qvals.get(i,_E3);

        double ion_charge_dens = rho_i/ion_mass;
        double J1_i = M1_i/ion_mass;
        double J2_i = M2_i/ion_mass;
        double J3_i = M3_i/ion_mass;
        double J3 = 2*J3_i;

        // compute the temperature to determine
        // the equilibration rate
        double local_slowing_rate = 0.0;
         if(slowing_rate_slope > EPSILON)
         {
             const double mui3 = ion_mass*abs(M3_i)/rho_i;
             if(mui3>trigger_mui3)
             {
               //cout << "at location ("<<x<<","<<y<<") u3_i="<<u3_i<<endl;
               // should probably use something that is smooth,
               // not just piecewise continuous.
               local_slowing_rate = slowing_rate_slope*(mui3-trigger_mui3);
             }
         }
         else if(slowing_period > EPSILON)
         {
             local_slowing_rate = 1./slowing_period;
             // const double energy_i = qvals.get(i,5);
             // const double KE_i = M1_i*M1_i/(2.*rho_i);
             // const double thermal_energy_i = energy_i - KE_i;
             // const double press_i = (2./3.)*thermal_energy_i;
             // const double n_i = rho_i/ion_mass;
             // const double T_i = press_i/n_i;
             // const double T_3 = T_i*T_i*T_i;
             // const double T_3_2 = sqrt(T_3);
             // local_slowing_rate *= (n_i/T_3_2);
         }
        
        // Source term
        source.set(i,_rho_i, 0.0 );
        source.set(i,_M1_i , (-J3_i*B2));
        source.set(i,_M2_i , (J3_i*B1));
        source.set(i,_M3_i , ion_charge_dens*E3 + (J1_i*B2 - J2_i*B1)
                               - local_slowing_rate*M3_i);
        source.set(i,_N_i  , J3_i*E3);
        source.set(i,_B1   , 0.0 );
        source.set(i,_B2   , 0.0 );
        source.set(i,_E3   , -J3*one_over_epsilon);
        const int& meqn = source.getsize(2);
        if(_psi > meqn) continue;
        const double psi = qvals.get(i,_psi);
        source.set(i,_psi,  -eps_Bcp_func(psi,time));
        if(_y_i > meqn) continue;
        source.set(i,_x_i, 0.0);
        source.set(i,_y_i, 0.0);
    }
}
