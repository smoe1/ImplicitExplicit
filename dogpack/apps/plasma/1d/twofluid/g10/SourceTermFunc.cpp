#include "tensors.h"
#include "Params.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& source)
{
    const int numpts=xpts.getsize();
    //
    double relax_11_i = 0.0;
    double relax_12_i = 0.0;
    double relax_13_i = 0.0;
    double relax_22_i = 0.0;
    double relax_23_i = 0.0;
    double relax_33_i = 0.0;
    //
    double relax_11_e = 0.0;
    double relax_12_e = 0.0;
    double relax_13_e = 0.0;
    double relax_22_e = 0.0;
    double relax_23_e = 0.0;
    double relax_33_e = 0.0;
    
    // Parameters
    const double& ion_mass      = appParams.ion_mass;
    const double& elc_mass      = appParams.elc_mass;
    const double& debye         = appParams.debye;
    const double& larmor_radius = appParams.larmor_radius;
    const double  ion_gyroradius = larmor_radius*ion_mass;
    const double  elc_gyroradius = larmor_radius*elc_mass;
    const double& ion_rlx_period = appParams.ion_rlx_period;
    const double& elc_rlx_period = appParams.elc_rlx_period;
    
    const double epsilon = debye*debye*larmor_radius;

    // Loop over each point
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);

        // Variables
        const double& rho_i = qvals.get(i,1);
        const double& M1_i  = qvals.get(i,2);
        const double& M2_i  = qvals.get(i,3);
        const double& M3_i  = qvals.get(i,4);
        const double  u1_i  = M1_i/rho_i;
        const double  u2_i  = M2_i/rho_i;
        const double  u3_i  = M3_i/rho_i;
        const double& N11_i = qvals.get(i,5);
        const double& N12_i = qvals.get(i,6);
        const double& N13_i = qvals.get(i,7);
        const double& N22_i = qvals.get(i,8);
        const double& N23_i = qvals.get(i,9);
        const double& N33_i = qvals.get(i,10);
        
        const double& rho_e = qvals.get(i,11);
        const double& M1_e  = qvals.get(i,12);
        const double& M2_e  = qvals.get(i,13);
        const double& M3_e  = qvals.get(i,14);
        const double  u1_e  = M1_e/rho_e;
        const double  u2_e  = M2_e/rho_e;
        const double  u3_e  = M3_e/rho_e;
        const double& N11_e = qvals.get(i,15);
        const double& N12_e = qvals.get(i,16);
        const double& N13_e = qvals.get(i,17);
        const double& N22_e = qvals.get(i,18);
        const double& N23_e = qvals.get(i,19);
        const double& N33_e = qvals.get(i,20);
        
        const double& B1    = qvals.get(i,21);
        const double& B2    = qvals.get(i,22);
        const double& B3    = qvals.get(i,23);
        const double& E1    = qvals.get(i,24);
        const double& E2    = qvals.get(i,25);
        const double& E3    = qvals.get(i,26);

        // if the relaxation period is 0
        // then we will instantaneously relax instead
        //
        // (Perhaps we should change this to instanteously
        // relax if the relaxation period is enough
        // shorter than a time step to make the numerical
        // method unstable.)
        //
        if(ion_rlx_period > 0.0)
        {
           const double  P11_i = N11_i-M1_i*u1_i;
           const double  P12_i = N12_i-M1_i*u2_i;
           const double  P13_i = N13_i-M1_i*u3_i;
           const double  P22_i = N22_i-M2_i*u2_i;
           const double  P23_i = N23_i-M2_i*u3_i;
           const double  P33_i = N33_i-M3_i*u3_i;
           const double  p_i   = (P11_i + P22_i + P33_i)/3.0;
           //
           relax_11_i = (p_i-P11_i)/ion_rlx_period;
           relax_12_i = (   -P12_i)/ion_rlx_period;
           relax_13_i = (   -P13_i)/ion_rlx_period;
           relax_22_i = (p_i-P22_i)/ion_rlx_period;
           relax_23_i = (   -P23_i)/ion_rlx_period;
           relax_33_i = (p_i-P33_i)/ion_rlx_period;
        }
        //
        if(elc_rlx_period > 0.0)
        {
           const double  P11_e = N11_e-M1_e*u1_e;
           const double  P12_e = N12_e-M1_e*u2_e;
           const double  P13_e = N13_e-M1_e*u3_e;
           const double  P22_e = N22_e-M2_e*u2_e;
           const double  P23_e = N23_e-M2_e*u3_e;
           const double  P33_e = N33_e-M3_e*u3_e;
           const double  p_e   = (P11_e + P22_e + P33_e)/3.0;
           //
           relax_11_e = (p_e-P11_e)/elc_rlx_period;
           relax_12_e = (   -P12_e)/elc_rlx_period;
           relax_13_e = (   -P13_e)/elc_rlx_period;
           relax_22_e = (p_e-P22_e)/elc_rlx_period;
           relax_23_e = (   -P23_e)/elc_rlx_period;
           relax_33_e = (p_e-P33_e)/elc_rlx_period;
        }
        
        // Source term -- ion fluid equations
        source.set(i,1,   0.0 );
        source.set(i,2,   rho_i*(E1 + u2_i*B3 - u3_i*B2)/ion_gyroradius );
        source.set(i,3,   rho_i*(E2 + u3_i*B1 - u1_i*B3)/ion_gyroradius );
        source.set(i,4,   rho_i*(E3 + u1_i*B2 - u2_i*B1)/ion_gyroradius );
        source.set(i,5,   relax_11_i + 2.0*(M1_i*E1 + N12_i*B3 - N13_i*B2)/ion_gyroradius );
        source.set(i,6,   relax_12_i + (M1_i*E2 + M2_i*E1 + N13_i*B1 - N11_i*B3 + N22_i*B3 - N23_i*B2)/ion_gyroradius );
        source.set(i,7,   relax_13_i + (M1_i*E3 + M3_i*E1 + N11_i*B2 - N12_i*B1 + N23_i*B3 - N33_i*B2)/ion_gyroradius );
        source.set(i,8,   relax_22_i + 2.0*(M2_i*E2 + N23_i*B1 - N12_i*B3)/ion_gyroradius );
        source.set(i,9,   relax_23_i + (M2_i*E3 + M3_i*E2 + N12_i*B2 - N22_i*B1 + N33_i*B1 - N13_i*B3)/ion_gyroradius );
        source.set(i,10,  relax_33_i + 2.0*(M3_i*E3 + N13_i*B2 - N23_i*B1)/ion_gyroradius );

        // Source term -- electron fluid equations
        source.set(i,11,   0.0 );
        source.set(i,12,  -rho_e*(E1 + u2_e*B3 - u3_e*B2)/elc_gyroradius );
        source.set(i,13,  -rho_e*(E2 + u3_e*B1 - u1_e*B3)/elc_gyroradius );
        source.set(i,14,  -rho_e*(E3 + u1_e*B2 - u2_e*B1)/elc_gyroradius );
        source.set(i,15,  relax_11_e - 2.0*(M1_e*E1 + N12_e*B3 - N13_e*B2)/elc_gyroradius );
        source.set(i,16,  relax_11_e - (M1_e*E2 + M2_e*E1 + N13_e*B1 - N11_e*B3 + N22_e*B3 - N23_e*B2)/elc_gyroradius );
        source.set(i,17,  relax_11_e - (M1_e*E3 + M3_e*E1 + N11_e*B2 - N12_e*B1 + N23_e*B3 - N33_e*B2)/elc_gyroradius );
        source.set(i,18,  relax_11_e - 2.0*(M2_e*E2 + N23_e*B1 - N12_e*B3)/elc_gyroradius );
        source.set(i,19,  relax_11_e - (M2_e*E3 + M3_e*E2 + N12_e*B2 - N22_e*B1 + N33_e*B1 - N13_e*B3)/elc_gyroradius );
        source.set(i,20,  relax_11_e - 2.0*(M3_e*E3 + N13_e*B2 - N23_e*B1)/elc_gyroradius );

        // Source term -- Maxwell equations
        source.set(i,21,  0.0 );
        source.set(i,22,  0.0 );
        source.set(i,23,  0.0 );
        source.set(i,24,  (M1_e/elc_mass - M1_i/ion_mass)/epsilon );
        source.set(i,25,  (M2_e/elc_mass - M2_i/ion_mass)/epsilon );
        source.set(i,26,  (M3_e/elc_mass - M3_i/ion_mass)/epsilon );
    }
}
