#include "tensors.h"
#include "Maxwell.h"
#include "MaxwellParams.h"
#include "PlasmaParams.h"
#include "Relax.h"
#include "Components.h"
#include "DogSolver.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc( const dTensor2& xpts, const dTensor2& qvals,
     const dTensor2& auxvals, dTensor2& source)
{
  source.setall(0.);

  const int numpts=xpts.getsize(1);
  
  // Parameters
  const double one_over_epsilon = maxwellParams.get_cs_light_squared();
  const double cp_speed_squared = maxwellParams.get_cp_speed_squared();
  double charge_density_coef = cp_speed_squared*one_over_epsilon;
  if(!maxwellParams.get_clean_E_field()) charge_density_coef = 0.;
  // Parameters
  const double q_mi_inv  = plasmaParams.get_ion_q_over_m();
  const double q_me_inv  = plasmaParams.get_elc_q_over_m();

  const double time = DogSolver::get_time_hack();
  DogSolver& dogSolver = DogSolver::fetch_solver();
  
  // Loop over each point
  for (int i=1; i<=numpts; i++)
  {
    double x = xpts.get(i,1);
    double y = xpts.get(i,2);

    const double& rho_i = qvals.get(i,_rho_i);
    const double& rho_e = qvals.get(i,_rho_e);
    //
    const double ion_charge_dens = rho_i*q_mi_inv;
    const double elc_charge_dens = rho_e*q_me_inv;
    const double charge_density = ion_charge_dens + elc_charge_dens;
    const int& meqn = source.getsize(2);
    if(meqn >= _psi)
    {
      const double psi = qvals.get(i,_psi);
      const double psi_source = -eps_Bcp_func(psi,time);
      source.set(i,_psi,  psi_source);
    }
    if(meqn >= _phi && charge_density_coef)
    {
      const double phi = qvals.get(i,_phi);
      source.set(i,_phi, charge_density_coef*charge_density
                         -eps_Ecp_func(phi,time));
    }

    // Variables
    const double rho_i_inv = 1./rho_i;
    const double& M1_i  = qvals.get(i,_M1_i );
    const double& M2_i  = qvals.get(i,_M2_i );
    const double& M3_i  = qvals.get(i,_M3_i );
    const double  u1_i  = M1_i*rho_i_inv;
    const double  u2_i  = M2_i*rho_i_inv;
    const double  u3_i  = M3_i*rho_i_inv;
    const double& N11_i = qvals.get(i,_N11_i);
    const double& N12_i = qvals.get(i,_N12_i);
    const double& N13_i = qvals.get(i,_N13_i);
    const double& N22_i = qvals.get(i,_N22_i);
    const double& N23_i = qvals.get(i,_N23_i);
    const double& N33_i = qvals.get(i,_N33_i);
    const double& QN111_i = qvals.get(i,_Q111_i);
    const double& QN112_i = qvals.get(i,_Q112_i);
    const double& QN113_i = qvals.get(i,_Q113_i);
    const double& QN122_i = qvals.get(i,_Q122_i);
    const double& QN123_i = qvals.get(i,_Q123_i);
    const double& QN133_i = qvals.get(i,_Q133_i);
    const double& QN222_i = qvals.get(i,_Q222_i);
    const double& QN223_i = qvals.get(i,_Q223_i);
    const double& QN233_i = qvals.get(i,_Q233_i);
    const double& QN333_i = qvals.get(i,_Q333_i);
    
    const double rho_e_inv = 1./rho_e;
    const double& M1_e  = qvals.get(i,_M1_e );
    const double& M2_e  = qvals.get(i,_M2_e );
    const double& M3_e  = qvals.get(i,_M3_e );
    const double  u1_e  = M1_e*rho_e_inv;
    const double  u2_e  = M2_e*rho_e_inv;
    const double  u3_e  = M3_e*rho_e_inv;
    const double& N11_e = qvals.get(i,_N11_e);
    const double& N12_e = qvals.get(i,_N12_e);
    const double& N13_e = qvals.get(i,_N13_e);
    const double& N22_e = qvals.get(i,_N22_e);
    const double& N23_e = qvals.get(i,_N23_e);
    const double& N33_e = qvals.get(i,_N33_e);
    const double& QN111_e = qvals.get(i,_Q111_e);
    const double& QN112_e = qvals.get(i,_Q112_e);
    const double& QN113_e = qvals.get(i,_Q113_e);
    const double& QN122_e = qvals.get(i,_Q122_e);
    const double& QN123_e = qvals.get(i,_Q123_e);
    const double& QN133_e = qvals.get(i,_Q133_e);
    const double& QN222_e = qvals.get(i,_Q222_e);
    const double& QN223_e = qvals.get(i,_Q223_e);
    const double& QN233_e = qvals.get(i,_Q233_e);
    const double& QN333_e = qvals.get(i,_Q333_e);
    
    const double& B1    = qvals.get(i,_B1);
    const double& B2    = qvals.get(i,_B2);
    const double& B3    = qvals.get(i,_B3);
    const double& E1    = qvals.get(i,_E1);
    const double& E2    = qvals.get(i,_E2);
    const double& E3    = qvals.get(i,_E3);

    // solve the electro-momentum system
    //
    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE);
    {
      const double J1_i = M1_i*q_mi_inv;
      const double J2_i = M2_i*q_mi_inv;
      const double J3_i = M3_i*q_mi_inv;
      const double J1_e = M1_e*q_me_inv;
      const double J2_e = M2_e*q_me_inv;
      const double J3_e = M3_e*q_me_inv;
      const double J1 = J1_i + J1_e;
      const double J2 = J2_i + J2_e;
      const double J3 = J3_i + J3_e;
      source.fetch(i, _M1_i) += ion_charge_dens*E1 + (J2_i*B3 - J3_i*B2);
      source.fetch(i, _M2_i) += ion_charge_dens*E2 + (J3_i*B1 - J1_i*B3);
      source.fetch(i, _M3_i) += ion_charge_dens*E3 + (J1_i*B2 - J2_i*B1);
      source.fetch(i, _M1_e) += elc_charge_dens*E1 + (J2_e*B3 - J3_e*B2);
      source.fetch(i, _M2_e) += elc_charge_dens*E2 + (J3_e*B1 - J1_e*B3);
      source.fetch(i, _M3_e) += elc_charge_dens*E3 + (J1_e*B2 - J2_e*B1);
      source.fetch(i, _E1  ) -= J1*one_over_epsilon;
      source.fetch(i, _E2  ) -= J2*one_over_epsilon;
      source.fetch(i, _E3  ) -= J3*one_over_epsilon;
    }

    // evolve the energy tensor

    // decompose energy into kinetic and thermal components
    //
    const double  K11_i = M1_i*u1_i;
    const double  K12_i = M1_i*u2_i;
    const double  K13_i = M1_i*u3_i;
    const double  K22_i = M2_i*u2_i;
    const double  K23_i = M2_i*u3_i;
    const double  K33_i = M3_i*u3_i;
    //
    const double  P11_i = N11_i-K11_i;
    const double  P12_i = N12_i-K12_i;
    const double  P13_i = N13_i-K13_i;
    const double  P22_i = N22_i-K22_i;
    const double  P23_i = N23_i-K23_i;
    const double  P33_i = N33_i-K33_i;
    //
    const double  K11_e = M1_e*u1_e;
    const double  K12_e = M1_e*u2_e;
    const double  K13_e = M1_e*u3_e;
    const double  K22_e = M2_e*u2_e;
    const double  K23_e = M2_e*u3_e;
    const double  K33_e = M3_e*u3_e;
    //
    const double  P11_e = N11_e-K11_e;
    const double  P12_e = N12_e-K12_e;
    const double  P13_e = N13_e-K13_e;
    const double  P22_e = N22_e-K22_e;
    const double  P23_e = N23_e-K23_e;
    const double  P33_e = N33_e-K33_e;

    // collisions should be handled by operator splitting
    //
    assert(plasmaParams.get_source_type_isoP() != SourceType::SIMPLE);

    // non-collisional evolution of energy tensor
    //
    // (Evolve kinetic energy tensor and rotate the pressure tensor.)
    // Note that exact electromomentum evolution
    // implies exact kinetic energy evolution.
    //
    // evolve kinetic energy
    //
    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE);
    // (q/m)*Sym2(E*M)
    const double Sym2EM_11_i = 2.0*(M1_i*E1);
    const double Sym2EM_12_i = (M1_i*E2+M2_i*E1);
    const double Sym2EM_13_i = (M1_i*E3+M3_i*E1);
    const double Sym2EM_22_i = 2.0*(M2_i*E2);
    const double Sym2EM_23_i = (M2_i*E3+M3_i*E2);
    const double Sym2EM_33_i = 2.0*(M3_i*E3);
    //
    const double Sym2EM_11_e = 2.0*(M1_e*E1);
    const double Sym2EM_12_e = (M1_e*E2+M2_e*E1);
    const double Sym2EM_13_e = (M1_e*E3+M3_e*E1);
    const double Sym2EM_22_e = 2.0*(M2_e*E2);
    const double Sym2EM_23_e = (M2_e*E3+M3_e*E2);
    const double Sym2EM_33_e = 2.0*(M3_e*E3);

    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE
    && plasmaParams.get_source_type_rotP() == SourceType::SIMPLE);

    // (q/m)*Sym2(NxB)
    const double Sym2NxB_11_i = 2.0*(N12_i*B3-N13_i*B2);
    const double Sym2NxB_12_i = (N13_i*B1-N11_i*B3+N22_i*B3-N23_i*B2);
    const double Sym2NxB_13_i = (N11_i*B2-N12_i*B1+N23_i*B3-N33_i*B2);
    const double Sym2NxB_22_i = 2.0*(N23_i*B1-N12_i*B3);
    const double Sym2NxB_23_i = (N12_i*B2-N22_i*B1+N33_i*B1-N13_i*B3);
    const double Sym2NxB_33_i = 2.0*(N13_i*B2-N23_i*B1);
    //
    const double Sym2NxB_11_e = 2.0*(N12_e*B3-N13_e*B2);
    const double Sym2NxB_12_e = (N13_e*B1-N11_e*B3+N22_e*B3-N23_e*B2);
    const double Sym2NxB_13_e = (N11_e*B2-N12_e*B1+N23_e*B3-N33_e*B2);
    const double Sym2NxB_22_e = 2.0*(N23_e*B1-N12_e*B3);
    const double Sym2NxB_23_e = (N12_e*B2-N22_e*B1+N33_e*B1-N13_e*B3);
    const double Sym2NxB_33_e = 2.0*(N13_e*B2-N23_e*B1);
    //
    //
    source.fetch(i,_N11_i)+=q_mi_inv*(Sym2EM_11_i + Sym2NxB_11_i);
    source.fetch(i,_N12_i)+=q_mi_inv*(Sym2EM_12_i + Sym2NxB_12_i);
    source.fetch(i,_N13_i)+=q_mi_inv*(Sym2EM_13_i + Sym2NxB_13_i);
    source.fetch(i,_N22_i)+=q_mi_inv*(Sym2EM_22_i + Sym2NxB_22_i);
    source.fetch(i,_N23_i)+=q_mi_inv*(Sym2EM_23_i + Sym2NxB_23_i);
    source.fetch(i,_N33_i)+=q_mi_inv*(Sym2EM_33_i + Sym2NxB_33_i);
    //
    source.fetch(i,_N11_e)+=q_me_inv*(Sym2EM_11_e + Sym2NxB_11_e);
    source.fetch(i,_N12_e)+=q_me_inv*(Sym2EM_12_e + Sym2NxB_12_e);
    source.fetch(i,_N13_e)+=q_me_inv*(Sym2EM_13_e + Sym2NxB_13_e);
    source.fetch(i,_N22_e)+=q_me_inv*(Sym2EM_22_e + Sym2NxB_22_e);
    source.fetch(i,_N23_e)+=q_me_inv*(Sym2EM_23_e + Sym2NxB_23_e);
    source.fetch(i,_N33_e)+=q_me_inv*(Sym2EM_33_e + Sym2NxB_33_e);

    // Sym3ENijk = Ei*Njk+Ej*Nik+Ek*Nij;
    //
    const double Sym3EN_111_i = 3*E1*N11_i;
    const double Sym3EN_112_i = 2*E1*N12_i + E2*N11_i;
    const double Sym3EN_113_i = 2*E1*N13_i + E3*N11_i;
    const double Sym3EN_122_i = E1*N22_i + 2*E2*N12_i;
    const double Sym3EN_123_i = E1*N23_i + E2*N13_i + E3*N12_i;
    const double Sym3EN_133_i = E1*N33_i + 2*E3*N13_i;
    const double Sym3EN_222_i = 3*E2*N22_i;
    const double Sym3EN_223_i = 2*E2*N23_i + E3*N22_i;
    const double Sym3EN_233_i = E2*N33_i + 2*E3*N23_i;
    const double Sym3EN_333_i = 3*E3*N33_i;
    //
    // Sym3QNxB_ijk
    //
    const double Sym3QNxB_111_i = 3*(QN112_i*B3 - QN113_i*B2);
    const double Sym3QNxB_112_i =    QN113_i*B1 - QN111_i*B3
                                + 2*(QN122_i*B3 - QN123_i*B2);
    const double Sym3QNxB_113_i =    QN111_i*B2 - QN112_i*B1
                                + 2*(QN123_i*B3 - QN133_i*B2);
    const double Sym3QNxB_122_i = 2*(QN123_i*B1 - QN112_i*B3)
                                  + (QN222_i*B3 - QN223_i*B2);
    const double Sym3QNxB_123_i =   (QN112_i*B2 - QN122_i*B1)
                                  + (QN133_i*B1 - QN113_i*B3)
                                  + (QN223_i*B3 - QN233_i*B2);
    const double Sym3QNxB_133_i = 2*(QN113_i*B2 - QN123_i*B1)
                                  + (QN233_i*B3 - QN333_i*B2);
    const double Sym3QNxB_222_i = 3*(QN223_i*B1 - QN122_i*B3);
    const double Sym3QNxB_223_i =   (QN122_i*B2 - QN222_i*B1)
                                + 2*(QN233_i*B1 - QN123_i*B3);
    const double Sym3QNxB_233_i = 2*(QN123_i*B2 - QN223_i*B1)
                                  + (QN333_i*B1 - QN133_i*B3);
    const double Sym3QNxB_333_i = 3*(QN133_i*B2 - QN233_i*B1);

    source.set(i,_Q111_i, q_mi_inv*(Sym3EN_111_i + Sym3QNxB_111_i));
    source.set(i,_Q112_i, q_mi_inv*(Sym3EN_112_i + Sym3QNxB_112_i));
    source.set(i,_Q113_i, q_mi_inv*(Sym3EN_113_i + Sym3QNxB_113_i));
    source.set(i,_Q122_i, q_mi_inv*(Sym3EN_122_i + Sym3QNxB_122_i));
    source.set(i,_Q123_i, q_mi_inv*(Sym3EN_123_i + Sym3QNxB_123_i));
    source.set(i,_Q133_i, q_mi_inv*(Sym3EN_133_i + Sym3QNxB_133_i));
    source.set(i,_Q222_i, q_mi_inv*(Sym3EN_222_i + Sym3QNxB_222_i));
    source.set(i,_Q223_i, q_mi_inv*(Sym3EN_223_i + Sym3QNxB_223_i));
    source.set(i,_Q233_i, q_mi_inv*(Sym3EN_233_i + Sym3QNxB_233_i));
    source.set(i,_Q333_i, q_mi_inv*(Sym3EN_333_i + Sym3QNxB_333_i));

    // Sym3ENijk = Ei*Njk+Ej*Nik+Ek*Nij;
    //
    const double Sym3EN_111_e = 3*E1*N11_e;
    const double Sym3EN_112_e = 2*E1*N12_e + E2*N11_e;
    const double Sym3EN_113_e = 2*E1*N13_e + E3*N11_e;
    const double Sym3EN_122_e = E1*N22_e + 2*E2*N12_e;
    const double Sym3EN_123_e = E1*N23_e + E2*N13_e + E3*N12_e;
    const double Sym3EN_133_e = E1*N33_e + 2*E3*N13_e;
    const double Sym3EN_222_e = 3*E2*N22_e;
    const double Sym3EN_223_e = 2*E2*N23_e + E3*N22_e;
    const double Sym3EN_233_e = E2*N33_e + 2*E3*N23_e;
    const double Sym3EN_333_e = 3*E3*N33_e;
    //
    // Sym3QNxB_ijk
    //
    const double Sym3QNxB_111_e = 3*(QN112_e*B3 - QN113_e*B2);
    const double Sym3QNxB_112_e =    QN113_e*B1 - QN111_e*B3
                                + 2*(QN122_e*B3 - QN123_e*B2);
    const double Sym3QNxB_113_e =    QN111_e*B2 - QN112_e*B1
                                + 2*(QN123_e*B3 - QN133_e*B2);
    const double Sym3QNxB_122_e = 2*(QN123_e*B1 - QN112_e*B3)
                                  + (QN222_e*B3 - QN223_e*B2);
    const double Sym3QNxB_123_e =   (QN112_e*B2 - QN122_e*B1)
                                  + (QN133_e*B1 - QN113_e*B3)
                                  + (QN223_e*B3 - QN233_e*B2);
    const double Sym3QNxB_133_e = 2*(QN113_e*B2 - QN123_e*B1)
                                  + (QN233_e*B3 - QN333_e*B2);
    const double Sym3QNxB_222_e = 3*(QN223_e*B1 - QN122_e*B3);
    const double Sym3QNxB_223_e =   (QN122_e*B2 - QN222_e*B1)
                                + 2*(QN233_e*B1 - QN123_e*B3);
    const double Sym3QNxB_233_e = 2*(QN123_e*B2 - QN223_e*B1)
                                  + (QN333_e*B1 - QN133_e*B3);
    const double Sym3QNxB_333_e = 3*(QN133_e*B2 - QN233_e*B1);

    source.set(i,_Q111_e, q_me_inv*(Sym3EN_111_e + Sym3QNxB_111_e));
    source.set(i,_Q112_e, q_me_inv*(Sym3EN_112_e + Sym3QNxB_112_e));
    source.set(i,_Q113_e, q_me_inv*(Sym3EN_113_e + Sym3QNxB_113_e));
    source.set(i,_Q122_e, q_me_inv*(Sym3EN_122_e + Sym3QNxB_122_e));
    source.set(i,_Q123_e, q_me_inv*(Sym3EN_123_e + Sym3QNxB_123_e));
    source.set(i,_Q133_e, q_me_inv*(Sym3EN_133_e + Sym3QNxB_133_e));
    source.set(i,_Q222_e, q_me_inv*(Sym3EN_222_e + Sym3QNxB_222_e));
    source.set(i,_Q223_e, q_me_inv*(Sym3EN_223_e + Sym3QNxB_223_e));
    source.set(i,_Q233_e, q_me_inv*(Sym3EN_233_e + Sym3QNxB_233_e));
    source.set(i,_Q333_e, q_me_inv*(Sym3EN_333_e + Sym3QNxB_333_e));
  }
}
