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
    
    const double& B1    = qvals.get(i,_B1);
    const double& B2    = qvals.get(i,_B2);
    const double& B3    = qvals.get(i,_B3);
    const double& E1    = qvals.get(i,_E1);
    const double& E2    = qvals.get(i,_E2);
    const double& E3    = qvals.get(i,_E3);

    // solve the electro-momentum system
    //
    if(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE)
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

    // isotropize pressure tensors
    //
    if(plasmaParams.get_source_type_isoP() == SourceType::SIMPLE)
    {
      // the isotropization rate should not be infinity
      // (or the user should be instantaneously relaxing)
      assert_ne(plasmaParams.get_ion_base_iso_period(), 0.);
      assert_ne(plasmaParams.get_elc_base_iso_period(), 0.);
      if(plasmaParams.get_ion_base_iso_rate() > 0.0)
      {
        const double  p_i   = (P11_i + P22_i + P33_i)*(1./3.);
        //
        double ion_iso_rate = Relax::get_iso_rate(
          plasmaParams.get_iso_period_type(),
          plasmaParams.get_ion_base_iso_rate(),
          plasmaParams.get_ion_mass_inv(), rho_i,
          p_i, P11_i, P12_i, P13_i, P22_i, P23_i, P33_i);
        source.fetch(i, _N11_i) += (p_i-P11_i)*ion_iso_rate;
        source.fetch(i, _N12_i) += (   -P12_i)*ion_iso_rate;
        source.fetch(i, _N13_i) += (   -P13_i)*ion_iso_rate;
        source.fetch(i, _N22_i) += (p_i-P22_i)*ion_iso_rate;
        source.fetch(i, _N23_i) += (   -P23_i)*ion_iso_rate;
        source.fetch(i, _N33_i) += (p_i-P33_i)*ion_iso_rate;
      }
      //
      if(plasmaParams.get_elc_base_iso_rate() > 0.0)
      {
        const double  p_e   = (P11_e + P22_e + P33_e)*(1./3.);
        double elc_iso_rate = Relax::get_iso_rate(
          plasmaParams.get_iso_period_type(),
          plasmaParams.get_elc_base_iso_rate(),
          plasmaParams.get_elc_mass_inv(), rho_e,
          p_e, P11_e, P12_e, P13_e, P22_e, P23_e, P33_e);
        source.fetch(i, _N11_e) += (p_e-P11_e)*elc_iso_rate;
        source.fetch(i, _N12_e) += (   -P12_e)*elc_iso_rate;
        source.fetch(i, _N13_e) += (   -P13_e)*elc_iso_rate;
        source.fetch(i, _N22_e) += (p_e-P22_e)*elc_iso_rate;
        source.fetch(i, _N23_e) += (   -P23_e)*elc_iso_rate;
        source.fetch(i, _N33_e) += (p_e-P33_e)*elc_iso_rate;
      }
    }

    // non-collisional evolution of energy tensor
    //
    // (Evolve kinetic energy tensor and rotate the pressure tensor.)
    // Note that exact electromomentum evolution
    // implies exact kinetic energy evolution.
    //
    // evolve kinetic energy
    //
    if(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE)
    {
      source.fetch(i,_N11_i)+=q_mi_inv*2.0*(M1_i*E1);
      source.fetch(i,_N12_i)+=q_mi_inv*(M1_i*E2+M2_i*E1);
      source.fetch(i,_N13_i)+=q_mi_inv*(M1_i*E3+M3_i*E1);
      source.fetch(i,_N22_i)+=q_mi_inv*2.0*(M2_i*E2);
      source.fetch(i,_N23_i)+=q_mi_inv*(M2_i*E3+M3_i*E2);
      source.fetch(i,_N33_i)+=q_mi_inv*2.0*(M3_i*E3);
      //
      source.fetch(i,_N11_e)+=q_me_inv*2.0*(M1_e*E1);
      source.fetch(i,_N12_e)+=q_me_inv*(M1_e*E2+M2_e*E1);
      source.fetch(i,_N13_e)+=q_me_inv*(M1_e*E3+M3_e*E1);
      source.fetch(i,_N22_e)+=q_me_inv*2.0*(M2_e*E2);
      source.fetch(i,_N23_e)+=q_me_inv*(M2_e*E3+M3_e*E2);
      source.fetch(i,_N33_e)+=q_me_inv*2.0*(M3_e*E3);
    }

    // if using exact non collisional energy evolution continue
    if(plasmaParams.get_source_type_Emom() != SourceType::SIMPLE
    && plasmaParams.get_source_type_rotP() != SourceType::SIMPLE)
      continue;

    // declare component of energy tensors not evolved exactly
    double i11=0., e11=0.,
           i12=0., e12=0.,
           i13=0., e13=0.,
           i22=0., e22=0.,
           i23=0., e23=0.,
           i33=0., e33=0.;
    if(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE
    && plasmaParams.get_source_type_rotP() == SourceType::SIMPLE)
    {
      i11=N11_i; e11=N11_e;
      i12=N12_i; e12=N12_e;
      i13=N13_i; e13=N13_e;
      i22=N22_i; e22=N22_e;
      i23=N23_i; e23=N23_e;
      i33=N33_i; e33=N33_e;
    }
    else if(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE
         && plasmaParams.get_source_type_rotP() != SourceType::SIMPLE)
    // the exact solver will handle thermal energy rotation
    // but the source term must handle the contribution of the
    // kinetic energy evolution to the total energy evolution
    {
      i11=K11_i; e11=K11_e;
      i12=K12_i; e12=K12_e;
      i13=K13_i; e13=K13_e;
      i22=K22_i; e22=K22_e;
      i23=K23_i; e23=K23_e;
      i33=K33_i; e33=K33_e;
    }
    else if(plasmaParams.get_source_type_Emom() != SourceType::SIMPLE
         && plasmaParams.get_source_type_rotP() == SourceType::SIMPLE)
    // the exact solver will handle kinetic energy evolution
    // but the source term would need to handle the contribution of the
    // thermal energy to the total energy evolution.
    {
      i11=P11_i; e11=P11_e;
      i12=P12_i; e12=P12_e;
      i13=P13_i; e13=P13_e;
      i22=P22_i; e22=P22_e;
      i23=P23_i; e23=P23_e;
      i33=P33_i; e33=P33_e;
    }
    else
    {
      eprintf("should not get here");
    }
    //
    source.fetch(i,_N11_i)+=q_mi_inv*2.0*(i12*B3-i13*B2);
    source.fetch(i,_N12_i)+=q_mi_inv*(i13*B1-i11*B3+i22*B3-i23*B2);
    source.fetch(i,_N13_i)+=q_mi_inv*(i11*B2-i12*B1+i23*B3-i33*B2);
    source.fetch(i,_N22_i)+=q_mi_inv*2.0*(i23*B1-i12*B3);
    source.fetch(i,_N23_i)+=q_mi_inv*(i12*B2-i22*B1+i33*B1-i13*B3);
    source.fetch(i,_N33_i)+=q_mi_inv*2.0*(i13*B2-i23*B1);
    //
    source.fetch(i,_N11_e)+=q_me_inv*2.0*(e12*B3-e13*B2);
    source.fetch(i,_N12_e)+=q_me_inv*(e13*B1-e11*B3+e22*B3-e23*B2);
    source.fetch(i,_N13_e)+=q_me_inv*(e11*B2-e12*B1+e23*B3-e33*B2);
    source.fetch(i,_N22_e)+=q_me_inv*2.0*(e23*B1-e12*B3);
    source.fetch(i,_N23_e)+=q_me_inv*(e12*B2-e22*B1+e33*B1-e13*B3);
    source.fetch(i,_N33_e)+=q_me_inv*2.0*(e13*B2-e23*B1);
    //
    //source.fetch(i,_N11_i)+=q_mi_inv*2.0*(M1_i*E1+i12*B3-i13*B2);
    //source.fetch(i,_N12_i)+=q_mi_inv*(M1_i*E2+M2_i*E1+i13*B1-i11*B3+i22*B3-i23*B2);
    //source.fetch(i,_N13_i)+=q_mi_inv*(M1_i*E3+M3_i*E1+i11*B2-i12*B1+i23*B3-i33*B2);
    //source.fetch(i,_N22_i)+=q_mi_inv*2.0*(M2_i*E2+i23*B1-i12*B3);
    //source.fetch(i,_N23_i)+=q_mi_inv*(M2_i*E3+M3_i*E2+i12*B2-i22*B1+i33*B1-i13*B3);
    //source.fetch(i,_N33_i)+=q_mi_inv*2.0*(M3_i*E3+i13*B2-i23*B1);
    ////
    //source.fetch(i,_N11_e)+=q_me_inv*2.0*(M1_e*E1+e12*B3-e13*B2);
    //source.fetch(i,_N12_e)+=q_me_inv*(M1_e*E2+M2_e*E1+e13*B1-e11*B3+e22*B3-e23*B2);
    //source.fetch(i,_N13_e)+=q_me_inv*(M1_e*E3+M3_e*E1+e11*B2-e12*B1+e23*B3-e33*B2);
    //source.fetch(i,_N22_e)+=q_me_inv*2.0*(M2_e*E2+e23*B1-e12*B3);
    //source.fetch(i,_N23_e)+=q_me_inv*(M2_e*E3+M3_e*E2+e12*B2-e22*B1+e33*B1-e13*B3);
    //source.fetch(i,_N33_e)+=q_me_inv*2.0*(M3_e*E3+e13*B2-e23*B1);
  }
}
