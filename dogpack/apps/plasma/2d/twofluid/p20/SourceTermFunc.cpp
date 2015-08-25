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
  //const double q_me_inv  = plasmaParams.get_elc_q_over_m();

  const double time = DogSolver::get_time_hack();
  const double dt_inv = 1./DogSolver::get_dt_hack();
  DogSolver& dogSolver = DogSolver::fetch_solver();
  
  // Loop over each point
  for (int i=1; i<=numpts; i++)
  {
    double x = xpts.get(i,1);
    double y = xpts.get(i,2);

    const double& rho_i = qvals.get(i,_rho_i);
    //const double& rho_e = qvals.get(i,_rho_e);
    //
    const double ion_charge_dens = rho_i*q_mi_inv;
    //const double elc_charge_dens = rho_e*q_me_inv;
    //const double charge_density = ion_charge_dens + elc_charge_dens;
    const int& meqn = source.getsize(2);
    if(meqn >= _psi)
    {
      const double psi = qvals.get(i,_psi);
      const double psi_source = -eps_Bcp_func(psi,time);
      source.set(i,_psi,  psi_source);
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
    
    const double& B1    = qvals.get(i,_B1);
    const double& B2    = qvals.get(i,_B2);
    const double& E3    = qvals.get(i,_E3);

    // solve the electro-momentum system
    //
    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE);
    {
      const double J1_i = M1_i*q_mi_inv;
      const double J2_i = M2_i*q_mi_inv;
      const double J3_i = M3_i*q_mi_inv;
      //const double J1_e = M1_e*q_me_inv;
      //const double J2_e = M2_e*q_me_inv;
      //const double J3_e = M3_e*q_me_inv;
      //const double J1 = J1_i + J1_e;
      //const double J2 = J2_i + J2_e;
      const double J3 = 2*J3_i; // + J3_e;
      source.fetch(i, _M1_i) -= J3_i*B2;
      source.fetch(i, _M2_i) += J3_i*B1;
      source.fetch(i, _M3_i) += ion_charge_dens*E3 + (J1_i*B2 - J2_i*B1);
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

    // damp deviatoric pressure and heat flux
    //
    // Collisions are better handled by being split off and solved exactly
    if(plasmaParams.get_source_type_isoP() == SourceType::SIMPLE)
    {
      assert_ge(plasmaParams.get_ion_base_iso_period(), 0.);
      if(plasmaParams.get_ion_base_iso_rate() > 0.0)
      {
        const double  p_i   = (P11_i + P22_i + P33_i)*(1./3.);
        //
        double ion_iso_rate = Relax::get_iso_rate(
          plasmaParams.get_iso_period_type(),
          plasmaParams.get_ion_base_iso_rate(),
          plasmaParams.get_ion_mass_inv(), rho_i,
          p_i, P11_i, P12_i, P13_i, P22_i, P23_i, P33_i);
        // the isotropization period should not exceed the time step
        // (or the user should split this off and solve exactly)
        assert_ge(dt_inv, ion_iso_rate);
        const double minus_ion_iso_rate = -ion_iso_rate;
        const double CP11_i = (P11_i-p_i)*minus_ion_iso_rate;
        const double CP12_i = (P12_i    )*minus_ion_iso_rate;
        const double CP13_i = (P13_i    )*minus_ion_iso_rate;
        const double CP22_i = (P22_i-p_i)*minus_ion_iso_rate;
        const double CP23_i = (P23_i    )*minus_ion_iso_rate;
        const double CP33_i = (P33_i-p_i)*minus_ion_iso_rate;
        source.fetch(i, _N11_i) += CP11_i;
        source.fetch(i, _N12_i) += CP12_i;
        source.fetch(i, _N13_i) += CP13_i;
        source.fetch(i, _N22_i) += CP22_i;
        source.fetch(i, _N23_i) += CP23_i;
        source.fetch(i, _N33_i) += CP33_i;

        // relaxing heat flux unfortunately requires computing
        // heat flux from conserved variables.
        //
        const double Sym3uN111_i = u1_i*N11_i*3;
        const double Sym3uN112_i = u1_i*N12_i*2 + u2_i*N11_i;
        const double Sym3uN113_i = u1_i*N13_i*2 + u3_i*N11_i;
        const double Sym3uN122_i = u1_i*N22_i + 2*u2_i*N12_i;
        const double Sym3uN123_i = u1_i*N23_i + u2_i*N13_i + u3_i*N12_i;
        const double Sym3uN133_i = u1_i*N33_i + 2*u3_i*N13_i;
        const double Sym3uN222_i = u2_i*N22_i*3;
        const double Sym3uN223_i = u2_i*N23_i*2 + u3_i*N22_i;
        const double Sym3uN233_i = u2_i*N33_i + 2*u3_i*N23_i;
        const double Sym3uN333_i = u3_i*N33_i*3;
        //
        // Q = QN - Sym3uN + 2*rho*uuu;
        const double Q111_i = QN111_i - Sym3uN111_i + 2*u1_i*K11_i;
        const double Q112_i = QN112_i - Sym3uN112_i + 2*u1_i*K12_i;
        const double Q113_i = QN113_i - Sym3uN113_i + 2*u1_i*K13_i;
        const double Q122_i = QN122_i - Sym3uN122_i + 2*u1_i*K22_i;
        const double Q123_i = QN123_i - Sym3uN123_i + 2*u1_i*K23_i;
        const double Q133_i = QN133_i - Sym3uN133_i + 2*u1_i*K33_i;
        const double Q222_i = QN222_i - Sym3uN222_i + 2*u2_i*K22_i;
        const double Q223_i = QN223_i - Sym3uN223_i + 2*u2_i*K23_i;
        const double Q233_i = QN233_i - Sym3uN233_i + 2*u2_i*K33_i;
        const double Q333_i = QN333_i - Sym3uN333_i + 2*u3_i*K33_i;
        //
        // CQ = -Q/tau
        const double CQ111_i = Q111_i*minus_ion_iso_rate;
        const double CQ112_i = Q112_i*minus_ion_iso_rate;
        const double CQ113_i = Q113_i*minus_ion_iso_rate;
        const double CQ122_i = Q122_i*minus_ion_iso_rate;
        const double CQ123_i = Q123_i*minus_ion_iso_rate;
        const double CQ133_i = Q133_i*minus_ion_iso_rate;
        const double CQ222_i = Q222_i*minus_ion_iso_rate;
        const double CQ223_i = Q223_i*minus_ion_iso_rate;
        const double CQ233_i = Q233_i*minus_ion_iso_rate;
        const double CQ333_i = Q333_i*minus_ion_iso_rate;

        // Sym3uCPijk = ui*CPjk+uj*CPik+uk*CPij;
        const double Sym3uCP111_i = 3*u1_i*CP11_i;
        const double Sym3uCP112_i = 2*u1_i*CP12_i+u2_i*CP11_i;
        const double Sym3uCP113_i = 2*u1_i*CP13_i+u3_i*CP11_i;
        const double Sym3uCP122_i = u1_i*CP22_i + 2*u2_i*CP12_i;
        const double Sym3uCP123_i = u1_i*CP23_i + u2_i*CP13_i + u3_i*CP12_i;
        const double Sym3uCP133_i = u1_i*CP33_i + 2*u3_i*CP13_i;
        const double Sym3uCP222_i = 3*u2_i*CP22_i;
        const double Sym3uCP223_i = 2*u2_i*CP23_i+u3_i*CP22_i;
        const double Sym3uCP233_i = u2_i*CP33_i + 2*u3_i*CP23_i;
        const double Sym3uCP333_i = 3*u3_i*CP33_i;

        // CQN = CQ + Sym3uCP
        source.fetch(i, _Q111_i) += CQ111_i + Sym3uCP111_i;
        source.fetch(i, _Q112_i) += CQ112_i + Sym3uCP112_i;
        source.fetch(i, _Q113_i) += CQ113_i + Sym3uCP113_i;
        source.fetch(i, _Q122_i) += CQ122_i + Sym3uCP122_i;
        source.fetch(i, _Q123_i) += CQ123_i + Sym3uCP123_i;
        source.fetch(i, _Q133_i) += CQ133_i + Sym3uCP133_i;
        source.fetch(i, _Q222_i) += CQ222_i + Sym3uCP222_i;
        source.fetch(i, _Q223_i) += CQ223_i + Sym3uCP223_i;
        source.fetch(i, _Q233_i) += CQ233_i + Sym3uCP233_i;
        source.fetch(i, _Q333_i) += CQ333_i + Sym3uCP333_i;
      }
    }

    // non-collisional evolution of energy tensor
    //
    // (Evolve kinetic energy tensor and rotate the pressure tensor.)
    //
    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE);
    {
      // (q/m)*Sym2(E*M); recall that E1=E2=0.
      //
      //source.fetch(i,_N11_i)+=q_mi_inv*2.0*(M1_i*E1);
      //source.fetch(i,_N12_i)+=q_mi_inv*(M1_i*E2+M2_i*E1);
      source.fetch(i,_N13_i)+=q_mi_inv*(M1_i*E3);
      //source.fetch(i,_N22_i)+=q_mi_inv*2.0*(M2_i*E2);
      source.fetch(i,_N23_i)+=q_mi_inv*(M2_i*E3);
      source.fetch(i,_N33_i)+=q_mi_inv*2.0*(M3_i*E3);
    }

    assert(plasmaParams.get_source_type_Emom() == SourceType::SIMPLE
    && plasmaParams.get_source_type_rotP() == SourceType::SIMPLE);
    //
    // (q/m)Sym2(NxB)
    // (q/m)(NxB)_li = (q/m)N
    source.fetch(i,_N11_i)+=q_mi_inv*2.0*(-N13_i*B2);
    source.fetch(i,_N12_i)+=q_mi_inv*(N13_i*B1-N23_i*B2);
    source.fetch(i,_N13_i)+=q_mi_inv*(N11_i*B2-N12_i*B1-N33_i*B2);
    source.fetch(i,_N22_i)+=q_mi_inv*2.0*(N23_i*B1);
    source.fetch(i,_N23_i)+=q_mi_inv*(N12_i*B2-N22_i*B1+N33_i*B1);
    source.fetch(i,_N33_i)+=q_mi_inv*2.0*(N13_i*B2-N23_i*B1);
    //
    //source.fetch(i,_N11_i)+=q_mi_inv*2.0*(M1_i*E1+N12_i*B3-N13_i*B2);
    //source.fetch(i,_N12_i)+=q_mi_inv*(M1_i*E2+M2_i*E1+N13_i*B1-N11_i*B3+N22_i*B3-N23_i*B2);
    //source.fetch(i,_N13_i)+=q_mi_inv*(M1_i*E3+M3_i*E1+N11_i*B2-N12_i*B1+N23_i*B3-N33_i*B2);
    //source.fetch(i,_N22_i)+=q_mi_inv*2.0*(M2_i*E2+N23_i*B1-N12_i*B3);
    //source.fetch(i,_N23_i)+=q_mi_inv*(M2_i*E3+M3_i*E2+N12_i*B2-N22_i*B1+N33_i*B1-N13_i*B3);
    //source.fetch(i,_N33_i)+=q_mi_inv*2.0*(M3_i*E3+N13_i*B2-N23_i*B1);
    //  QxB_ij1= Q_ij2*B3-Q_ij3*B2;
    //  QxB_ij2= Q_ij3*B1-Q_ij1*B3;
    //  QxB_ij3= Q_ij1*B2-Q_ij2*B1;

    // Sym3ENijk = Ei*Njk+Ej*Nik+Ek*Nij;
    //
    const double Sym3EN_111_i = 0.;
    const double Sym3EN_112_i = 0.;
    const double Sym3EN_113_i = E3*N11_i;
    const double Sym3EN_122_i = 0.;
    const double Sym3EN_123_i = E3*N12_i;
    const double Sym3EN_133_i = 2*E3*N13_i;
    const double Sym3EN_222_i = 0.;
    const double Sym3EN_223_i = E3*N22_i;
    const double Sym3EN_233_i = 2*E3*N23_i;
    const double Sym3EN_333_i = 3*E3*N33_i;
    //
    // Sym3QNxB_ijk
    //
    const double Sym3QNxB_111_i = -3*QN113_i*B2;
    const double Sym3QNxB_112_i = QN113_i*B1 - 2*QN123_i*B2;
    const double Sym3QNxB_113_i = QN111_i*B2 - QN112_i*B1 - 2*QN133_i*B2;
    const double Sym3QNxB_122_i = 2*QN123_i*B1 - QN223_i*B2;
    const double Sym3QNxB_123_i = QN112_i*B2 - QN122_i*B1 + QN133_i*B1 - QN233_i*B2;
    const double Sym3QNxB_133_i = 2*(QN113_i*B2 - QN123_i*B1) - QN333_i*B2;
    const double Sym3QNxB_222_i = 3*QN223_i*B1;
    const double Sym3QNxB_223_i = QN122_i*B2 - QN222_i*B1 + 2*QN233_i*B1;
    const double Sym3QNxB_233_i = 2*(QN123_i*B2 - QN223_i*B1) + QN333_i*B1;
    const double Sym3QNxB_333_i = 3*(QN133_i*B2 - QN233_i*B1);

    source.set(i,_Q111_i, q_mi_inv*(               Sym3QNxB_111_i));
    source.set(i,_Q112_i, q_mi_inv*(               Sym3QNxB_112_i));
    source.set(i,_Q113_i, q_mi_inv*(Sym3EN_113_i + Sym3QNxB_113_i));
    source.set(i,_Q122_i, q_mi_inv*(               Sym3QNxB_122_i));
    source.set(i,_Q123_i, q_mi_inv*(Sym3EN_123_i + Sym3QNxB_123_i));
    source.set(i,_Q133_i, q_mi_inv*(Sym3EN_133_i + Sym3QNxB_133_i));
    source.set(i,_Q222_i, q_mi_inv*(               Sym3QNxB_222_i));
    source.set(i,_Q223_i, q_mi_inv*(Sym3EN_223_i + Sym3QNxB_223_i));
    source.set(i,_Q233_i, q_mi_inv*(Sym3EN_233_i + Sym3QNxB_233_i));
    source.set(i,_Q333_i, q_mi_inv*(Sym3EN_333_i + Sym3QNxB_333_i));
  }
}
