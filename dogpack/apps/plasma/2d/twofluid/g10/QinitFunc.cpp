#include <cmath>
#include "dogdefs.h"
#include "FiveMomentParams.h"
#include "DogParamsCart2.h"
#include "PlasmaParams.h"
#include "GEMparams.h"
#include "GEM.h"
#include "gas05.h"
#include "Components.h"
#include "Params.h"
#include "constants.h" // for pi

#include <stdlib.h> // for rand, srand
#include <time.h> // to seed

static void seed_rand()
{
  srand((unsigned)time(NULL)); 
}

static double rand_double()
{
  // produces a number in [0, 1)
  double scale_inv = 1./(RAND_MAX+1.);
  //out = double(rand())/(RAND_MAX+1.);
  int outi = rand();
  //dprint(outi)
  //dprint(scale_inv)
  double out = double(outi)*scale_inv;
  //dprint(out)
  return out;
}

// return a double in the interval [inc, exc);
// may have exc < inc.
static double rand_double(double inc, double exc)
{
  double diff = exc-inc;
  double shift = rand_double();
  double retval = inc+shift*diff;
  return retval;
}

// create a random positive definite pressure tensor
static void get_random_pressure(
  double& P11,
  double& P12,
  double& P13,
  double& P22,
  double& P23,
  double& P33)
{
  // create a random matrix
  double r[4][4];
  for(int i=1;i<=3;i++)
  for(int j=1;j<=3;j++)
  {
    r[i][j] = rand_double(-1.,1);
  }
  // Let P = transpose(r).r (so that P is positive definite)
  P11 = 0; P12 = 0; P13 = 0; P22 = 0; P23 = 0; P33 = 0;
  for(int k=1;k<=3;k++)
  {
    P11 += r[k][1]*r[k][1];
    P12 += r[k][1]*r[k][2];
    P13 += r[k][1]*r[k][3];
    P22 += r[k][2]*r[k][2];
    P23 += r[k][2]*r[k][3];
    P33 += r[k][3]*r[k][3];
  }
  // confirm that pressure tensor is positive definite
  const double Adj33 = P11*P22 - P12*P12;
  const double Adj13 = P12*P23 - P22*P13;
  const double Adj23 = P12*P13 - P11*P23;
  const double detP = Adj13*P13 + Adj23*P23 + Adj33*P33;
  assert_ge(P11,0.);
  assert_ge(Adj33,0.);
  assert_ge(detP,0.);
  const double p = (1./3.)*(P11+P22+P33);
  dprint(p)
  dprint(P11);
  dprint(Adj33);
  dprint(detP);
}

void QinitODE(dTensor2& qvals)
{
  const int meqn = qvals.getsize(2);

  // Brio-Wu parameters
  //
  const double& ion_mass = plasmaParams.get_ion_mass();
  const double& elc_mass = plasmaParams.get_elc_mass();

  //
  // initial values
  //
  const double ni = 2.;//+rand_double(-1.,1);
  const double ne = 2.;//+rand_double(-1.,1);
  const double rho_i   = ni*ion_mass;
  const double u1_i    = 0.; // rand_double(-1.,1); //
  const double u2_i    = 0.; // rand_double(-1.,1); //
  const double u3_i    = 0.; // rand_double(-1.,1); //
  const double press_i = 0.5;
        //
  const double rho_e   = ne*elc_mass;
  const double u1_e    = 0.; // rand_double(-1.,1); //
  const double u2_e    = 0.; // rand_double(-1.,1); //
  const double u3_e    = 0.; // rand_double(-1.,1); //
  const double press_e = 0.5;
        //
  const double B1 = 0.; // rand_double(-1.,1); //
  const double B2 = 0.; // rand_double(-1.,1); //
  const double B3 = 1.; // rand_double(-1.,1); //
  const double E1 = 1.; // rand_double(-1.,1); //
  const double E2 = 0.; // rand_double(-1.,1); //
  const double E3 = 1.; // rand_double(-1.,1); //

  //
  // derived values
  //
  double const M1_i = rho_i*u1_i;
  double const M2_i = rho_i*u2_i;
  double const M3_i = rho_i*u3_i;
  //
  double const M1_e = rho_e*u1_e;
  double const M2_e = rho_e*u2_e;
  double const M3_e = rho_e*u3_e;
  //
  double P11_i, P11_e,
         P12_i, P12_e,
         P13_i, P13_e,
         P22_i, P22_e,
         P23_i, P23_e,
         P33_i, P33_e;

  P11_i = press_i+ press_i*(2./3.);
  P12_i = 0.;
  P13_i = 0.;
  P22_i = press_i- press_i*(1./3.);
  P23_i = 0.;
  P33_i = press_i- press_i*(1./3.);

  P11_e = press_e+ press_e*(2./3.);
  P12_e = 0.;
  P13_e = 0.;
  P22_e = press_e- press_e*(1./3.);
  P23_e = 0.;
  P33_e = press_e- press_e*(1./3.);

  //get_random_pressure(P11_i, P12_i, P13_i, P22_i, P23_i, P33_i);
  //get_random_pressure(P11_e, P12_e, P13_e, P22_e, P23_e, P33_e);
  //
  const double N11_i = M1_i*u1_i + P11_i;
  const double N12_i = M1_i*u2_i + P12_i;
  const double N13_i = M1_i*u3_i + P13_i;
  const double N22_i = M2_i*u2_i + P22_i;
  const double N23_i = M2_i*u3_i + P23_i;
  const double N33_i = M3_i*u3_i + P33_i;
  //
  const double N11_e = M1_e*u1_e + P11_e;
  const double N12_e = M1_e*u2_e + P12_e;
  const double N13_e = M1_e*u3_e + P13_e;
  const double N22_e = M2_e*u2_e + P22_e;
  const double N23_e = M2_e*u3_e + P23_e;
  const double N33_e = M3_e*u3_e + P33_e;
  //
  if(debug2){
    const double one_over_epsilon = maxwellParams.get_cs_light_squared();
    const double Oi2 = (ni/ion_mass)*one_over_epsilon;
    const double Oe2 = (ne/elc_mass)*one_over_epsilon;
    const double Op2 = Oi2+Oe2;
    const double Oi = sqrt(Oi2);
    const double Oe = sqrt(Oe2);
    const double Op = sqrt(Op2);
    println("plasma variables:")
    printvc(ni); printvc(ne); printf("\n");
    printvc(rho_i); printvc(rho_e); printf("\n");
    println("frequencies and periods:");
    printvc(Oi); printvn(2*pi/Oi)
    printvc(Oe); printvn(2*pi/Oe)
    printvc(Op); printvn(2*pi/Op)
  }
  if(debug2){
    println("initial conditions:");
    printvn(rho_i);
    printvn(M1_i );
    printvn(M2_i );
    printvn(M3_i );
    printvn(N11_i);
    printvn(N12_i);
    printvn(N13_i);
    printvn(N22_i);
    printvn(N23_i);
    printvn(N33_i);

    printvn(rho_e);
    printvn(M1_e );
    printvn(M2_e );
    printvn(M3_e );
    printvn(N11_e);
    printvn(N12_e);
    printvn(N13_e);
    printvn(N22_e);
    printvn(N23_e);
    printvn(N33_e);

    printvn(B1);
    printvn(B2);
    printvn(B3);
    printvn(E1);
    printvn(E2);
    printvn(E3);
  }
  // compute the eigenvalues
  {
  }
  // Initial conditions
  for (int i=1; i<=qvals.getsize(1); i++)
  {
    qvals.set(i,_rho_i, rho_i );
    qvals.set(i,_M1_i , M1_i );
    qvals.set(i,_M2_i , M2_i );
    qvals.set(i,_M3_i , M3_i );
    qvals.set(i,_N11_i, N11_i);
    qvals.set(i,_N12_i, N12_i);
    qvals.set(i,_N13_i, N13_i);
    qvals.set(i,_N22_i, N22_i);
    qvals.set(i,_N23_i, N23_i);
    qvals.set(i,_N33_i, N33_i);

    qvals.set(i,_rho_e, rho_e );
    qvals.set(i,_M1_e , M1_e );
    qvals.set(i,_M2_e , M2_e );
    qvals.set(i,_M3_e , M3_e );
    qvals.set(i,_N11_e, N11_e);
    qvals.set(i,_N12_e, N12_e);
    qvals.set(i,_N13_e, N13_e);
    qvals.set(i,_N22_e, N22_e);
    qvals.set(i,_N23_e, N23_e);
    qvals.set(i,_N33_e, N33_e);

    qvals.set(i,_B1 , B1 );
    qvals.set(i,_B2 , B2 );
    qvals.set(i,_B3 , B3 );
    qvals.set(i,_E1 , E1 );
    qvals.set(i,_E2 , E2 );
    qvals.set(i,_E3 , E3 );
    if(_psi > meqn) continue;
    qvals.set(i,_psi  , 0. );
    if(_phi > meqn) continue;
    qvals.set(i,_phi  , 0. );
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
    const double& left_P11_i = left_press_i;
    const double  left_P12_i = 0.0;
    const double  left_P13_i = 0.0;
    const double& left_P22_i = left_press_i;
    const double  left_P23_i = 0.0;
    const double& left_P33_i = left_press_i;
    const double  left_entropy_i=get_entropy_per_vol05(left_rho_i,left_press_i);
    //
    const double& left_P11_e = left_press_e;
    const double  left_P12_e = 0.0;
    const double  left_P13_e = 0.0;
    const double& left_P22_e = left_press_e;
    const double  left_P23_e = 0.0;
    const double& left_P33_e = left_press_e;
    const double  left_entropy_e=get_entropy_per_vol05(left_rho_e,left_press_e);

    double const rght_M1_i = rght_rho_i*rght_u1_i;
    double const rght_M2_i = rght_rho_i*rght_u2_i;
    double const rght_M3_i = rght_rho_i*rght_u3_i;
    //
    double const rght_M1_e = rght_rho_e*rght_u1_e;
    double const rght_M2_e = rght_rho_e*rght_u2_e;
    double const rght_M3_e = rght_rho_e*rght_u3_e;
    //
    const double& rght_P11_i = rght_press_i;
    const double  rght_P12_i = 0.0;
    const double  rght_P13_i = 0.0;
    const double& rght_P22_i = rght_press_i;
    const double  rght_P23_i = 0.0;
    const double& rght_P33_i = rght_press_i;
    const double  rght_entropy_i=get_entropy_per_vol05(rght_rho_i,rght_press_i);
    //
    const double& rght_P11_e = rght_press_e;
    const double  rght_P12_e = 0.0;
    const double  rght_P13_e = 0.0;
    const double& rght_P22_e = rght_press_e;
    const double  rght_P23_e = 0.0;
    const double& rght_P33_e = rght_press_e;
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
            qvals.set(i,_rho_i, left_rho_i );
            qvals.set(i,_M1_i , left_M1_i );
            qvals.set(i,_M2_i , left_M2_i );
            qvals.set(i,_M3_i , left_M3_i );
            qvals.set(i,_N11_i, left_M1_i*left_u1_i + left_P11_i);
            qvals.set(i,_N12_i, left_M1_i*left_u2_i + left_P12_i);
            qvals.set(i,_N13_i, left_M1_i*left_u3_i + left_P13_i);
            qvals.set(i,_N22_i, left_M2_i*left_u2_i + left_P22_i);
            qvals.set(i,_N23_i, left_M2_i*left_u3_i + left_P23_i);
            qvals.set(i,_N33_i, left_M3_i*left_u3_i + left_P33_i);

            qvals.set(i,_rho_e, left_rho_e );
            qvals.set(i,_M1_e , left_M1_e );
            qvals.set(i,_M2_e , left_M2_e );
            qvals.set(i,_M3_e , left_M3_e );
            qvals.set(i,_N11_e, left_M1_e*left_u1_e + left_P11_e);
            qvals.set(i,_N12_e, left_M1_e*left_u2_e + left_P12_e);
            qvals.set(i,_N13_e, left_M1_e*left_u3_e + left_P13_e);
            qvals.set(i,_N22_e, left_M2_e*left_u2_e + left_P22_e);
            qvals.set(i,_N23_e, left_M2_e*left_u3_e + left_P23_e);
            qvals.set(i,_N33_e, left_M3_e*left_u3_e + left_P33_e);

            qvals.set(i,_B1 , left_B1 );
            qvals.set(i,_B2 , left_B2 );
            qvals.set(i,_B3 , left_B3 );
            qvals.set(i,_E1 , left_E1 );
            qvals.set(i,_E2 , left_E2 );
            qvals.set(i,_E3 , left_E3 );
            qvals.set(i,_psi, 0.0 );
            qvals.set(i,_phi, 0.0 );
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
            qvals.set(i,_N11_i,  rght_M1_i*rght_u1_i + rght_P11_i);
            qvals.set(i,_N12_i,  rght_M1_i*rght_u2_i + rght_P12_i);
            qvals.set(i,_N13_i,  rght_M1_i*rght_u3_i + rght_P13_i);
            qvals.set(i,_N22_i,  rght_M2_i*rght_u2_i + rght_P22_i);
            qvals.set(i,_N23_i,  rght_M2_i*rght_u3_i + rght_P23_i);
            qvals.set(i,_N33_i,  rght_M3_i*rght_u3_i + rght_P33_i);

            qvals.set(i,_rho_e, rght_rho_e );
            qvals.set(i,_M1_e , rght_M1_e );
            qvals.set(i,_M2_e , rght_M2_e );
            qvals.set(i,_M3_e , rght_M3_e );
            qvals.set(i,_N11_e, rght_M1_e*rght_u1_e + rght_P11_e);
            qvals.set(i,_N12_e, rght_M1_e*rght_u2_e + rght_P12_e);
            qvals.set(i,_N13_e, rght_M1_e*rght_u3_e + rght_P13_e);
            qvals.set(i,_N22_e, rght_M2_e*rght_u2_e + rght_P22_e);
            qvals.set(i,_N23_e, rght_M2_e*rght_u3_e + rght_P23_e);
            qvals.set(i,_N33_e, rght_M3_e*rght_u3_e + rght_P33_e);

            qvals.set(i,_B1 , rght_B1 );
            qvals.set(i,_B2 , rght_B2 );
            qvals.set(i,_B3 , rght_B3 );
            qvals.set(i,_E1 , rght_E1 );
            qvals.set(i,_E2 , rght_E2 );
            qvals.set(i,_E3 , rght_E3 );
            qvals.set(i,_psi, 0.0 );
            qvals.set(i,_phi, 0.0 );
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

        qvals.set(i,_rho_e, rho_e);
        qvals.set(i,_M1_e , 0.0  );
        qvals.set(i,_M2_e , 0.0  );
        qvals.set(i,_M3_e , M3_e );
        qvals.set(i,_N11_e, press_e);
        qvals.set(i,_N12_e, 0.);
        qvals.set(i,_N13_e, 0.);
        qvals.set(i,_N22_e, press_e);
        qvals.set(i,_N23_e, 0.);
        qvals.set(i,_N33_e, M3_e*M3_e/rho_e + press_e );

        qvals.set(i,_B1   , B1 );
        qvals.set(i,_B2   , B2 );
        qvals.set(i,_B3   , B_guide );
        qvals.set(i,_E1   , 0. );
        qvals.set(i,_E2   , 0. );
        qvals.set(i,_E3   , 0. );

        if(_psi > meqn) continue;
        qvals.set(i,_psi  , 0. );
        if(_phi > meqn) continue;
        qvals.set(i,_phi  , 0. );

        if(_entropy_i > meqn) continue;
        qvals.set(i,_entropy_i, get_entropy_per_vol05(rho_i,press_i));
        if(_entropy_e > meqn) continue;
        qvals.set(i,_entropy_e, get_entropy_per_vol05(rho_e,press_e));
    }
}

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  using namespace Params;
  switch(get_problem())
  {
    case Params::GEM:
    case Params::STEADY:
      QinitGEM(xpts, qvals);
      break;
    case Params::SHOCK:
      QinitShock(xpts,  qvals);
      break;
    case Params::ODE:
      QinitODE(qvals);
      break;
    default:
      invalid_value_error(Params::get_problem());
  }
}

