#include <cmath>
#include "assert.h"
#include "Source2Fluid.h"
#include "Polynomial.h"
#include "gas20.h"
#include "Relax.h"
#include "constants.h" // for pi

void set_B_orth_system_sympair(double a[4][4],
  double B1,double B2,double Bmag)
{
  const double Bmag_inv = 1./Bmag;
  const double b1 = B1*Bmag_inv;
  const double b2 = B2*Bmag_inv;
  const double b1sq = b1*b1;
  const double b2sq = b2*b2;
  //double a1[4] = a[1];
  //double a2[4] = a[2];
  // These are row vectors. e.g. a[3][2] = a3[2]
  double * a1 = a[1];
  double *  e = a[2];
  double * a3 = a[3];
  a3[1] = b1;  // a[3][1]
  a3[2] = b2;  // a[3][2]
  a3[3] = 0.;  // a[3][3]
  // rotate in the plane by 90 degrees
  a1[1] = -b2; // a[1][1]
  a1[2] = b1;  // a[1][2]
  a1[3] = 0.;  // a[1][3]
  // third vector (a2) is out-of-plane
  e[1] = 0.;   // a[2][1]
  e[2] = 0.;   // a[2][2]
  e[3] = 1.;   // a[2][3]
  // check that a is an orthonormal system
  if(false)
  {
    for(int k1=1;k1<=3;k1++)
    {
      for(int k2=1;k2<=k1;k2++)
      {
        // compute a[k1].a[k2]
        double acc=0.;
        for(int j=1;j<=3;j++)
        {
          acc += a[k1][j]*a[k2][j];
        }
        if(k1==k2)
          assert_almost_eq(acc, 1.);
        else
          assert_almost_eq(acc, 0.);
      }
    }
  }
}

void set_B_orth_system(double a[4][4],
  double B1,double B2,double B3,double Bmag)
{
  const double Bmag_inv = 1./Bmag;
  const double b1 = B1*Bmag_inv;
  const double b2 = B2*Bmag_inv;
  const double b3 = B3*Bmag_inv;
  const double b1sq = b1*b1;
  const double b2sq = b2*b2;
  const double b3sq = b3*b3;
  //double a1[4] = a[1];
  //double a2[4] = a[2];
  // These are row vectors. e.g. a[3][2] = a3[2]
  double * a1 = a[1];
  double *  e = a[2];
  double * a3 = a[3];
  a3[1] = b1;
  a3[2] = b2;
  a3[3] = b3;
  // find an elementary basis vector ej for which
  // ej.b is small and use it to construct the desired basis.
  if(b1sq<.34) // e1
  {
    double norm_sq = b2sq+b3sq;
    double norm = sqrt(norm_sq);
    double norm_inv = 1./norm;
    a1[1] = 0.;
    a1[2] = -b3*norm_inv;
    a1[3] =  b2*norm_inv;
    e[1] = norm;
    e[2] = -b1*a1[3];
    e[3] =  b1*a1[2];
  }
  else if(b2sq<.34) // e2
  {
    double norm_sq = b3sq+b1sq;
    double norm = sqrt(norm_sq);
    double norm_inv = 1./sqrt(norm_sq);
    a1[2] = 0.;
    a1[3] = -b1*norm_inv;
    a1[1] =  b3*norm_inv;
    e[2] = norm;
    e[3] = -b2*a1[1];
    e[1] =  b2*a1[3];
  }
  else
  {
    assert_lt(b3sq,.34); // e3
    double norm_sq = b1sq+b2sq;
    double norm = sqrt(norm_sq);
    double norm_inv = 1./sqrt(norm_sq);
    a1[3] = 0.;
    a1[1] = -b2*norm_inv;
    a1[2] =  b1*norm_inv;
    e[3] = norm;
    e[1] = -b3*a1[2];
    e[2] =  b3*a1[1];
  }
  // check that a is an orthonormal system
  if(false)
  {
    for(int k1=1;k1<=3;k1++)
    {
      for(int k2=1;k2<=k1;k2++)
      {
        // compute a[k1].a[k2]
        double acc=0.;
        for(int j=1;j<=3;j++)
        {
          acc += a[k1][j]*a[k2][j];
        }
        if(k1==k2)
          assert_almost_eq(acc, 1.);
        else
          assert_almost_eq(acc, 0.);
      }
    }
  }
}

// in the symmetric pair plasma case
// we need only solve the perpendicular system
void SourceSolver::solve_electro_momentum_sympair(
  dTensor2& prim,
  const double d_time,
  int _rho_i,
  int _B1, int _B2, int _E3)
{
  const double mi_inv = plasmaParams.get_ion_mass_inv();
  const double me_inv = plasmaParams.get_elc_mass_inv();
  assert_eq(mi_inv,me_inv);
  const double one_over_eps = maxwellParams.get_one_over_epsilon();
  const double Oi2_factor = one_over_eps*mi_inv;

  for(int i=1;i<=prim.getsize(1);i++)
  {
    // get the magnetic field
    const double B1 = prim.get(i,_B1);
    const double B2 = prim.get(i,_B2);
    const double B1sq = B1*B1;
    const double B2sq = B2*B2;
    const double Bmag2 = B1sq + B2sq;
    const double Bmag = sqrt(Bmag2);
    // We would need to handle this degenerate case separately. For
    // this degenerate case we can simply solve the parallel system
    // along each principal axis. For my discretization of the
    // GEM problem the degenerate case should never hold at a
    // point where the source term is evaluated.
    assert_ne(Bmag,0.);

    // get the densities
    //
    const double rho_i = prim.get(i, _rho_i);
    const double ni = rho_i*mi_inv;

    // compute the plasma frequencies
    //
    const double Oi2 = ni*Oi2_factor;
    const double Op2 = 2*Oi2;
    const double Op = sqrt(Op2);
    if(debug3)
    {
      const double Oi = sqrt(Oi2);
      println("plasma frequencies and periods:")
      printvc(Oi); printvn(2*pi/Oi);
      printvc(Op); printvn(2*pi/Op);
    }
    //
    // compute the gyrofrequences
    //
    const double gi = Bmag*mi_inv;
    if(debug3)
    {
      println("gyrofrequencies and gyroperiods:")
      printvc(gi); printvn(2*pi/gi)
    }
    //
    // compute the eigenfrequencies
    //
    // eigenfrequencies om must satisfy
    // 0 = om^3 - (gi^2+Op2)*om, so
    // om = +/-sqrt(gi^2+Op2) or om=0.
    double om[3];
    double gi2 = gi*gi;
    om[2]=sqrt(gi2+Op2);
    om[1]=0.;
    om[0]=-om[2];
    if(debug3)
    {
      println("frequencies and periods:")
      printvc(om[0]); printvn(2*pi/om[0]);
      printvc(om[1]); printvn(2*pi/om[1]);
      printvc(om[2]); printvn(2*pi/om[2]);
    }

    // determine a positively oriented orthonormal system
    // exb, e, b, that is,
    // a1, a2, a3 (where a1 cross a2 = a3)
    //
    double a[4][4];
    set_B_orth_system_sympair(a,B1,B2,Bmag);

    // rescale the velocities as needed to have an antisymmetric system.
    //
    const double ui0_inv = sqrt(rho_i*one_over_eps);
    double E[4],ui[4];//,ue[4]
    ui[1] = prim.get(i,_rho_i+1)*ui0_inv;
    ui[2] = prim.get(i,_rho_i+2)*ui0_inv;
    ui[3] = prim.get(i,_rho_i+3)*ui0_inv;
    E[1] = 0.;
    E[2] = 0.;
    E[3] = prim.get(i,_E3);

    // project components of E, ui, and ue onto the orthonormal system
    double ui_a[4], ue_a[4], E_a[4];
    ui_a[1] = a[1][1]*ui[1] + a[1][2]*ui[2];
    ui_a[2] = ui[3];
    ui_a[3] = a[3][1]*ui[1] + a[3][2]*ui[2];
    ue_a[1] = ui_a[1];
    ue_a[2] = -ui_a[2];
    ue_a[3] = ui_a[3];
    E_a[1] = 0.;
    E_a[2] = E[3];
    E_a[3] = 0.;
    if(debug3)
    {
      double ue[4];
      ue[1] = ui[1];
      ue[2] = ui[2];
      ue[3] =-ui[3];
      assert_almost_eq(ui_a[1], a[1][1]*ui[1] + a[1][2]*ui[2] + a[1][3]*ui[3]);
      assert_almost_eq(ui_a[2], a[2][1]*ui[1] + a[2][2]*ui[2] + a[2][3]*ui[3]);
      assert_almost_eq(ui_a[3], a[3][1]*ui[1] + a[3][2]*ui[2] + a[3][3]*ui[3]);
      assert_almost_eq(ue_a[1], a[1][1]*ue[1] + a[1][2]*ue[2] + a[1][3]*ue[3]);
      assert_almost_eq(ue_a[2], a[2][1]*ue[1] + a[2][2]*ue[2] + a[2][3]*ue[3]);
      assert_almost_eq(ue_a[3], a[3][1]*ue[1] + a[3][2]*ue[2] + a[3][3]*ue[3]);
      assert_almost_eq( E_a[1], a[1][1]* E[1] + a[1][2]* E[2] + a[1][3]* E[3]);
      assert_almost_eq( E_a[2], a[2][1]* E[1] + a[2][2]* E[2] + a[2][3]* E[3]);
      assert_almost_eq( E_a[3], a[3][1]* E[1] + a[3][2]* E[2] + a[3][3]* E[3]);
    }

    // compute eigensolution components
    // for each eigenvalue for subspace orthogonal to b
    //
    double bibe[3], Oibe[3], Oebi[3], norm_sq_inv[3];
    const double Oi = sqrt(Oi2);
    {
      bibe[0] = -Op2;
      bibe[1] = gi2;
      bibe[2] = -Op2;
      //
      Oebi[0] = Oibe[2] = Oi*(gi-om[2]);
      Oebi[1] = Oibe[1] = Oi*gi;
      Oebi[2] = Oibe[0] = Oi*(gi+om[2]);
      //
      assert_almost_eq(Oibe[0], Oi*(gi-om[0]));
      assert_almost_eq(Oibe[1], Oi*(gi-om[1]));
      assert_almost_eq(Oibe[2], Oi*(gi-om[2]));

      // normalize the eigenvectors
      double norm_sq[3], norm_inv[3];
      double bibe2[3], Oibe2[3], Oebi2[3];
      bibe2[0] = bibe[0]*bibe[0];
      bibe2[1] = bibe[1]*bibe[1];
      bibe2[2] = bibe2[0];
      Oibe2[0] = Oibe[0]*Oibe[0];
      Oibe2[1] = Oibe[1]*Oibe[1];
      Oibe2[2] = Oibe[2]*Oibe[2];
      Oebi2[0] = Oibe2[2];
      Oebi2[1] = Oibe2[1];
      Oebi2[2] = Oibe2[0];
      norm_sq[0] = bibe[0]*bibe[0] + Oibe[0]*Oibe[0] + Oebi[0]*Oebi[0];
      norm_sq[1] = bibe[1]*bibe[1] + Oibe[1]*Oibe[1] + Oebi[1]*Oebi[1];
      norm_sq[2] = norm_sq[0];
      assert_almost_eq(norm_sq[2], bibe[2]*bibe[2] + Oibe[2]*Oibe[2] + Oebi[2]*Oebi[2]);

      norm_inv[0] = 1./sqrt(norm_sq[0]);
      norm_inv[1] = 1./sqrt(norm_sq[1]);
      norm_inv[2] = norm_inv[0];
      assert_almost_eq(norm_inv[2], 1./sqrt(norm_sq[2]));

      for(int k=0;k<3;k++)
      {
        assert_almost_eq(norm_sq[k],
            bibe[k]*bibe[k]
          + Oibe[k]*Oibe[k]
          + Oebi[k]*Oebi[k]);
        bibe[k] *= norm_inv[k];
        Oibe[k] *= norm_inv[k];
        Oebi[k] *= norm_inv[k];
        // we have renormalized, so:
        norm_sq_inv[k] = 1.;
      }
    }

    // project the initial conditions onto the initial eigensolutions
    // to determine the coefficients of the expansion
    //
    double c[3][3];
    double* c0=c[0]; // use for parallel component
    double* c1=c[1]; // 1st perpendicular component
    double* c2=c[2]; // 2nd perpendicular component
    for(int k=0;k<3;k++)
    {
      c1[k] =(bibe[k]* E_a[2]
              + (Oibe[k] + Oebi[k])*ui_a[1])*norm_sq_inv[k];
      c2[k] = ((-Oibe[k] + Oebi[k])*ui_a[2])*norm_sq_inv[k];
      assert_almost_eq(c1[k],
             (bibe[k]* E_a[2]
            + Oibe[k]*ui_a[1]
            + Oebi[k]*ue_a[1])*norm_sq_inv[k]);
      assert_almost_eq(c2[k],
             (bibe[k]* E_a[1]
            - Oibe[k]*ui_a[2]
            - Oebi[k]*ue_a[2])*norm_sq_inv[k]);
    }
    //const double Op = sqrt(Op2);
    const double Op_inv = 1./Op; // normalize
    const double Oi_Op = Oi*Op_inv;
    c0[0] = 2*Oi_Op*ui_a[3];
    c0[1] = 0.;
    c0[2] = 0.;
    {
      const double Oe_Op = Oi_Op; // const double Oe_Op = Oi*Op_inv;
      assert_almost_eq(c0[0], Oe_Op*ui_a[3] + Oi_Op*ue_a[3]);
      assert_almost_eq(c0[1], Oi_Op*ui_a[3] - Oe_Op*ue_a[3]);
      assert_almost_eq(c0[2], E_a[3]);
    }

    // evaluate the eigensolutions at time d_time
    //
    // declare solutions in basis a at time t
    double E_t_a [4], ui_t_a[4], ue_t_a[4];
    //
    // the perpendicular component
    //
    {
      double om_dt[3], sin_om_dt[3], cos_om_dt[3];
      double rotc1[3], rotc2[3];
      const double om_dt_2 = om[2]*d_time;
      sin_om_dt[2] = sin(om_dt_2);
      cos_om_dt[2] = cos(om_dt_2);
      sin_om_dt[1] = 0.;
      cos_om_dt[1] = 1.;
      sin_om_dt[0] = -sin_om_dt[2];
      cos_om_dt[0] = cos_om_dt[2];
      for(int j=0;j<3;j++)
      {
        rotc1[j] =  c1[j]*cos_om_dt[j] + c2[j]*sin_om_dt[j];
        rotc2[j] = -c1[j]*sin_om_dt[j] + c2[j]*cos_om_dt[j];
      }
      E_t_a[1] =  bibe[0]*rotc2[0] + bibe[1]*rotc2[1] + bibe[2]*rotc2[2];
      E_t_a[2] =  bibe[0]*rotc1[0] + bibe[1]*rotc1[1] + bibe[2]*rotc1[2];
      ui_t_a[1] = Oibe[0]*rotc1[0] + Oibe[1]*rotc1[1] + Oibe[2]*rotc1[2];
      ui_t_a[2]=-(Oibe[0]*rotc2[0] + Oibe[1]*rotc2[1] + Oibe[2]*rotc2[2]);
      ue_t_a[1] = ui_t_a[1];
      ue_t_a[2] =-ui_t_a[2];
      assert_almost_eq(ue_t_a[1],
        Oebi[0]*rotc1[0] + Oebi[1]*rotc1[1] + Oebi[2]*rotc1[2]);
      assert_almost_eq(ue_t_a[2],
        -(Oebi[0]*rotc2[0] + Oebi[1]*rotc2[1] + Oebi[2]*rotc2[2]));
      //
      //E_t_a[1] =  bibe[0]*rotc2[0] + bibe[1]*rotc2[1] + bibe[2]*rotc2[2];
      //E_t_a[2] =  bibe[0]*rotc1[0] + bibe[1]*rotc1[1] + bibe[2]*rotc1[2];
      //ui_t_a[1] = Oibe[0]*rotc1[0] + Oibe[1]*rotc1[1] + Oibe[2]*rotc1[2];
      //ui_t_a[2]=-(Oibe[0]*rotc2[0] + Oibe[1]*rotc2[1] + Oibe[2]*rotc2[2]);
      //ue_t_a[1] = Oebi[0]*rotc1[0] + Oebi[1]*rotc1[1] + Oebi[2]*rotc1[2];
      //ue_t_a[2]=-(Oebi[0]*rotc2[0] + Oebi[1]*rotc2[1] + Oebi[2]*rotc2[2]);
    }
    //
    // the parallel component
    //
    // const double Op_dt = Op*d_time;
    // const double sin_Op_dt = sin(Op_dt);
    // const double cos_Op_dt = cos(Op_dt);
    // const double rotc1 = 0; // c0[1]*cos_Op_dt + c0[2]*sin_Op_dt;
    // const double rotc2 = 0; //-c0[1]*sin_Op_dt + c0[2]*cos_Op_dt;
    E_t_a [3] = 0; // rotc2;
    ui_t_a[3] = 0; // Oi_Op*rotc1;
    ue_t_a[3] = 0; //-Oe_Op*rotc1;

    // project components of E_t, ui_t, and ue_t back onto the standard basis
    double E_t[4], ui_t[4];//, ue_t[4],;
    ui_t[1] = a[1][1]*ui_t_a[1] + a[3][1]*ui_t_a[3];
    ui_t[2] = a[1][2]*ui_t_a[1] + a[3][2]*ui_t_a[3];
    ui_t[3] = ui_t_a[2];
    E_t[3] = E_t_a[2];
    // written out in full
    assert_almost_eq(ui_t[1],a[1][1]*ui_t_a[1]+a[2][1]*ui_t_a[2]+a[3][1]*ui_t_a[3]);
    assert_almost_eq(ui_t[2],a[1][2]*ui_t_a[1]+a[2][2]*ui_t_a[2]+a[3][2]*ui_t_a[3]);
    assert_almost_eq(ui_t[3],a[1][3]*ui_t_a[1]+a[2][3]*ui_t_a[2]+a[3][3]*ui_t_a[3]);
    assert_almost_eq( E_t[3],a[1][3]* E_t_a[1]+a[2][3]* E_t_a[2]+a[3][3]* E_t_a[3]);

    // copy out the solution, restoring the velocity scale
    const double ui0 = 1./ui0_inv;
    const double ui_new1 = ui_t[1]*ui0;
    const double ui_new2 = ui_t[2]*ui0;
    const double ui_new3 = ui_t[3]*ui0;
    if(debug3)
    {
      println("\n ui:")
      printvn(ui_new1);
      printvn(ui_new2);
      printvn(ui_new3);
      println(" E:")
      printvn(E_t[3]);
    }

    prim.set(i,_rho_i+1, ui_new1);
    prim.set(i,_rho_i+2, ui_new2);
    prim.set(i,_rho_i+3, ui_new3);
    prim.set(i,_E3, E_t[3]);
  }
}

int gsl_poly_solve_cubic (const double* p, double *roots);
void SourceSolver::solve_electro_momentum(
  dTensor2& prim,
  const double d_time,
  int _rho_i, int _rho_e,
  int _B1, int _B2, int _B3,
  int _E1, int _E2, int _E3)
{
  const double mi_inv = plasmaParams.get_ion_mass_inv();
  const double me_inv = plasmaParams.get_elc_mass_inv();
  const double one_over_eps = maxwellParams.get_one_over_epsilon();
  const double Oi2_factor = one_over_eps*mi_inv;
  const double Oe2_factor = one_over_eps*me_inv;

  for(int i=1;i<=prim.getsize(1);i++)
  {
    // get the magnetic field
    const double B1 = prim.get(i,_B1);
    const double B2 = prim.get(i,_B2);
    const double B3 = _B3 ? prim.get(i,_B3) : 0.;
    const double B1sq = B1*B1;
    const double B2sq = B2*B2;
    const double B3sq = B3*B3;
    const double Bmag2 = B1sq + B2sq + B3sq;
    const double Bmag = sqrt(Bmag2);
    // We would need to handle this degenerate case separately. For
    // this degenerate case we can simply solve the parallel system
    // along each principal axis. For my discretization of the
    // GEM problem the degenerate case should never hold at a
    // point where the source term is evaluated.
    assert_ne(Bmag,0.);

    // get the densities
    //
    const double rho_i = prim.get(i, _rho_i);
    const double rho_e = prim.get(i, _rho_e);
    const double ni = rho_i*mi_inv;
    const double ne = rho_e*me_inv;

    // compute the plasma frequencies
    //
    const double Oi2 = ni*Oi2_factor;
    const double Oe2 = ne*Oe2_factor;
    const double Op2 = Oi2+Oe2;
    const double Op = sqrt(Op2);
    //if(debug3 && printf("plasma angles:"),printf("\n"),true)
    if(debug3)
    {
      const double Oi = sqrt(Oi2);
      const double Oe = sqrt(Oe2);
      //println("plasma angles:")
      //printvn(Oi*d_time)
      //printvn(Oe*d_time)
      println("plasma frequencies and periods:")
      printvc(Oi); printvn(2*pi/Oi);
      printvc(Oe); printvn(2*pi/Oe);
      printvc(Op); printvn(2*pi/Op);
    }
    //
    // compute the gyrofrequences
    //
    const double gi = Bmag*mi_inv;
    const double ge = Bmag*me_inv;
    if(debug3)
    {
      //println("gyro angles:")
      //printvn(gi*d_time)
      //printvn(ge*d_time)
      println("gyrofrequencies and gyroperiods:")
      printvc(gi); printvn(2*pi/gi)
      printvc(ge); printvn(2*pi/ge)
    }
    //
    // compute the eigenfrequencies
    //
    // eigenfrequencies om must satisfy
    // 0 = om^3 + (gi-ge)*om^2 - (gi*ge+Oi2+Oe2)*om + (Oi2*ge-Oe2*gi)
    double om[3];
    double coef[4];
    coef[0] = Oi2*ge-Oe2*gi;
    coef[1] = -(gi*ge+Op2); // note that Op2 = Oi2+Oe2
    coef[2] = (gi-ge);
    coef[3] = 1; // this is actually assumed by the root-finder.
    if(debug3)
    {
      // compute bounds for eigenvalues
      // (based on wikipedia page on properties of polynomial roots)
      const double radicand = coef[2]*coef[2]-3*coef[1];
      assert_ge(radicand,0.);
      const double radical = sqrt(radicand);
      // const double omega_m = (1./3.)*(-coef[2] - 2*radical);
      // const double omega_p = (1./3.)*(-coef[2] + 2*radical);
      // println("frequency interval: (%24.16e, %24.16e)",omega_m,omega_p);
      const double max_omega_size = (1./3.)*(fabs(coef[2]) + fabs(2*radical));
      println("max_omega_size = %24.16e",max_omega_size);
      double max_angle = max_omega_size*d_time;
      // RK3 approximates cos and sin with "cos" and "sin",
      // defined to be the first two nonzero terms of their
      // Taylor series. Setting "cos"^2+"sin"^2=1 reveals that
      // max_angle <= sqrt(3) for RK3 stability. (<=2*sqrt(2) for RK4).
      println("electro-momentum cfl = %24.16e",max_angle*(1./sqrt(3.)));
    }
    gsl_poly_solve_cubic(coef,om);
    if(true) // polish roots
    {
      double OM[3]; for(int i=0;i<3;i++) OM[i]=om[i];
      Polynomial poly(coef, 3);
      poly.polish_roots(om, 3);
      // polishing seems to consistently help only the
      // middle root (which is smaller by two orders of magnitude)
      if(false)
      {
        for(int i=0;i<3;i++)
        {
          double pom = poly.eval(om[i]);
          double pOM = poly.eval(OM[i]);
          double ratio = pom/pOM;
          //const char * success = fcmp(pOM,0.,1e-10) ? "yes!" : "no";
          const char * success = fabs(ratio) < 1. ? "yes!" : "no";
          //dprintf("for OM[%d] =%24.16e p(OM) = %24.16e",
          //  i,OM[i],poly.eval(OM[i]));
          dprintf("for om[%d] =%24.16e "
                  "pOM = %10.2e, "
                  "pom = %10.2e, "
                  "pom/pOM=%6.3f %s",
            i,om[i],pOM,pom,ratio, success);
        }
      }
    }
    if(debug3)
    {
      //println("angles:")
      //printvn(om[0]*d_time)
      //printvn(om[1]*d_time)
      //printvn(om[2]*d_time)
      println("frequencies and periods:")
      printvc(om[0]); printvn(2*pi/om[0]);
      printvc(om[1]); printvn(2*pi/om[1]);
      printvc(om[2]); printvn(2*pi/om[2]);
    }

    // determine a positively oriented orthonormal system
    // exb, e, b, that is,
    // a1, a2, a3 (where a1 cross a2 = a3)
    //
    double a[4][4];
    set_B_orth_system(a,B1,B2,B3,Bmag);

    // rescale the velocities as needed to have an antisymmetric system.
    //
    const double ui0_inv = sqrt(rho_i*one_over_eps);
    const double ue0_inv = sqrt(rho_e*one_over_eps);
    double ui[4],ue[4],E[4];
    ui[1] = prim.get(i,_rho_i+1)*ui0_inv;
    ui[2] = prim.get(i,_rho_i+2)*ui0_inv;
    ui[3] = prim.get(i,_rho_i+3)*ui0_inv;
    ue[1] = prim.get(i,_rho_e+1)*ue0_inv;
    ue[2] = prim.get(i,_rho_e+2)*ue0_inv;
    ue[3] = prim.get(i,_rho_e+3)*ue0_inv;
    E[1] = _E1 ? prim.get(i,_E1) : 0.;
    E[2] = _E2 ? prim.get(i,_E2) : 0.;
    E[3] = prim.get(i,_E3);

    // project components of E, ui, and ue onto the orthonormal system
    //
    // Theory:
    // (u = u_a[j]*a[j] = u_e[j]*e[j]).a[k] says
    //  u_a[k] = u_e[j]*(e[j].a[k])/(a[k].a[k])
    //         = u_e[j]*a[k][j]
    double ui_a[4], ue_a[4], E_a[4];
    //for(int k=1;k<=3;k++)
    //{
    //  ui_a[k]=0.; ue_a[k]=0.; E_a[k]=0.;
    //  for(int j=1;j<=3;j++)
    //  {
    //    ui_a[k] += a[k][j]*ui[j];
    //    ue_a[k] += a[k][j]*ue[j];
    //    E_a[k]  += a[k][j]*E[j];
    //  }
    //}
    // writing it out is more efficient and can be optimized
    ui_a[1] = a[1][1]*ui[1] + a[1][2]*ui[2] + a[1][3]*ui[3];
    ui_a[2] = a[2][1]*ui[1] + a[2][2]*ui[2] + a[2][3]*ui[3];
    ui_a[3] = a[3][1]*ui[1] + a[3][2]*ui[2] + a[3][3]*ui[3];
    ue_a[1] = a[1][1]*ue[1] + a[1][2]*ue[2] + a[1][3]*ue[3];
    ue_a[2] = a[2][1]*ue[1] + a[2][2]*ue[2] + a[2][3]*ue[3];
    ue_a[3] = a[3][1]*ue[1] + a[3][2]*ue[2] + a[3][3]*ue[3];
     E_a[1] = a[1][1]* E[1] + a[1][2]* E[2] + a[1][3]* E[3];
     E_a[2] = a[2][1]* E[1] + a[2][2]* E[2] + a[2][3]* E[3];
     E_a[3] = a[3][1]* E[1] + a[3][2]* E[2] + a[3][3]* E[3];

    // compute eigensolution components
    // for each eigenvalue for subspace orthogonal to b
    //
    double bibe[3], Oibe[3], Oebi[3], norm_sq_inv[3];
    // I should rewrite this to reduce square roots
    const double Oi = sqrt(Oi2);
    const double Oe = sqrt(Oe2);
    for(int k=0;k<3;k++)
    {
      const double bi = gi + om[k];
      const double be = ge - om[k];
      bibe[k] = bi*be;
      Oibe[k] = Oi*be;
      Oebi[k] = Oe*bi;
      //
      // normalizing the eigenvectors is expensive because of sqrt,
      // but I'm afraid of overflow/underflow otherwise
      //
      const double norm_sq = 
          bibe[k]*bibe[k]
        + Oibe[k]*Oibe[k]
        + Oebi[k]*Oebi[k];
      const double norm = sqrt(norm_sq);
      const double norm_inv = 1./norm;
      bibe[k] *= norm_inv;
      Oibe[k] *= norm_inv;
      Oebi[k] *= norm_inv;
      // we have renormalized, so:
      norm_sq_inv[k] = 1.;
    }

    // project the initial conditions onto the initial eigensolutions
    // to determine the coefficients of the expansion
    //
    double c[3][3];
    double* c0=c[0]; // use for parallel component
    double* c1=c[1]; // 1st perpendicular component
    double* c2=c[2]; // 2nd perpendicular component
    for(int k=0;k<3;k++)
    {
      c1[k] =(bibe[k]* E_a[2]
            + Oibe[k]*ui_a[1]
            + Oebi[k]*ue_a[1])*norm_sq_inv[k];
      c2[k] =(bibe[k]* E_a[1]
            - Oibe[k]*ui_a[2]
            - Oebi[k]*ue_a[2])*norm_sq_inv[k];
    }
    //const double Op = sqrt(Op2);
    const double Op_inv = 1./Op; // normalize
    const double Oe_Op = Oe*Op_inv;
    const double Oi_Op = Oi*Op_inv;
    c0[0] = Oe_Op*ui_a[3] + Oi_Op*ue_a[3];
    c0[1] = Oi_Op*ui_a[3] - Oe_Op*ue_a[3];
    c0[2] = E_a[3];

    // evaluate the eigensolutions at time d_time
    //
    // declare and initialize solutions in basis a at time t
    double E_t_a [4], ui_t_a[4], ue_t_a[4];
    for(int k=1;k<=3;k++)
    {
      E_t_a [k]=0.;
      ui_t_a[k]=0.;
      ue_t_a[k]=0.;
    }
    // the perpendicular component
    for(int j=0;j<3;j++)
    {
      const double om_dt = om[j]*d_time;
      const double sin_om_dt = sin(om_dt);
      const double cos_om_dt = cos(om_dt);
      const double rotc1 =  c1[j]*cos_om_dt + c2[j]*sin_om_dt;
      const double rotc2 = -c1[j]*sin_om_dt + c2[j]*cos_om_dt;
      E_t_a [1] += bibe[j]*rotc2;
      E_t_a [2] += bibe[j]*rotc1;
      ui_t_a[1] += Oibe[j]*rotc1;
      ui_t_a[2] -= Oibe[j]*rotc2;
      ue_t_a[1] += Oebi[j]*rotc1;
      ue_t_a[2] -= Oebi[j]*rotc2;
    }
    //
    // the parallel component
    //
    const double Op_dt = Op*d_time;
    const double sin_Op_dt = sin(Op_dt);
    const double cos_Op_dt = cos(Op_dt);
    const double rotc1 =  c0[1]*cos_Op_dt + c0[2]*sin_Op_dt;
    const double rotc2 = -c0[1]*sin_Op_dt + c0[2]*cos_Op_dt;
    E_t_a [3] =                     rotc2;
    ui_t_a[3] = c0[0]*Oe_Op + Oi_Op*rotc1;
    ue_t_a[3] = c0[0]*Oi_Op - Oe_Op*rotc1;

    // project components of E_t, ui_t, and ue_t back onto the standard basis
    //
    // Theory:
    // (u = u_a[j]*a[j] = u_e[j]*e[j]).e[k] says
    //  u_e[k] = u_a[j]*(a[j].e[k])
    //         = u_a[j]*a[j][k]
    double ui_t[4], ue_t[4], E_t[4];
    //for(int k=1;k<=3;k++)
    //{
    //  ui_t[k]=0.; ue_t[k]=0.; E_t[k]=0.;
    //  for(int j=1;j<=3;j++)
    //  {
    //    ui_t[k] += a[j][k]*ui_t_a[j];
    //    ue_t[k] += a[j][k]*ue_t_a[j];
    //    E_t[k]  += a[j][k]*E_t_a[j];
    //  }
    //}
    // write it out so it can be optimized
    // (or is the optimizer possibly smart enough to do that automatically?)
    ui_t[1] = a[1][1]*ui_t_a[1]+a[2][1]*ui_t_a[2]+a[3][1]*ui_t_a[3];
    ui_t[2] = a[1][2]*ui_t_a[1]+a[2][2]*ui_t_a[2]+a[3][2]*ui_t_a[3];
    ui_t[3] = a[1][3]*ui_t_a[1]+a[2][3]*ui_t_a[2]+a[3][3]*ui_t_a[3];
    ue_t[1] = a[1][1]*ue_t_a[1]+a[2][1]*ue_t_a[2]+a[3][1]*ue_t_a[3];
    ue_t[2] = a[1][2]*ue_t_a[1]+a[2][2]*ue_t_a[2]+a[3][2]*ue_t_a[3];
    ue_t[3] = a[1][3]*ue_t_a[1]+a[2][3]*ue_t_a[2]+a[3][3]*ue_t_a[3];
     E_t[1] = a[1][1]* E_t_a[1]+a[2][1]* E_t_a[2]+a[3][1]* E_t_a[3];
     E_t[2] = a[1][2]* E_t_a[1]+a[2][2]* E_t_a[2]+a[3][2]* E_t_a[3];
     E_t[3] = a[1][3]* E_t_a[1]+a[2][3]* E_t_a[2]+a[3][3]* E_t_a[3];

    // copy out the solution, restoring the velocity scale
    const double ui0 = 1./ui0_inv;
    const double ue0 = 1./ue0_inv;
    const double ui_new1 = ui_t[1]*ui0;
    const double ui_new2 = ui_t[2]*ui0;
    const double ui_new3 = ui_t[3]*ui0;
    const double ue_new1 = ue_t[1]*ue0;
    const double ue_new2 = ue_t[2]*ue0;
    const double ue_new3 = ue_t[3]*ue0;
    if(debug3)
    {
      println("\n ui:")
      printvn(ui_new1);
      printvn(ui_new2);
      printvn(ui_new3);
      println(" ue:")
      printvn(ue_new1);
      printvn(ue_new2);
      printvn(ue_new3);
      println(" E:")
      printvn(E_t[1]);
      printvn(E_t[2]);
      printvn(E_t[3]);
    }

    prim.set(i,_rho_i+1, ui_new1);
    prim.set(i,_rho_i+2, ui_new2);
    prim.set(i,_rho_i+3, ui_new3);
    prim.set(i,_rho_e+1, ue_new1);
    prim.set(i,_rho_e+2, ue_new2);
    prim.set(i,_rho_e+3, ue_new3);
    if(_E1) prim.set(i,_E1, E_t[1]); else assert_almost_eq(E_t[1],0.);
    if(_E2) prim.set(i,_E2, E_t[2]); else assert_almost_eq(E_t[2],0.);
    prim.set(i,_E3, E_t[3]);
  }
}

// rotate pressure tensor using exact solution.
//
// trapezoid rule would be less accurate and probably about as expensive;
//
// cost of cos, sin, and sqrt is each about cost of 20 multiplications
// cost of division is like 10 multiplications
//
void SourceSolver::rotate_pressure_tensor(
  dTensor2& prim,
  double d_time,
  double q_over_m,
  int moffset,
  const int _B1_s,
  const int _B2_s,
  const int _B3_s)
{
  using namespace TenMomentComponentID;
  const int _P11_s = _P11 + moffset;
  const int _P12_s = _P12 + moffset;
  const int _P13_s = _P13 + moffset;
  const int _P22_s = _P22 + moffset;
  const int _P23_s = _P23 + moffset;
  const int _P33_s = _P33 + moffset;

  for(int i=1;i<=prim.getsize(1);i++)
  {
    // get the magnetic field
    const double B1 = prim.get(i,_B1_s);
    const double B2 = prim.get(i,_B2_s);
    const double B3 = _B3_s>0 ? prim.get(i,_B3_s) : 0.;
    const double Bmag = sqrt(B1*B1 + B2*B2 + B3*B3);
    if(Bmag==0.) return; // nothing to do

    // determine the angle to rotate.
    const double theta_over_Bmag = q_over_m*d_time;
    const double theta = theta_over_Bmag*Bmag;

    if(debug3){
      const double gyrofrequency = q_over_m*Bmag;
      printf("\n");
      printvn(d_time)
      printvn(gyrofrequency);
      printvn(gyrofrequency*d_time);
      printvn(theta)
    }

    // To relate rotated to unrotated coordinates
    // we rotate the standard basis vectors.

    // Theory: Suppose we wish to rotate the vector u by R = theta*b.
    //
    // Let u_para = u.b*b.
    // Let u_perp = u-u_para.
    // So u = u_perp + u_para.
    // Let uxb = u cross b.
    // Let th = theta = Rmag.
    //
    // The rotated vector is
    //   u = u_para + (cos th)*u_perp + (sin th)*uxb.
    //     = u*(cos th) + (1-cos th)* u_para + (sin th)*uxb.
    // That is,
    //   u = u*(cos th) + 2(sin^2(th/2))*u.b*b + 2*cos(th/2)*uxb*sin(th/2)
    //     = u*(cos th) + 0.5(sinc^2(th/2))*u.R*R + cos(th/2)*uxR*sinc(th/2).
    //
    double c_th, r_scale, r_factor, rr_factor;
    if(true)
    // r = b = rotation direction vector
    // issue: when Bmag is small will (b1, b2, b3) still be of unit length?
    {
      const double s_th = sin(theta);
      const double Bmag_inv = 1./Bmag;
    //
      c_th = cos(theta);
      r_scale = Bmag_inv;
      r_factor = s_th;
      rr_factor = (1-c_th);
    }
    else
    // use sinc to avoid normalizing
    // r = R = rotation vector
    {
      const double th_2 = theta*0.5;
      const double s_th_2 = sin(th_2);
      const double c_th_2 = cos(th_2);
      const double sc_th_2 = s_th_2/th_2;

      c_th = c_th_2*c_th_2-s_th_2*s_th_2;
      r_scale = theta_over_Bmag;
      r_factor = sc_th_2*c_th_2;
      rr_factor = 0.5*sc_th_2*sc_th_2;
    }
    const double r1 = B1*r_scale;
    const double r2 = B2*r_scale;
    const double r3 = B3*r_scale;

    // Let e_i  be the i-th elementary basis vector, and
    // let e_j' be the j-th rotated elementary basis vector.
    double Rot[4][4]; // Rot[i][j] = e_i.e_j'
    //
    const double rot1 = r_factor*r1;
    const double rot2 = r_factor*r2;
    const double rot3 = r_factor*r3;
    //
    const double rr_term11 = rr_factor*r1*r1;
    const double rr_term22 = rr_factor*r2*r2;
    const double rr_term33 = rr_factor*r3*r3;
    //
    const double rr_term23 = rr_factor*r2*r3;
    const double rr_term31 = rr_factor*r3*r1;
    const double rr_term12 = rr_factor*r1*r2;
    //
    Rot[1][1] = rr_term11 + c_th;
    Rot[2][2] = rr_term22 + c_th;
    Rot[3][3] = rr_term33 + c_th;
    //
    Rot[2][3] = rr_term23 + rot1;
    Rot[3][1] = rr_term31 + rot2;
    Rot[1][2] = rr_term12 + rot3;
    //
    Rot[3][2] = rr_term23 - rot1;
    Rot[1][3] = rr_term31 - rot2;
    Rot[2][1] = rr_term12 - rot3;
    //
    // To rotate the pressure tensor we use that
    // the components are invariant in the rotating
    // frame.
    //
    //   P_rot = P_{i,j}*e_i'*e_j', so
    //   P_rot_{m,n} = (e_m.e_i')*P_{i,j}*(e_j'.e_n)

    const double P11 = prim.get(i, _P11_s);
    const double P12 = prim.get(i, _P12_s);
    const double P13 = prim.get(i, _P13_s);
    const double P22 = prim.get(i, _P22_s);
    const double P23 = prim.get(i, _P23_s);
    const double P33 = prim.get(i, _P33_s);
    double P[4][4];
    P[1][1]=P11; P[1][2]=P12; P[1][3]=P13; 
    P[2][1]=P12; P[2][2]=P22; P[2][3]=P23; 
    P[3][1]=P13; P[3][2]=P23; P[3][3]=P33; 

    // This costs 6*(9+3) = 72 multiplications.
    double PR[4][4];
    for(int n=1;n<=3;n++)
    for(int m=1;m<=n;m++) // we only use these components
    // for(int m=1;m<=3;m++)
    {
      PR[m][n]=0.;
      for(int i=1;i<=3;i++)
      {
        double acc=0.;
        for(int j=1;j<=3;j++)
        {
          acc+=P[i][j]*Rot[n][j];
          // PR[m][n]+=Rot[m][i]*P[i][j]*Rot[n][j];
        }
        // factored this out for computational efficiency
        PR[m][n]+=Rot[m][i]*acc;
      }
    }
    prim.set(i, _P11_s, PR[1][1]);
    prim.set(i, _P12_s, PR[1][2]);
    prim.set(i, _P13_s, PR[1][3]);
    prim.set(i, _P22_s, PR[2][2]);
    prim.set(i, _P23_s, PR[2][3]);
    prim.set(i, _P33_s, PR[3][3]);
  }
}

bool SourceSolver::relax_gas(
  dTensor2& prim,
  double d_time,
  GasModelType::Enum gasModelType,
  SpeciesType::Enum speciesType,
  int moffset,
  const int _B1_s,
  const int _B2_s,
  const int _B3_s)
{
  // if(gasModelType==GasModelType::gas05) return false;
  assert(gasModelType!=GasModelType::gas05);

  // determine whether/how to isotropize
  //
  IsoPeriodType::Enum iso_period_type
    = plasmaParams.get_iso_period_type();
  //double in_iso_period;
  //double pcl_mass;
  double in_iso_rate;
  double pcl_mass_inv;
  double prandtl_number;
  switch(speciesType)
  {
    case SpeciesType::ION:
      in_iso_rate = plasmaParams.get_ion_base_iso_rate();
      prandtl_number = plasmaParams.get_ion_prandtl_number();
      pcl_mass_inv = plasmaParams.get_ion_mass_inv();
      break;
    case SpeciesType::ELC:
      in_iso_rate = plasmaParams.get_elc_base_iso_rate();
      prandtl_number = plasmaParams.get_elc_prandtl_number();
      pcl_mass_inv = plasmaParams.get_elc_mass_inv();
      break;
    default:
      unsupported_value_error(speciesType);
  }
  const bool isotropize = (in_iso_rate > 0.);
  if(!isotropize)
  {
    dprintf("not isotropizing");
    return false;
  }

  assert(gasModelType==GasModelType::gas10
      || gasModelType==GasModelType::gas20);
  using namespace TenMomentComponentID;
  const int _rho_s = _rho + moffset;
  const int _u1_s  =  _u1 + moffset;
  const int _u2_s  =  _u2 + moffset;
  const int _u3_s  =  _u3 + moffset;
  const int _P11_s = _P11 + moffset;
  const int _P12_s = _P12 + moffset;
  const int _P13_s = _P13 + moffset;
  const int _P22_s = _P22 + moffset;
  const int _P23_s = _P23 + moffset;
  const int _P33_s = _P33 + moffset;

  //const bool physically_isotropize =
  //   (iso_period_type==IsoPeriodType::det
  //  ||iso_period_type==IsoPeriodType::trace);

  // instantaneously relax toward isotropy
  //
  // Loop over each point
  const int numpts = prim.getsize(1);
  for (int i=1; i<=numpts; i++)
  {
      const double& rho = prim.get(i, _rho_s);
      // change this: should not make a numerical assert
      // statement; instead should set an error flag which
      // is later checked and results in an exit with
      // proper behavior (e.g. dumping current state)
      // (probably easiest would be to throw and catch an error)
      assert_gt(rho,0.);
      const double u1  = prim.get(i, _u1_s);
      const double u2  = prim.get(i, _u2_s);
      const double u3  = prim.get(i, _u3_s);
            double P11 = prim.get(i, _P11_s);
            double P12 = prim.get(i, _P12_s);
            double P13 = prim.get(i, _P13_s);
            double P22 = prim.get(i, _P22_s);
            double P23 = prim.get(i, _P23_s);
            double P33 = prim.get(i, _P33_s);
      //
      double p = (P11 + P22 + P33)*(1./3.);

      // determine portion to keep
      double iso_rate = Relax::get_iso_rate(
        iso_period_type, in_iso_rate, pcl_mass_inv,
        rho, p, P11, P12, P13, P22, P23, P33);
      double portion_to_keep = Relax::get_portion_to_keep(iso_rate,d_time);

      // relax toward isotropic pressure
      {
          const double portion_to_isotropize = 1.-portion_to_keep;
          Relax::isotropizePressure(portion_to_isotropize, p,
               P11, P12, P13, P22, P23, P33,
               P11, P12, P13, P22, P23, P33);
          const bool laterally_isotropize = false;
          if(laterally_isotropize && in_iso_rate >= 0. && in_iso_rate<=DBL_MAX)
          {
              const double& B1 = prim.get(i,_B1_s);
              const double& B2 = prim.get(i,_B2_s);
              const double& B3 = _B3_s>0 ? prim.get(i,_B3_s) : 0.;
              const double B = sqrt(B1*B1 + B2*B2 + B3*B3);
              const double gyrofreq = B*pcl_mass_inv;
              // const double base_lat_iso_timescale
              //   = plasmaParams.get_base_lat_iso_timescale();
              const double base_lat_iso_timescale = 1.;
              // it is not clear to me what lat_iso_rate should be
              const double lat_iso_rate = gyrofreq*
                (base_lat_iso_timescale*iso_rate);
              const double portion_to_gyrotropize = 1.-exp(-d_time*lat_iso_rate);
              
              if(portion_to_gyrotropize > 0.)
              {
                  Relax::laterallyIsotropizePressure(
                       portion_to_gyrotropize,
                       B1,B2,B3,B,
                       P11, P12, P13, P22, P23, P33,
                       P11, P12, P13, P22, P23, P33);
              }
          }
      }

      prim.set(i, _P11_s, P11);
      prim.set(i, _P12_s, P12);
      prim.set(i, _P13_s, P13);
      prim.set(i, _P22_s, P22);
      prim.set(i, _P23_s, P23);
      prim.set(i, _P33_s, P33);

      if(gasModelType==GasModelType::gas20)
      {
          // Should modify portion_to_keep as necessary to keep
          // states in the region of hyperbolicity.
          //
          // Want to enforce |Qxxx| <= alpha*rho*sqrt(Theta_xx)^3 for some alpha.
          //
          // Use that tr(adj(Theta))/det(Theta) = tr(Theta\inverse) >= Theta_xx\inv
          // and that Qxxx^2 <= |tr(Q)|^2.
          //
          // So it is enough to enforce that 
          //   |tr(Q)| <= alpha*rho*sqrt(det(Theta)/tr(adj(Theta)))^3
          // Crudely: compute both sides and multiply Q by
          //   a smoothed version of the ratio if appropriate.
          //
          // always use the determinant (to maintain positivity)
          if(iso_period_type!= IsoPeriodType::det)
          {
              iso_rate = prandtl_number*Relax::get_iso_rate(
                IsoPeriodType::det, in_iso_rate, pcl_mass_inv,
                rho, p, P11, P12, P13, P22, P23, P33);
              portion_to_keep = Relax::get_portion_to_keep(iso_rate,d_time);
          }
          else if(prandtl_number!=1.)
          {
              iso_rate = prandtl_number*iso_rate;
              portion_to_keep = Relax::get_portion_to_keep(iso_rate,d_time);
          }
          using namespace TwentyMomentComponentID;
          const int _Q111_s = _Q111 + moffset;
          const int _Q112_s = _Q112 + moffset;
          const int _Q113_s = _Q113 + moffset;
          const int _Q122_s = _Q122 + moffset;
          const int _Q123_s = _Q123 + moffset;
          const int _Q133_s = _Q133 + moffset;
          const int _Q222_s = _Q222 + moffset;
          const int _Q223_s = _Q223 + moffset;
          const int _Q233_s = _Q233 + moffset;
          const int _Q333_s = _Q333 + moffset;
          prim.fetch(i, _Q111_s) *= portion_to_keep;
          prim.fetch(i, _Q112_s) *= portion_to_keep;
          prim.fetch(i, _Q113_s) *= portion_to_keep;
          prim.fetch(i, _Q122_s) *= portion_to_keep;
          prim.fetch(i, _Q123_s) *= portion_to_keep;
          prim.fetch(i, _Q133_s) *= portion_to_keep;
          prim.fetch(i, _Q222_s) *= portion_to_keep;
          prim.fetch(i, _Q223_s) *= portion_to_keep;
          prim.fetch(i, _Q233_s) *= portion_to_keep;
          prim.fetch(i, _Q333_s) *= portion_to_keep;
      }
  }
  return true;
}

bool SourceSolver::isotropize_pressure(
  dTensor2& prim,
  double d_time,
  SpeciesType::Enum speciesType,
  int moffset,
  const int _B1_s,
  const int _B2_s,
  const int _B3_s)
{
  return SourceSolver::relax_gas(prim, d_time,
    GasModelType::gas10, speciesType, moffset, _B1_s, _B2_s, _B3_s);
}

