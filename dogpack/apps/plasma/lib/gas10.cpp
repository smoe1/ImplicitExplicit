//#define CHECK_BOUNDS
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "PlasmaParams.h"
#include "gas10.h"
#include "Limiters.h" // for PositivityReturnType

const bool check_eigenstructure = true;

// ten-moment gas routines

// cons and prim may reference the same array
void gas10_convert_cons_to_prim(int moffset,
  const dTensor1& cons, dTensor1& prim)
{
  using namespace TenMomentComponentID;
  const int rho_ = moffset+_rho;
  const int M1_ = moffset+_M1;
  const int M2_ = moffset+_M2;
  const int M3_ = moffset+_M3;
  const int P11_ = moffset+_N11;
  const int P12_ = moffset+_N12;
  const int P13_ = moffset+_N13;
  const int P22_ = moffset+_N22;
  const int P23_ = moffset+_N23;
  const int P33_ = moffset+_N33;

  const double rho = cons.get(rho_);
  const double rho_inv = 1./rho;
  const double M1  = cons.get(M1_);
  const double M2  = cons.get(M2_);
  const double M3  = cons.get(M3_);
  const double u1  = M1*rho_inv;
  const double u2  = M2*rho_inv;
  const double u3  = M3*rho_inv;
  const double N11 = cons.get(P11_);
  const double N12 = cons.get(P12_);
  const double N13 = cons.get(P13_);
  const double N22 = cons.get(P22_);
  const double N23 = cons.get(P23_);
  const double N33 = cons.get(P33_);

  prim.set(rho_, rho);
  prim.set(M1_, u1);
  prim.set(M2_, u2);
  prim.set(M3_, u3);
  prim.set(P11_, N11 - M1*u1);
  prim.set(P12_, N12 - M1*u2);
  prim.set(P13_, N13 - M1*u3);
  prim.set(P22_, N22 - M2*u2);
  prim.set(P23_, N23 - M2*u3);
  prim.set(P33_, N33 - M3*u3);
}

// cons and prim may reference the same array
void gas10_convert_cons_to_prim(int moffset,
  const dTensor2& cons, dTensor2& prim)
{
  using namespace TenMomentComponentID;
  const int rho_ = moffset+_rho;
  const int M1_ = moffset+_M1;
  const int M2_ = moffset+_M2;
  const int M3_ = moffset+_M3;
  const int P11_ = moffset+_N11;
  const int P12_ = moffset+_N12;
  const int P13_ = moffset+_N13;
  const int P22_ = moffset+_N22;
  const int P23_ = moffset+_N23;
  const int P33_ = moffset+_N33;

  const int numpts = prim.getsize(1);
  assert_eq(cons.getsize(1),numpts);
  for(int i=1;i<=numpts;i++)
  {
    const double rho = cons.get(i, rho_);
    const double rho_inv = 1./rho;
    const double M1  = cons.get(i, M1_);
    const double M2  = cons.get(i, M2_);
    const double M3  = cons.get(i, M3_);
    const double u1  = M1*rho_inv;
    const double u2  = M2*rho_inv;
    const double u3  = M3*rho_inv;
    const double N11 = cons.get(i, P11_);
    const double N12 = cons.get(i, P12_);
    const double N13 = cons.get(i, P13_);
    const double N22 = cons.get(i, P22_);
    const double N23 = cons.get(i, P23_);
    const double N33 = cons.get(i, P33_);

    prim.set(i, rho_, rho);
    prim.set(i, M1_, u1);
    prim.set(i, M2_, u2);
    prim.set(i, M3_, u3);
    prim.set(i, P11_, N11 - M1*u1);
    prim.set(i, P12_, N12 - M1*u2);
    prim.set(i, P13_, N13 - M1*u3);
    prim.set(i, P22_, N22 - M2*u2);
    prim.set(i, P23_, N23 - M2*u3);
    prim.set(i, P33_, N33 - M3*u3);
  }
}

// recall formula for characteristic polynomial:
//
// det(I*c-A) = c^3-c^2*(Adj11+Adj22+Adj33)+c*tr(A)-det(A),
//
//  where Adj11 is the determinant of the minor obtained by
//  deleting from A the first row and first column, etc.
//
void gas10_get_pressure_evals(double* evals,
    double P11,
    double P12,
    double P13,
    double P22,
    double P23,
    double P33)
{
  const double trP = P11+P22+P33;

  // compute the trace of the adjugate
  const double Adj11 = P22*P33 - P23*P23;
  const double Adj22 = P33*P11 - P13*P13;
  const double Adj33 = P11*P22 - P12*P12;
  const double trAdj = Adj11+Adj22+Adj33;

  // compute the determinant
  const double Adj13 = P12*P23 - P22*P13;
  const double Adj23 = P12*P13 - P11*P23;
  const double detP = Adj13*P13 + Adj23*P23 + Adj33*P33;

  // compute the eigenvalues
  int gsl_poly_solve_cubic (const double* p, double *roots);
  double coef[4];
  // det(I*c-P) = c^3-trAdj*c^2+trP*c-det(P)
  coef[0] = -detP;
  coef[1] = trP;
  coef[2] = -trAdj;
  //coef[3] = 1; // this is assumed by the root-finder.
  // x^3 + coef[2]*x^2 + coef[1]*x + coef[0]
  gsl_poly_solve_cubic(coef,evals);
  assert_le(evals[0],evals[1]);
  assert_le(evals[1],evals[2]);
  // if there is a repeated root then it should also be a root
  // of the derivative
  // 3*x^2 + 2*coef[2]*x + coef[1]
  if(debug3)
  {
    const double disc = coef[2]*coef[2]-3*coef[1];
    double sq_disc = 0.;
    if(disc > 0) sq_disc = sqrt(disc);
    const double r1 = (1./3.)*(-coef[2] - sq_disc);
    const double r2 = (1./3.)*(-coef[2] + sq_disc);
    println("roots of derivative:");
    printvn(r1)
    printvn(r2)
  }
}

// eliminate? really just prints pressure evals.
void gas10_print_prim(int moffset, const dTensor2& prim)
{
  using namespace TenMomentComponentID;
  //const int rho_ = moffset+_rho;
  //const int M1_ = moffset+_M1;
  //const int M2_ = moffset+_M2;
  //const int M3_ = moffset+_M3;
  const int N11_ = moffset+_N11;
  const int N12_ = moffset+_N12;
  const int N13_ = moffset+_N13;
  const int N22_ = moffset+_N22;
  const int N23_ = moffset+_N23;
  const int N33_ = moffset+_N33;

  const int numpts = prim.getsize(1);
  for(int i=1;i<=numpts;i++)
  {
    //const double rho = prim.get(i, _rho);
    //const double u1  = prim.get(i, M1_);
    //const double u2  = prim.get(i, M2_);
    //const double u3  = prim.get(i, M3_);
    const double P11 = prim.get(i, N11_);
    const double P12 = prim.get(i, N12_);
    const double P13 = prim.get(i, N13_);
    const double P22 = prim.get(i, N22_);
    const double P23 = prim.get(i, N23_);
    const double P33 = prim.get(i, N33_);

    double evals[3];
    gas10_get_pressure_evals(evals, P11,P12,P13,P22,P23,P33);
    printvn(evals[0])
    printvn(evals[1])
    printvn(evals[2])
    double p = (1./3.)*(evals[0]+evals[1]+evals[2]);
    printvn(p)
  }
}

// prim and cons may reference the same array
void gas10_convert_prim_to_cons(int moffset,
  const dTensor1& prim, dTensor1& cons)
{
  using namespace TenMomentComponentID;
  const int rho_ = moffset+_rho;
  const int M1_ = moffset+_M1;
  const int M2_ = moffset+_M2;
  const int M3_ = moffset+_M3;
  const int N11_ = moffset+_N11;
  const int N12_ = moffset+_N12;
  const int N13_ = moffset+_N13;
  const int N22_ = moffset+_N22;
  const int N23_ = moffset+_N23;
  const int N33_ = moffset+_N33;

  const double rho = prim.get(rho_);
  const double u1  = prim.get(M1_);
  const double u2  = prim.get(M2_);
  const double u3  = prim.get(M3_);
  const double M1  = u1*rho;
  const double M2  = u2*rho;
  const double M3  = u3*rho;
  const double P11 = prim.get(N11_);
  const double P12 = prim.get(N12_);
  const double P13 = prim.get(N13_);
  const double P22 = prim.get(N22_);
  const double P23 = prim.get(N23_);
  const double P33 = prim.get(N33_);
  //
  cons.set(rho_, rho);
  cons.set(M1_, M1);
  cons.set(M2_, M2);
  cons.set(M3_, M3);
  cons.set(N11_, P11 + M1*u1);
  cons.set(N12_, P12 + M1*u2);
  cons.set(N13_, P13 + M1*u3);
  cons.set(N22_, P22 + M2*u2);
  cons.set(N23_, P23 + M2*u3);
  cons.set(N33_, P33 + M3*u3);
}

// prim and cons may reference the same array
void gas10_convert_prim_to_cons(int moffset,
  const dTensor2& prim, dTensor2& cons)
{
  using namespace TenMomentComponentID;
  const int rho_ = moffset+_rho;
  const int M1_ = moffset+_M1;
  const int M2_ = moffset+_M2;
  const int M3_ = moffset+_M3;
  const int N11_ = moffset+_N11;
  const int N12_ = moffset+_N12;
  const int N13_ = moffset+_N13;
  const int N22_ = moffset+_N22;
  const int N23_ = moffset+_N23;
  const int N33_ = moffset+_N33;

  const int numpts = prim.getsize(1);
  assert_eq(cons.getsize(1),numpts);
  for(int i=1;i<=numpts;i++)
  {
    const double rho = prim.get(i, rho_);
    const double u1  = prim.get(i, M1_);
    const double u2  = prim.get(i, M2_);
    const double u3  = prim.get(i, M3_);
    const double M1  = u1*rho;
    const double M2  = u2*rho;
    const double M3  = u3*rho;
    const double P11 = prim.get(i, N11_);
    const double P12 = prim.get(i, N12_);
    const double P13 = prim.get(i, N13_);
    const double P22 = prim.get(i, N22_);
    const double P23 = prim.get(i, N23_);
    const double P33 = prim.get(i, N33_);
    //
    cons.set(i, rho_, rho);
    cons.set(i, M1_, M1);
    cons.set(i, M2_, M2);
    cons.set(i, M3_, M3);
    cons.set(i, N11_, P11 + M1*u1);
    cons.set(i, N12_, P12 + M1*u2);
    cons.set(i, N13_, P13 + M1*u3);
    cons.set(i, N22_, P22 + M2*u2);
    cons.set(i, N23_, P23 + M2*u3);
    cons.set(i, N33_, P33 + M3*u3);
  }
}

int gas10_check_positivity_cons(int moffset, const dTensor1& Q)
{
  const double rho_epsilon = 0.5*LimiterParams::get_rho_epsilon();
  const double epsilon = 0.5*LimiterParams::get_pressure_epsilon();

  using namespace TenMomentComponentID;
  using namespace PositivityReturnType;
  const double rho = Q.get(moffset+_rho);
  if(rho<=rho_epsilon)
  {
    dprintf("negative density:"
      "\n  moffset = %d,"
      "\n  rho = %24.16e",
      moffset, rho);
    return NEG_DENSITY;
  }

  // confirm positivity indicators at cell average state
  double pos[4];

  const double M1  = Q.get(moffset+_M1);
  const double M2  = Q.get(moffset+_M2);
  const double M3  = Q.get(moffset+_M3);
  const double Ne11 = Q.get(moffset+_N11)-epsilon;
  const double Ne12 = Q.get(moffset+_N12);
  const double Ne13 = Q.get(moffset+_N13);
  const double Ne22 = Q.get(moffset+_N22)-epsilon;
  const double Ne23 = Q.get(moffset+_N23);
  const double Ne33 = Q.get(moffset+_N33)-epsilon;

  const double Pos11 = rho*Ne11 - M1*M1;
  const double Pos12 = rho*Ne12 - M1*M2;
  const double Pos13 = rho*Ne13 - M1*M3;
  const double Pos22 = rho*Ne22 - M2*M2;
  const double Pos23 = rho*Ne23 - M2*M3;
  const double Pos33 = rho*Ne33 - M3*M3;

  const double Adj13 = Pos12*Pos23 - Pos22*Pos13;
  const double Adj23 = Pos12*Pos13 - Pos11*Pos23;
  const double Adj33 = Pos11*Pos22 - Pos12*Pos12;

  const double detPos = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;

  // positive-definite means that all three of these are positive.
  //if(Pos11<0.||Adj33<0.||detPos<0.) retval = false;
  pos[0] = rho;
  pos[1] = Pos11;
  pos[2] = Adj33;
  pos[3] = detPos;
  if( pos[1]<=0. || pos[2]<=0. || pos[3]<=0.)
  {
    if(debug2)
    {
      const double Adj11 = Pos22*Pos33 - Pos23*Pos23;
      const double Adj22 = Pos33*Pos11 - Pos13*Pos13;
      println("positivity violation:");
      printvn(moffset);
      printvn(rho);
      printvn(Pos11);
      printvn(Pos22);
      printvn(Pos33);
      printvn(Adj11);
      printvn(Adj22);
      printvn(Adj33);
      printvn(detPos);
    }

    return NEG_PRESSURE;
  }
  return EVALUATIONS_CONFIRMED_POSITIVITY;
}


inline double _get_entropy_per_mass10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33)
{
  const double u1 = M1/rho;
  const double u2 = M2/rho;
  const double u3 = M3/rho;
  const double P11 = N11 - M1*u1;
  const double P12 = N12 - M1*u2;
  const double P13 = N13 - M1*u3;
  const double P22 = N22 - M2*u2;
  const double P23 = N23 - M2*u3;
  const double P33 = N33 - M3*u3;
  const double detP = det3x3symm(P11, P12, P13, P22, P23, P33);
  const double entropy = 0.5*log(detP/pow(rho,5));
  return entropy;
}
inline double _get_entropy_per_vol10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33)
{
  return rho*_get_entropy_per_mass10(rho,M1,M2,M3,N11,N12,N13,N22,N23,N33);
}
double get_entropy_per_mass10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33)
{
  return _get_entropy_per_mass10(rho,M1,M2,M3,N11,N12,N13,N22,N23,N33);
}
double get_entropy_per_vol10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33)
{
  return _get_entropy_per_vol10(rho,M1,M2,M3,N11,N12,N13,N22,N23,N33);
}

// returns zero on success
int gas10_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *gas_min_speed,
  double *gas_max_speed)
{
  using namespace TenMomentComponentID;
  const double& n1 = nhat.get(1);
  const double& n2 = nhat.get(2);

  // check for inexpensive special cases
  if(n2==0.)
  {
    // Left primitive state
    const double l_rho = Ql.get(moffset+_rho);
    assert_gt(l_rho,0.);
    const double l_rho_inv = 1./l_rho;
    const double l_M1  = Ql.get(moffset+_M1);
    const double l_u1  = l_M1*l_rho_inv;
    const double l_P11 = Ql.get(moffset+_N11) - l_M1*l_u1;
    
    // Left primitive state
    const double r_rho = Qr.get(moffset+_rho);
    assert_gt(r_rho,0.);
    const double r_rho_inv = 1./r_rho;
    const double r_M1  = Qr.get(moffset+_M1);
    const double r_u1  = r_M1*r_rho_inv;
    const double r_P11 = Qr.get(moffset+_N11) - r_M1*r_u1;
    
    // Average primitive states
    const double rho = 0.5e0*(l_rho + r_rho);
    const double u1  = 0.5e0*(l_u1  + r_u1 );
    const double P11 = 0.5e0*(l_P11 + r_P11);
    
    // Sound speeds
    // 
    // speeds for central state
    //
    const double Pnn = P11;
    const double Pnn_rho = Pnn/rho;
    const double cs = sqrt(Pnn_rho);
    const double cf = sq3*cs;
    //
    // speeds for left state
    //
    const double l_Pnn = l_P11;
    //if(l_Pnn<=0.) // debugging
    //{
    //  dprint1(l_Pnn)
    //  return -1;
    //}
    assert_gt(l_rho,0.);
    assert_gt(l_Pnn,0.);
    const double l_Pnn_rho = l_Pnn*l_rho_inv;
    const double l_cs = sqrt(l_Pnn_rho);
    const double l_cf = sq3*l_cs;
    //
    // speeds for right state
    //
    const double r_Pnn = r_P11;
    assert_gt(r_rho,0.);
    assert_gt(r_Pnn,0.);
    const double r_Pnn_rho = r_Pnn*r_rho_inv;
    const double r_cs = sqrt(r_Pnn_rho);
    const double r_cf = sq3*r_cs;
  
    // normal velocities
    const double un   = u1;
    const double l_un = l_u1;
    const double r_un = r_u1;

    (*gas_min_speed) = Min(l_un-l_cf, un-cf);
    (*gas_max_speed) = Max(r_un+r_cf, un+cf);
    return 0;
  }
  else if(n1==0.)
  {
    // Left primitive state
    const double l_rho = Ql.get(moffset+_rho);
    const double l_rho_inv = 1./l_rho;
    const double l_M2  = Ql.get(moffset+_M2);
    const double l_u2  = l_M2*l_rho_inv;
    const double l_P22 = Ql.get(moffset+_N22) - l_M2*l_u2;
    
    // Left primitive state
    const double r_rho = Qr.get(moffset+_rho);
    const double r_rho_inv = 1./r_rho;
    const double r_M2  = Qr.get(moffset+_M2);
    const double r_u2  = Qr.get(moffset+_M2)*r_rho_inv;
    const double r_P22 = Qr.get(moffset+_N22) - r_M2*r_u2;
    
    // Average primitive states
    const double rho = 0.5e0*(l_rho + r_rho);
    const double u2  = 0.5e0*(l_u2  + r_u2 );
    const double P22 = 0.5e0*(l_P22 + r_P22);
    
    // Sound speeds
    // 
    // speeds for central state
    //
    const double Pnn = P22;
    const double Pnn_rho = Pnn/rho;
    const double cs = sqrt(Pnn_rho);
    const double cf = sq3*cs;
    //
    // speeds for left state
    //
    const double l_Pnn = l_P22;
    assert_gt(l_rho,0.);
    //if(l_Pnn<=0.) // debugging
    //{
    //  dprint1(l_Pnn)
    //  return -1;
    //}
    assert_gt(l_Pnn,0.);
    const double l_Pnn_rho = l_Pnn*l_rho_inv;
    const double l_cs = sqrt(l_Pnn_rho);
    const double l_cf = sq3*l_cs;
    //
    // speeds for right state
    //
    const double r_Pnn = r_P22;
    assert_gt(r_rho,0.);
    assert_gt(r_Pnn,0.);
    const double r_Pnn_rho = r_Pnn*r_rho_inv;
    const double r_cs = sqrt(r_Pnn_rho);
    const double r_cf = sq3*r_cs;
  
    // normal velocities
    const double un   = u2;
    const double l_un = l_u2;
    const double r_un = r_u2;

    (*gas_min_speed) = Min(l_un-l_cf, un-cf);
    (*gas_max_speed) = Max(r_un+r_cf, un+cf);
    return 0;
  }
  else
  {
    const double nmag2 = n1*n1 + n2*n2;
    assert_almost_eq(nmag2,1.);
  
    // Left primitive state
    const double l_rho = Ql.get(moffset+_rho);
    const double l_rho_inv = 1./l_rho;
    const double l_M1  = Ql.get(moffset+_M1);
    const double l_M2  = Ql.get(moffset+_M2);
    const double l_u1  = l_M1*l_rho_inv;
    const double l_u2  = l_M2*l_rho_inv;
    const double l_P11 = Ql.get(moffset+_N11) - l_M1*l_u1;
    const double l_P12 = Ql.get(moffset+_N12) - l_M1*l_u2;
    const double l_P22 = Ql.get(moffset+_N22) - l_M2*l_u2;
    
    // Left primitive state
    const double r_rho = Qr.get(moffset+_rho);
    const double r_rho_inv = 1./r_rho;
    const double r_M1  = Qr.get(moffset+_M1);
    const double r_M2  = Qr.get(moffset+_M2);
    const double r_u1  = Qr.get(moffset+_M1)*r_rho_inv;
    const double r_u2  = Qr.get(moffset+_M2)*r_rho_inv;
    const double r_P11 = Qr.get(moffset+_N11) - r_M1*r_u1;
    const double r_P12 = Qr.get(moffset+_N12) - r_M1*r_u2;
    const double r_P22 = Qr.get(moffset+_N22) - r_M2*r_u2;
    
    // Average primitive states
    const double rho = 0.5e0*(l_rho + r_rho);
    const double u1  = 0.5e0*(l_u1  + r_u1 );
    const double u2  = 0.5e0*(l_u2  + r_u2 );
    const double P11 = 0.5e0*(l_P11 + r_P11);
    const double P12 = 0.5e0*(l_P12 + r_P12);
    const double P22 = 0.5e0*(l_P22 + r_P22);
    
    // Sound speeds
    // 
    // speeds for central state
    //
    const double Pnn = n1*n1*P11 + 2*n1*n2*P12 + n2*n2*P22;
    const double Pnn_rho = Pnn/rho;
    const double cs = sqrt(Pnn_rho);
    const double cf = sq3*cs;
    //
    // speeds for left state
    //
    const double l_Pnn = n1*n1*l_P11 + 2*n1*n2*l_P12 + n2*n2*l_P22;
    assert_gt(l_rho,0.);
    assert_gt(l_Pnn,0.);
    const double l_Pnn_rho = l_Pnn*l_rho_inv;
    const double l_cs = sqrt(l_Pnn_rho);
    const double l_cf = sq3*l_cs;
    //
    // speeds for right state
    //
    const double r_Pnn = n1*n1*r_P11 + 2*n1*n2*r_P12 + n2*n2*r_P22;
    assert_gt(r_rho,0.);
    assert_gt(r_Pnn,0.);
    const double r_Pnn_rho = r_Pnn*r_rho_inv;
    const double r_cs = sqrt(r_Pnn_rho);
    const double r_cf = sq3*r_cs;
  
    // normal velocities
    const double un   = n1*u1   + n2*u2;
    const double l_un = n1*l_u1 + n2*l_u2;
    const double r_un = n1*r_u1 + n2*r_u2;

    (*gas_min_speed) = Min(l_un-l_cf, un-cf);
    (*gas_max_speed) = Max(r_un+r_cf, un+cf);
    return 0;
  }
  eprintf("should not get here");
}

void set_entropy_per_vol10(dTensor2& q_new, const dTensor2& qvals,
    int firstIdx, int _entropy_s)
{
  // Loop over each point
  for (int i=1; i<=qvals.getsize(1); i++)
  {
    int idx=firstIdx;
    const double rho = qvals.get(i, idx++);
    const double M1  = qvals.get(i, idx++);
    const double M2  = qvals.get(i, idx++);
    const double M3  = qvals.get(i, idx++);
    const double N11 = qvals.get(i, idx++);
    const double N12 = qvals.get(i, idx++);
    const double N13 = qvals.get(i, idx++);
    const double N22 = qvals.get(i, idx++);
    const double N23 = qvals.get(i, idx++);
    const double N33 = qvals.get(i, idx++);
    const double entropy = _get_entropy_per_vol10(
      rho,M1,M2,M3,N11,N12,N13,N22,N23,N33);
    q_new.set(i, _entropy_s, entropy);
  }
}

inline void gas10_set_cons_rght_eig(double rr[11][11],
  const double cs,
  const double cf,
  const double rho,
  const double M1,
  const double M2,
  const double M3,
  const double P11,
  const double P12,
  const double P13,
  const double P22,
  const double P23,
  const double P33,
  const double u1,
  const double u2,
  const double u3,
  const double uu11,
  const double uu12,
  const double uu13,
  const double uu23,
  const double uu22,
  const double uu33)
{
  const int m1=1;
  const int m2=2;
  const int m3=3;
  const int m4=4;
  const int m5=5;
  const int m6=6;
  const int m7=7;
  const int m8=9; // these defs could be switched if also done at
  const int m9=8; // the_other_places_in_this_file_where_this_string_occurs.
  const int m10=10;

  // fast waves: 10 and 1

  const double drho = rho*P11;
  const double du1 = cf*P11;
  const double du2 = cf*P12;
  const double du3 = cf*P13;
  const double P11_t_3 = 3*P11;
  const double dP11 = P11_t_3*P11;
  const double dP12 = P11_t_3*P12;
  const double dP13 = P11_t_3*P13;
  const double dP23 = P11*P23 + 2*P12*P13;
  const double dP22 = P11*P22 + 2*P12*P12;
  const double dP33 = P11*P33 + 2*P13*P13;

  double a[11], b[11];
  a[m2] = drho*u1;
  a[m3] = drho*u2;
  a[m4] = drho*u3;
  a[m5] = drho*uu11 + dP11;
  a[m6] = drho*uu12 + dP12;
  a[m7] = drho*uu13 + dP13;
  a[m8] = drho*uu23 + dP23;
  a[m9] = drho*uu22 + dP22;
  a[m10]= drho*uu33 + dP33;

  b[m2] = rho*du1;
  b[m3] = rho*du2;
  b[m4] = rho*du3;
  b[m5] = 2*M1*du1;
  b[m6] = (M2*du1+M1*du2);
  b[m7] = (M3*du1+M1*du3);
  b[m8] = (M3*du2+M2*du3);
  b[m9] = 2*M2*du2;
  b[m10]= 2*M3*du3;

  // 10 - right eigenvector

  rr[m1][10]   = drho;
  rr[m2][10]   = a[m2] + b[m2];
  rr[m3][10]   = a[m3] + b[m3];
  rr[m4][10]   = a[m4] + b[m4];
  rr[m5][10]   = a[m5] + b[m5];
  rr[m6][10]   = a[m6] + b[m6];
  rr[m7][10]   = a[m7] + b[m7];
  rr[m8][10]   = a[m8] + b[m8];
  rr[m9][10]   = a[m9] + b[m9];
  rr[m10][10]  = a[m10]+ b[m10];

  // 1 - right eigenvector

  rr[m1][1]   = drho;
  rr[m2][1]   = a[m2] - b[m2];
  rr[m3][1]   = a[m3] - b[m3];
  rr[m4][1]   = a[m4] - b[m4];
  rr[m5][1]   = a[m5] - b[m5];
  rr[m6][1]   = a[m6] - b[m6];
  rr[m7][1]   = a[m7] - b[m7];
  rr[m8][1]   = a[m8] - b[m8];
  rr[m9][1]   = a[m9] - b[m9];
  rr[m10][1]  = a[m10]- b[m10];

  // slow waves

  // waves 9 and 2

  a[m9]=2*P12;

  b[m6]=M1*cs;
  b[m8]=M3*cs;
  b[m9]=2*M2*cs;

  const double rho_t_cs = rho*cs;
  // 9 - right eigenvector
  rr[m3][9]   = rho_t_cs;
  rr[m6][9]   = P11+b[m6];
  rr[m8][9]   = P13+b[m8];
  rr[m9][9]   = a[m9]+b[m9];

  // 2 - right eigenvector
  rr[m3][2]   = -rr[m3][9];
  rr[m6][2]   = P11-b[m6];
  rr[m8][2]   = P13-b[m8];
  rr[m9][2]   = a[m9]-b[m9];

  // waves 8 and 3

  b[m7] = M1*cs;
  b[m8] = M2*cs;
  b[m10]=2*M3*cs;
  // 8 - right eigenvector
  rr[m4][8]   = rho_t_cs;
  rr[m7][8]   = P11   + b[m7];
  rr[m8][8]   = P12   + b[m8];
  rr[m10][8]  = 2*P13 + b[m10];

  // 3 - right eigenvector
  rr[m4][3]   = -rr[m4][8];
  rr[m7][3]   = P11   - b[m7];
  rr[m8][3]   = P12   - b[m8];
  rr[m10][3]  = 2*P13 - b[m10];
 
  // waves with speed zero

  // 4 - right eigenvector
  rr[m1][4]   = 1.;
  rr[m2][4]   = u1;
  rr[m3][4]   = u2;
  rr[m4][4]   = u3;
  rr[m5][4]   = uu11;
  rr[m6][4]   = uu12;
  rr[m7][4]   = uu13;
  rr[m8][4]   = uu23;
  rr[m9][4]   = uu22;
  rr[m10][4]  = uu33;

  // 5 - right eigenvector
  rr[m8][5]   = 1.;

  // 6 - right eigenvector
  rr[m9][6]   = 1.;

  // 7 - right eigenvector
  rr[m10][7]  = 1.;

  // set components that are never actually used
  // (could initialize these to NAN)
  if(check_eigenstructure) // false to optimize
  {
    const int ZERO=0.;
    rr[m1][9]   = ZERO;
    rr[m2][9]   = ZERO;
    rr[m4][9]   = ZERO;
    rr[m5][9]   = ZERO;
    rr[m7][9]   = ZERO;
    rr[m10][9]  = ZERO;

    rr[m1][2]   = ZERO;
    rr[m2][2]   = ZERO;
    rr[m4][2]   = ZERO;
    rr[m5][2]   = ZERO;
    rr[m7][2]   = ZERO;
    rr[m10][2]  = ZERO;

    rr[m1][8]   = ZERO;
    rr[m2][8]   = ZERO;
    rr[m3][8]   = ZERO;
    rr[m5][8]   = ZERO;
    rr[m6][8]   = ZERO;
    rr[m9][8]   = ZERO;

    rr[m1][3]   = ZERO;
    rr[m2][3]   = ZERO;
    rr[m3][3]   = ZERO;
    rr[m5][3]   = ZERO;
    rr[m6][3]   = ZERO;
    rr[m9][3]   = ZERO;

    rr[m1][5]   = ZERO;
    rr[m2][5]   = ZERO;
    rr[m3][5]   = ZERO;
    rr[m4][5]   = ZERO;
    rr[m5][5]   = ZERO;
    rr[m6][5]   = ZERO;
    rr[m7][5]   = ZERO;
    rr[m9][5]   = ZERO;
    rr[m10][5]  = ZERO;

    rr[m1][6]   = ZERO;
    rr[m2][6]   = ZERO;
    rr[m3][6]   = ZERO;
    rr[m4][6]   = ZERO;
    rr[m5][6]   = ZERO;
    rr[m6][6]   = ZERO;
    rr[m7][6]   = ZERO;
    rr[m8][6]   = ZERO;
    rr[m10][6]  = ZERO;

    rr[m1][7]   = ZERO;
    rr[m2][7]   = ZERO;
    rr[m3][7]   = ZERO;
    rr[m4][7]   = ZERO;
    rr[m5][7]   = ZERO;
    rr[m6][7]   = ZERO;
    rr[m7][7]   = ZERO;
    rr[m8][7]   = ZERO;
    rr[m9][7]   = ZERO;
  }
}

inline void gas10_rescale_left_eigs(double rl[11][11],
  const double cs,
  const double cf,
  const double rho,
  const double P11,
  const double u1,
  const double uu11)
{
  const int m1=1;
  const int m2=2;
  const int m3=3;
  const int m4=4;
  const int m5=5;
  const int m6=6;
  const int m7=7;
  const int m8=9; // these defs could be switched if also done at
  const int m9=8; // the_other_places_in_this_file_where_this_string_occurs.
  const int m10=10;

  // we compute the portion of the right eigenvector components
  // actually used when taking dot products with corresponding
  // left eigenvectors
  double rr[11][11];
  double rlrr[11];
  double rlrr_inv[11];

  const double rho_t_cs = rho*cs;
  const double drho = rho*P11;
  const double dP11 = 3*P11*P11;
  const double du1 = cf*P11;

  // eigenvectors 10 and 1

  double a[11],b[11];
  a[m2] = drho*u1;
  a[m5] = drho*uu11 + dP11;

  b[m2] = rho*du1;
  b[m5] = rho*2*u1*du1;

  rr[m1][10]   = drho;
  rr[m2][10]   = a[m2] + b[m2];
  rr[m5][10]   = a[m5] + b[m5];
  rlrr[10] = rl[10][m1]*rr[m1][10]
           + rl[10][m2]*rr[m2][10]
           + rl[10][m5]*rr[m5][10];
  rlrr_inv[10] = 1./rlrr[10];
  rl[10][m1] *= rlrr_inv[10];
  rl[10][m2] *= rlrr_inv[10];
  rl[10][m5] *= rlrr_inv[10];

  rr[m1][1]   = drho;
  rr[m2][1]   = a[m2] - b[m2];
  rr[m5][1]   = a[m5] - b[m5];
  rlrr[1] = rl[1][m1]*rr[m1][1]
          + rl[1][m2]*rr[m2][1]
          + rl[1][m5]*rr[m5][1];
  rlrr_inv[1] = 1./rlrr[1];
  rl[1][m1] *= rlrr_inv[1];
  rl[1][m2] *= rlrr_inv[1];
  rl[1][m5] *= rlrr_inv[1];

  // eigenvectors 9 and 2

  b[m6]=rho*u1*(cs);

  rr[m3][9]   = rho_t_cs;
  rr[m6][9]   = P11+b[m6];
  rlrr[9] = rl[9][m3]*rr[m3][9]
          + rl[9][m6]*rr[m6][9];
  rlrr_inv[9] = 1./rlrr[9];
  rl[9][m1] *= rlrr_inv[9];
  rl[9][m2] *= rlrr_inv[9];
  rl[9][m3] *= rlrr_inv[9];
  rl[9][m5] *= rlrr_inv[9];
  rl[9][m6] *= rlrr_inv[9];

  rr[m3][2]   = -rr[m3][9];
  rr[m6][2]   = P11-b[m6];
  rlrr[2] = rl[2][m3]*rr[m3][2]
          + rl[2][m6]*rr[m6][2];
  rlrr_inv[2] = 1./rlrr[2];
  rl[2][m1] *= rlrr_inv[2];
  rl[2][m2] *= rlrr_inv[2];
  rl[2][m3] *= rlrr_inv[2];
  rl[2][m5] *= rlrr_inv[2];
  rl[2][m6] *= rlrr_inv[2];

  // eigenvectors 8 and 3
  //
  b[m7] = rho*u1*(cs);

  rr[m4][8]   = rho_t_cs;
  rr[m7][8]   = P11   + b[m7];
  rlrr[8] = rl[8][m4]*rr[m4][8]
          + rl[8][m7]*rr[m7][8];
  rlrr_inv[8] = 1./rlrr[8];
  rl[8][m1] *= rlrr_inv[8];
  rl[8][m2] *= rlrr_inv[8];
  rl[8][m4] *= rlrr_inv[8];
  rl[8][m5] *= rlrr_inv[8];
  rl[8][m7] *= rlrr_inv[8];

  rr[m4][3]   = -rr[m4][8];
  rr[m7][3]   = P11   - b[m7];
  rlrr[3] = rl[3][m4]*rr[m4][3]
          + rl[3][m7]*rr[m7][3];
  rlrr_inv[3] = 1./rlrr[3];
  rl[3][m1] *= rlrr_inv[3];
  rl[3][m2] *= rlrr_inv[3];
  rl[3][m4] *= rlrr_inv[3];
  rl[3][m5] *= rlrr_inv[3];
  rl[3][m7] *= rlrr_inv[3];

  // waves with speed zero

  // 4 - right eigenvector
  //rr[m1][4]   = 1.; // can leave out this factor
  rr[m2][4]   = u1;
  rr[m5][4]   = uu11;
  rlrr[4] = rl[4][m1]
          + rl[4][m2]*rr[m2][4]
          + rl[4][m5]*rr[m5][4];
  rlrr_inv[4] = 1./rlrr[4];
  rl[4][m1] *= rlrr_inv[4];
  rl[4][m2] *= rlrr_inv[4];
  rl[4][m5] *= rlrr_inv[4];

  // 5 - right eigenvector
  //rr[m8][5]   = 1.;
  rlrr_inv[5] = 1./rl[5][m8];
  rl[5][m1]  *= rlrr_inv[5];
  rl[5][m2]  *= rlrr_inv[5];
  rl[5][m3]  *= rlrr_inv[5];
  rl[5][m4]  *= rlrr_inv[5];
  rl[5][m5]  *= rlrr_inv[5];
  rl[5][m6]  *= rlrr_inv[5];
  rl[5][m7]  *= rlrr_inv[5];
  rl[5][m8]  = 1.; // not actually used

  // 6 - right eigenvector
  //rr[m9][6]   = 1.;
  rlrr_inv[6] = 1./rl[6][m9];
  rl[6][m1]  *= rlrr_inv[6];
  rl[6][m2]  *= rlrr_inv[6];
  rl[6][m3]  *= rlrr_inv[6];
  rl[6][m5]  *= rlrr_inv[6];
  rl[6][m6]  *= rlrr_inv[6];
  rl[6][m9]  = 1.; // not actually used

  // 7 - right eigenvector
  //rr[m10][7]  = 1.;
  rlrr_inv[7] = 1./rl[7][m10];
  rl[7][m1]  *= rlrr_inv[7];
  rl[7][m2]  *= rlrr_inv[7];
  rl[7][m4]  *= rlrr_inv[7];
  rl[7][m5]  *= rlrr_inv[7];
  rl[7][m7]  *= rlrr_inv[7];
  rl[7][m10] = 1.; // not actually used
}

// eigenstructure (translated from fortran)
//
void gas10_set_cons_left_eig(double rl[11][11],
  const double cs,
  const double cf,
  const double rho,
  const double rho_inv,
  const double M1,
  const double M2,
  const double M3,
  const double u1,
  const double u2,
  const double u3,
  const double P11,
  const double P12,
  const double P13,
  const double P22,
  const double P23,
  const double P33,
  const double u1u1,
  const double u1u2,
  const double u1u3,
  const double u2u3,
  const double u2u2,
  const double u3u3)
{
  const int m1=1;
  const int m2=2;
  const int m3=3;
  const int m4=4;
  const int m5=5;
  const int m6=6;
  const int m7=7;
  const int m8=9; // these defs could be switched if also done at
  const int m9=8; // the_other_places_in_this_file_where_this_string_occurs.
  const int m10=10;

  // fast waves: left eigenvectors 1 and 10
  const double cf_u1 = cf*u1;
  const double P11_t_3 = 3*P11;
  double a[11], b[11];
  b[1]=cf*u1u1;
  b[2]=2*cf_u1;
  a[2] = P11_t_3*rho_inv;
  a[1] = -a[2]*u1;

  rl[1][m1]  = a[1]-b[1];
  rl[1][m2]  = a[2]+b[2];
  rl[1][m5]  = -cf;

  rl[10][m1]  = a[1]+b[1];
  rl[10][m2]  = a[2]-b[2];
  rl[10][m5]  = cf;

  // slow waves: left eigenvectors 2 and 9, 3 and 8

  const double du=P11*P11;
  const double dP=cs*P11;

  // left eigenvectors 9 and 2

  const double du1b = -P12*P11;
  const double dP11b= -cs*P12;

  a[1]  = -(du1b*u1+du*u2)*rho_inv;
  a[2]  = du1b*rho_inv;
  a[3]  = du*rho_inv;

  b[1]  = u1u1*dP11b+u1u2*dP;
  b[2]  = -(2*u1*dP11b+u2*dP);
  b[3]  = -(u1*dP);

  rl[9][m1]  = a[1]+b[1];
  rl[9][m2]  = a[2]+b[2];
  rl[9][m3]  = a[3]+b[3];
  rl[9][m5]  = dP11b;
  rl[9][m6]  = dP;

  rl[2][m1]  = a[1]-b[1];
  rl[2][m2]  = a[2]-b[2];
  rl[2][m3]  = a[3]-b[3];
  rl[2][m5]  = -dP11b;
  rl[2][m6]  = -dP;

  // left eigenvectors 8 and 3

  const double du1c = -P13*P11;
  const double dP11c= -cs*P13;

  a[1]  = -(du1c*u1+du*u3)*rho_inv;
  a[2]  = du1c*rho_inv;
  a[4]  = a[3];

  b[1]  = u1u1*dP11c+u1u3*dP;
  b[2]  = -(2*u1*dP11c+u3*dP);
  b[4]  = b[3];

  rl[8][m1]  = a[1]+b[1];
  rl[8][m2]  = a[2]+b[2];
  rl[8][m4]  = a[4]+b[4];
  rl[8][m5]  = dP11c;
  rl[8][m7]  = dP;

  rl[3][m1]  = a[1]-b[1];
  rl[3][m2]  = a[2]-b[2];
  rl[3][m4]  = a[4]-b[4];
  rl[3][m5]  = -dP11c;
  rl[3][m7]  = -dP;


  // waves with speed zero
  // (note that waves 5 and 6 do not match up with James's 6 and 5,
  // which is entirely likely since they share a common eigenvalue.)

  // left eigenvector 4

  rl[4][m1]  = P11_t_3-rho*u1u1;
  rl[4][m2]  = 2*rho*u1;
  rl[4][m5]  = -rho;

  // waves 5, 6, and 7

  const double P11p2m3=P11_t_3*P11;

  // left eigenvector 5

  const double dP11e = 4*P12*P13-P23*P11;
  const double dP12e = -P11_t_3*P13;
  const double dP13e = -P11_t_3*P12;

  rl[5][m1]  = u1u1*dP11e+u1u2*dP12e+u1u3*dP13e+u2u3*P11p2m3;
  rl[5][m2]  = -(2*u1*dP11e+u2*dP12e+u3*dP13e);
  rl[5][m3]  = -(u1*dP12e+u3*P11p2m3);
  rl[5][m4]  = -(u1*dP13e+u2*P11p2m3);
  rl[5][m5]  = dP11e;
  rl[5][m6]  = dP12e;
  rl[5][m7]  = dP13e;
  rl[5][m8]  = P11p2m3;

  // left eigenvector 6

  const double dP11f = 4*P12*P12-P22*P11;
  const double dP12f = -6*P11*P12;

  rl[6][m1]  = u1u1*dP11f+u1u2*dP12f+u2u2*P11p2m3;
  rl[6][m2]  = -(2*u1*dP11f+u2*dP12f);
  rl[6][m3]  = -(u1*dP12f+2*u2*P11p2m3);
  rl[6][m5]  = dP11f;
  rl[6][m6]  = dP12f;
  rl[6][m9]  = P11p2m3;

  // left eigenvector 7

  const double dP11g = 4*P13*P13-P33*P11;
  const double dP13g = -6*P11*P13;

  rl[7][m1]  = u1u1*dP11g+u1u3*dP13g+u3u3*P11p2m3;
  rl[7][m2]  = -(2*u1*dP11g+u3*dP13g);
  rl[7][m4]  = -(u1*dP13g+2*u3*P11p2m3);
  rl[7][m5]  = dP11g;
  rl[7][m7]  = dP13g;
  rl[7][m10] = P11p2m3;

  // set components that are never actually used
  // (could initialize these to NAN)
  if(check_eigenstructure) // false to optimize
  {
    const int ZERO=0.;
    rl[1][m3]  = ZERO;
    rl[1][m4]  = ZERO;
    rl[1][m6]  = ZERO;
    rl[1][m7]  = ZERO;
    rl[1][m8]  = ZERO;
    rl[1][m9]  = ZERO;
    rl[1][m10] = ZERO;

    rl[10][m3] = ZERO;
    rl[10][m4] = ZERO;
    rl[10][m6] = ZERO;
    rl[10][m7] = ZERO;
    rl[10][m8] = ZERO;
    rl[10][m9] = ZERO;
    rl[10][m10]= ZERO;

    rl[9][m4]  = ZERO;
    rl[9][m7]  = ZERO;
    rl[9][m8]  = ZERO;
    rl[9][m9]  = ZERO;
    rl[9][m10] = ZERO;

    rl[2][m4]  = ZERO;
    rl[2][m7]  = ZERO;
    rl[2][m8]  = ZERO;
    rl[2][m9]  = ZERO;
    rl[2][m10] = ZERO;

    rl[8][m3]  = ZERO;
    rl[8][m6]  = ZERO;
    rl[8][m8]  = ZERO;
    rl[8][m9]  = ZERO;
    rl[8][m10] = ZERO;

    rl[3][m3]  = ZERO;
    rl[3][m6]  = ZERO;
    rl[3][m8]  = ZERO;
    rl[3][m9]  = ZERO;
    rl[3][m10] = ZERO;

    rl[4][m3]  = ZERO;
    rl[4][m4]  = ZERO;
    rl[4][m6]  = ZERO;
    rl[4][m7]  = ZERO;
    rl[4][m8]  = ZERO;
    rl[4][m9]  = ZERO;
    rl[4][m10] = ZERO;

    rl[5][m9]  = ZERO;
    rl[5][m10] = ZERO;

    rl[6][m4]  = ZERO;
    rl[6][m7]  = ZERO;
    rl[6][m8]  = ZERO;
    rl[6][m10] = ZERO;

    rl[7][m3]  = ZERO;
    rl[7][m6]  = ZERO;
    rl[7][m8]  = ZERO;
    rl[7][m9]  = ZERO;
  }

  gas10_rescale_left_eigs(rl,
    cs, cf, rho, P11, u1, u1u1);

  // check eigenstructure
  //
  if(check_eigenstructure) // false to optimize
  {
    double rr[11][11];
    double rlrr[11][11];
    gas10_set_cons_rght_eig(rr,
      cs, cf,
      rho,
      M1, M2, M3,
      P11, P12, P13, P22, P23, P33,
      u1, u2, u3,
      u1u1, u1u2, u1u3, u2u3, u2u2, u3u3);
    for(int i=1; i<=10;i++)
    for(int j=1; j<=10;j++)
    {
      rlrr[i][j] = 0.;
      for(int k=1;k<=10;k++)
      {
        rlrr[i][j] += rl[i][k]*rr[k][j];
        
      }
      if(i!=j)
      {
        if(fcmp(rlrr[i][j],0.,1e-14))
        {
          dprintf("rlrr[%d][%d]=%24.16e",i,j,rlrr[i][j]);
        }
        //assert_printf(!fcmp(rlrr[i][j],0.,1e-14),
        //  "rlrr[%d][%d]=%24.16e",i,j,rlrr[i][j]);
      }
      else
      {
        if(fcmp(rlrr[i][j],1.,1e-14))
        {
          dprintf("rlrr[%d][%d]=%24.16e",i,j,rlrr[i][j]);
        }
        //assert_printf(!fcmp(rlrr[i][j],1.,1e-13),
        //  "rlrr[%d][%d]=%24.16e",i,j,rlrr[i][j]);
      }
    }
  }
}

// W is presumed to hold the coefficients of
// a right eigenvector expansion of Q.
//
// use right eigenstructure of flux jacobian
// (derivative of flux with respect to conserved variables)
// to calculate each component of Q.
//
//
void ProjectRightEig_gas10(int ixy, int moffset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
  using namespace TenMomentComponentID;
  const int n_rho = moffset+_rho;
  const int n_M3  = moffset+_M3;
  const int n_N12 = moffset+_N12;
  const int n_N33 = moffset+_N33;

  int n_M1, n_M2, n_N11, n_N13, n_N22, n_N23;

  if(ixy==1)
  {
      n_M1  = moffset+_M1;
      n_M2  = moffset+_M2;

      n_N11 = moffset+_N11;
      n_N13 = moffset+_N13;
      n_N22 = moffset+_N22;
      n_N23 = moffset+_N23;
  }
  else
  {
      assert(ixy==2);

      n_M1  = moffset+_M2;
      n_M2  = moffset+_M1;

      n_N11 = moffset+_N22;
      n_N13 = moffset+_N23;
      n_N22 = moffset+_N11;
      n_N23 = moffset+_N13;
  }

  // get primitive variables.
  //
  const double rho = Q_ave.get(n_rho);
  const double rho_inv = 1./rho;
  const double M1  = Q_ave.get(n_M1);
  const double M2  = Q_ave.get(n_M2);
  const double M3  = Q_ave.get(n_M3);
  const double u1  = M1*rho_inv;
  const double u2  = M2*rho_inv;
  const double u3  = M3*rho_inv;
  const double P11 = Q_ave.get(n_N11) - M1*u1;
  const double P12 = Q_ave.get(n_N12) - M1*u2;
  const double P13 = Q_ave.get(n_N13) - M1*u3;
  const double P22 = Q_ave.get(n_N22) - M2*u2;
  const double P23 = Q_ave.get(n_N23) - M2*u3;
  const double P33 = Q_ave.get(n_N33) - M3*u3;

  // define variables used in calculation of right eigenstructure

  const double cs = sqrt(P11*rho_inv);
  const double cf = sqrt(3.)*cs;

  const double uu11 = u1*u1;
  const double uu12 = u1*u2;
  const double uu13 = u1*u3;
  const double uu23 = u2*u3;
  const double uu22 = u2*u2;
  const double uu33 = u3*u3;

  double rr[11][11];
  gas10_set_cons_rght_eig(rr,
    cs, cf,
    rho, M1, M2, M3,
    P11, P12, P13, P22, P23, P33,
    u1, u2, u3,
    uu11, uu12, uu13, uu23, uu22, uu33);

  // // Q_nk = rr
  // //
  // // make this faster by writing it out
  // // (so I don't have to test for zeros
  // // and so index offsets don't have to be
  // // computed at run time).
  // //
  // for (int m=1; m<=10; m++)
  // for (int k=1; k<=Q.getsize(2); k++)
  // {
  //   acc=0.;
  //   for (int j=1; j<=10; j++)
  //   {
  //     double rr_mj = rr[m,j];
  //     if(rr_mj) acc+=rr_mj*W.get(j+moffset,k);
  //   }
  //   Q.set(m+moffset,k, acc);
  // }

  // by definition
  // m8 is the N23 component and
  // m9 is the N22 component
  const int m1=1;
  const int m2=2;
  const int m3=3;
  const int m4=4;
  const int m5=5;
  const int m6=6;
  const int m7=7;
  const int m8=9; // these defs could be switched if also done at
  const int m9=8; // the_other_places_in_this_file_where_this_string_occurs.
  const int m10=10;

  const int n1 = 1 +moffset;
  const int n2 = 2 +moffset;
  const int n3 = 3 +moffset;
  const int n4 = 4 +moffset;
  const int n5 = 5 +moffset;
  const int n6 = 6 +moffset;
  const int n7 = 7 +moffset;
  const int n8 = 8 +moffset;
  const int n9 = 9 +moffset;
  const int n10= 10+moffset;

  for(int k=1; k<=Q.getsize(2); k++)
  {
    Q.set(n_rho,k,
        rr[m1][1] * W.get(n1,k)
      +             W.get(n4,k)
      + rr[m1][10]* W.get(n10,k));
    Q.set(n_M1,k,
        rr[m2][1] * W.get(n1,k)
      + rr[m2][4] * W.get(n4,k)
      + rr[m2][10]* W.get(n10,k));
    Q.set(n_M2,k,
        rr[m3][1] * W.get(n1,k)
      + rr[m3][2] * W.get(n2,k)
      + rr[m3][4] * W.get(n4,k)
      + rr[m3][9] * W.get(n9,k) // james: n8
      + rr[m3][10]* W.get(n10,k));
    Q.set(n_M3,k,
        rr[m4][1] * W.get(n1,k)
      + rr[m4][3] * W.get(n3,k)
      + rr[m4][4] * W.get(n4,k)
      + rr[m4][8] * W.get(n8,k) // james: n9
      + rr[m4][10]* W.get(n10,k));
    Q.set(n_N11,k,
        rr[m5][1] * W.get(n1,k)
      + rr[m5][4] * W.get(n4,k)
      + rr[m5][10]* W.get(n10,k));
    Q.set(n_N12,k,
        rr[m6][1] * W.get(n1,k)
      + rr[m6][2] * W.get(n2,k)
      + rr[m6][4] * W.get(n4,k)
      + rr[m6][9] * W.get(n9,k) // james: n8
      + rr[m6][10]* W.get(n10,k));
    Q.set(n_N13,k,
        rr[m7][1] * W.get(n1,k)
      + rr[m7][3] * W.get(n3,k)
      + rr[m7][4] * W.get(n4,k)
      + rr[m7][8] * W.get(n8,k) // james: n9
      + rr[m7][10]* W.get(n10,k));
    Q.set(n_N22,k, // by definition m9 is the N22 component
        rr[m9][1] * W.get(n1,k)
      + rr[m9][2] * W.get(n2,k)
      + rr[m9][4] * W.get(n4,k)
      +             W.get(n6,k) // james: n5
      + rr[m9][9] * W.get(n9,k) // james: n8
      + rr[m9][10]* W.get(n10,k));
    Q.set(n_N23,k, // by definition m8 is the N23 component
        rr[m8][1] * W.get(n1,k)
      + rr[m8][2] * W.get(n2,k)
      + rr[m8][3] * W.get(n3,k)
      + rr[m8][4] * W.get(n4,k)
      +             W.get(n5,k) // james: n6
      + rr[m8][8] * W.get(n8,k) // james: n9
      + rr[m8][9] * W.get(n9,k) // james: n8
      + rr[m8][10]* W.get(n10,k));
    Q.set(n_N33,k,
        rr[m10][1] * W.get(n1,k)
      + rr[m10][3] * W.get(n3,k)
      + rr[m10][4] * W.get(n4,k)
      +              W.get(n7,k)
      + rr[m10][8] * W.get(n8,k) // james: n9
      + rr[m10][10]* W.get(n10,k));
  }
}

void ProjectLeftEig_gas10(int ixy, int moffset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
  using namespace TenMomentComponentID;
  const int n_rho = moffset+_rho;
  const int n_M3  = moffset+_M3;
  const int n_N12 = moffset+_N12;
  const int n_N33 = moffset+_N33;

  int n_M1, n_M2, n_N11, n_N13, n_N22, n_N23;

  if(ixy==1)
  {
      n_M1  = moffset+_M1;
      n_M2  = moffset+_M2;

      n_N11 = moffset+_N11;
      n_N13 = moffset+_N13;
      n_N22 = moffset+_N22;
      n_N23 = moffset+_N23;
  }
  else
  {
      assert(ixy==2);

      n_M1  = moffset+_M2;
      n_M2  = moffset+_M1;

      n_N11 = moffset+_N22;
      n_N13 = moffset+_N23;
      n_N22 = moffset+_N11;
      n_N23 = moffset+_N13;
  }

  // get primitive variables.
  //
  const double rho = Q_ave.get(n_rho);
  const double rho_inv = 1./rho;
  const double M1  = Q_ave.get(n_M1);
  const double M2  = Q_ave.get(n_M2);
  const double M3  = Q_ave.get(n_M3);
  const double u1  = M1*rho_inv;
  const double u2  = M2*rho_inv;
  const double u3  = M3*rho_inv;
  const double P11 = Q_ave.get(n_N11) - M1*u1;
  const double P12 = Q_ave.get(n_N12) - M1*u2;
  const double P13 = Q_ave.get(n_N13) - M1*u3;
  const double P22 = Q_ave.get(n_N22) - M2*u2;
  const double P23 = Q_ave.get(n_N23) - M2*u3;
  const double P33 = Q_ave.get(n_N33) - M3*u3;

  // define variables used in calculation of left eigenstructure

  const double cs = sqrt(P11*rho_inv);
  const double cf = sqrt(3.)*cs;

  const double uu11 = u1*u1;
  const double uu12 = u1*u2;
  const double uu13 = u1*u3;
  const double uu23 = u2*u3;
  const double uu22 = u2*u2;
  const double uu33 = u3*u3;

  double rl[11][11];

  gas10_set_cons_left_eig(rl,
    cs, cf,
    rho, rho_inv,
    M1, M2, M3,
    u1, u2, u3,
    P11, P12, P13, P22, P23, P33,
    uu11, uu12, uu13, uu23, uu22, uu33);

  const int m1=1;
  const int m2=2;
  const int m3=3;
  const int m4=4;
  const int m5=5;
  const int m6=6;
  const int m7=7;
  const int m8=9; // these defs could be switched if also done at
  const int m9=8; // the_other_places_in_this_file_where_this_string_occurs.
  const int m10=10;

  const int n1 = moffset + 1;
  const int n2 = moffset + 2;
  const int n3 = moffset + 3;
  const int n4 = moffset + 4;
  const int n5 = moffset + 5;
  const int n6 = moffset + 6;
  const int n7 = moffset + 7;
  const int n8 = moffset + 8;
  const int n9 = moffset + 9;
  const int n10= moffset + 10;
  // project onto left eigenvectors
  for (int k=1; k<=W.getsize(2); k++)
  {
    W.set(n1,k,
        rl[1][m1]*Q.get(n1,k)
      + rl[1][m2]*Q.get(n2,k)
      + rl[1][m5]*Q.get(n5,k));
    W.set(n2,k,
        rl[2][m1]*Q.get(n1,k)
      + rl[2][m2]*Q.get(n2,k)
      + rl[2][m3]*Q.get(n3,k)
      + rl[2][m5]*Q.get(n5,k)
      + rl[2][m6]*Q.get(n6,k));
    W.set(n3,k,
        rl[3][m1]*Q.get(n1,k)
      + rl[3][m2]*Q.get(n2,k)
      + rl[3][m4]*Q.get(n4,k)
      + rl[3][m5]*Q.get(n5,k)
      + rl[3][m7]*Q.get(n7,k));
    W.set(n4,k,
        rl[4][m1]*Q.get(n1,k)
      + rl[4][m2]*Q.get(n2,k)
      + rl[4][m5]*Q.get(n5,k));
    W.set(n5,k, // this resembles James's 6th wave
        rl[5][m1]*Q.get(n1,k)
      + rl[5][m2]*Q.get(n2,k)
      + rl[5][m3]*Q.get(n3,k)
      + rl[5][m4]*Q.get(n4,k)
      + rl[5][m5]*Q.get(n5,k)
      + rl[5][m6]*Q.get(n6,k)
      + rl[5][m7]*Q.get(n7,k)
      +           Q.get(n8,k));
    W.set(n6,k, // this resembles James's 5th wave
        rl[6][m1]*Q.get(n1,k)
      + rl[6][m2]*Q.get(n2,k)
      + rl[6][m3]*Q.get(n3,k)
      + rl[6][m5]*Q.get(n5,k)
      + rl[6][m6]*Q.get(n6,k)
      +           Q.get(n9,k));
    W.set(n7,k,
      + rl[7][m1]*Q.get(n1,k)
      + rl[7][m2]*Q.get(n2,k)
      + rl[7][m4]*Q.get(n4,k)
      + rl[7][m5]*Q.get(n5,k)
      + rl[7][m7]*Q.get(n7,k)
      +           Q.get(n10,k));
    W.set(n8,k, // this is James's 9th wave
      + rl[8][m1]*Q.get(n1,k)
      + rl[8][m2]*Q.get(n2,k)
      + rl[8][m4]*Q.get(n4,k)
      + rl[8][m5]*Q.get(n5,k)
      + rl[8][m7]*Q.get(n7,k));
    W.set(n9,k, // this is James's 8th wave
        rl[9][m1]*Q.get(n1,k)
      + rl[9][m2]*Q.get(n2,k)
      + rl[9][m3]*Q.get(n3,k)
      + rl[9][m5]*Q.get(n5,k)
      + rl[9][m6]*Q.get(n6,k));
    W.set(n10,k,
        rl[10][m1]*Q.get(n1,k)
      + rl[10][m2]*Q.get(n2,k)
      + rl[10][m5]*Q.get(n5,k));
  }
}
