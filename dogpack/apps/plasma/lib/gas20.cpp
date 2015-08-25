//#define CHECK_BOUNDS
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "PlasmaParams.h"
#include "gas20.h"
//#include "Limiters.h" // for PositivityReturnType

// ten-moment gas routines

// cons and prim may reference the same array
void gas20_convert_cons_to_prim(int moffset,
  const dTensor1& cons, dTensor1& prim)
{
  using namespace TwentyMomentComponentID;
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
  const int Q111_ = moffset+_Q111;
  const int Q112_ = moffset+_Q112;
  const int Q113_ = moffset+_Q113;
  const int Q122_ = moffset+_Q122;
  const int Q123_ = moffset+_Q123;
  const int Q133_ = moffset+_Q133;
  const int Q222_ = moffset+_Q222;
  const int Q223_ = moffset+_Q223;
  const int Q233_ = moffset+_Q233;
  const int Q333_ = moffset+_Q333;

  const double rho = cons.get(rho_);
  const double rho_inv = 1./rho;
  const double M1  = cons.get(M1_);
  const double M2  = cons.get(M2_);
  const double M3  = cons.get(M3_);
  const double u1  = M1*rho_inv;
  const double u2  = M2*rho_inv;
  const double u3  = M3*rho_inv;
  //
  // the energy tensor
  //
  const double N11 = cons.get(P11_);
  const double N12 = cons.get(P12_);
  const double N13 = cons.get(P13_);
  const double N22 = cons.get(P22_);
  const double N23 = cons.get(P23_);
  const double N33 = cons.get(P33_);
  //
  const double QN111 = cons.get(Q111_);
  const double QN112 = cons.get(Q112_);
  const double QN113 = cons.get(Q113_);
  const double QN122 = cons.get(Q122_);
  const double QN123 = cons.get(Q123_);
  const double QN133 = cons.get(Q133_);
  const double QN222 = cons.get(Q222_);
  const double QN223 = cons.get(Q223_);
  const double QN233 = cons.get(Q233_);
  const double QN333 = cons.get(Q333_);
  //
  // the kinetic energy tensor
  //
  const double Mu11 = M1*u1;
  const double Mu12 = M1*u2;
  const double Mu13 = M1*u3;
  const double Mu22 = M2*u2;
  const double Mu23 = M2*u3;
  const double Mu33 = M3*u3;
  //
  const double Sym3uN111 = u1*N11*3;
  const double Sym3uN112 = u1*N12*2 + u2*N11;
  const double Sym3uN113 = u1*N13*2 + u3*N11;
  const double Sym3uN122 = u1*N22 + 2*u2*N12;
  const double Sym3uN123 = u1*N23 + u2*N13 + u3*N12;
  const double Sym3uN133 = u1*N33 + 2*u3*N13;
  const double Sym3uN222 = u2*N22*3;
  const double Sym3uN223 = u2*N23*2 + u3*N22;
  const double Sym3uN233 = u2*N33 + 2*u3*N23;
  const double Sym3uN333 = u3*N33*3;

  prim.set(rho_, rho);
  prim.set(M1_, u1);
  prim.set(M2_, u2);
  prim.set(M3_, u3);
  prim.set(P11_, N11 - Mu11);
  prim.set(P12_, N12 - Mu12);
  prim.set(P13_, N13 - Mu13);
  prim.set(P22_, N22 - Mu22);
  prim.set(P23_, N23 - Mu23);
  prim.set(P33_, N33 - Mu33);
  prim.set(Q111_, QN111 - Sym3uN111 + 2*u1*Mu11);
  prim.set(Q112_, QN112 - Sym3uN112 + 2*u1*Mu12);
  prim.set(Q113_, QN113 - Sym3uN113 + 2*u1*Mu13);
  prim.set(Q122_, QN122 - Sym3uN122 + 2*u1*Mu22);
  prim.set(Q123_, QN123 - Sym3uN123 + 2*u1*Mu23);
  prim.set(Q133_, QN133 - Sym3uN133 + 2*u1*Mu33);
  prim.set(Q222_, QN222 - Sym3uN222 + 2*u2*Mu22);
  prim.set(Q223_, QN223 - Sym3uN223 + 2*u2*Mu23);
  prim.set(Q233_, QN233 - Sym3uN233 + 2*u2*Mu33);
  prim.set(Q333_, QN333 - Sym3uN333 + 2*u3*Mu33);
}

// cons and prim may reference the same array
void gas20_convert_cons_to_prim(int moffset,
  const dTensor2& cons, dTensor2& prim)
{
  using namespace TwentyMomentComponentID;
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
  const int Q111_ = moffset+_Q111;
  const int Q112_ = moffset+_Q112;
  const int Q113_ = moffset+_Q113;
  const int Q122_ = moffset+_Q122;
  const int Q123_ = moffset+_Q123;
  const int Q133_ = moffset+_Q133;
  const int Q222_ = moffset+_Q222;
  const int Q223_ = moffset+_Q223;
  const int Q233_ = moffset+_Q233;
  const int Q333_ = moffset+_Q333;

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
    const double QN111 = cons.get(i, Q111_);
    const double QN112 = cons.get(i, Q112_);
    const double QN113 = cons.get(i, Q113_);
    const double QN122 = cons.get(i, Q122_);
    const double QN123 = cons.get(i, Q123_);
    const double QN133 = cons.get(i, Q133_);
    const double QN222 = cons.get(i, Q222_);
    const double QN223 = cons.get(i, Q223_);
    const double QN233 = cons.get(i, Q233_);
    const double QN333 = cons.get(i, Q333_);
    const double Mu11 = M1*u1;
    const double Mu12 = M1*u2;
    const double Mu13 = M1*u3;
    const double Mu22 = M2*u2;
    const double Mu23 = M2*u3;
    const double Mu33 = M3*u3;
    const double Sym3uN111 = u1*N11*3;
    const double Sym3uN112 = u1*N12*2 + u2*N11;
    const double Sym3uN113 = u1*N13*2 + u3*N11;
    const double Sym3uN122 = u1*N22 + 2*u2*N12;
    const double Sym3uN123 = u1*N23 + u2*N13 + u3*N12;
    const double Sym3uN133 = u1*N33 + 2*u3*N13;
    const double Sym3uN222 = u2*N22*3;
    const double Sym3uN223 = u2*N23*2 + u3*N22;
    const double Sym3uN233 = u2*N33 + 2*u3*N23;
    const double Sym3uN333 = u3*N33*3;

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
    prim.set(i, Q111_, QN111 - Sym3uN111 + 2*u1*Mu11);
    prim.set(i, Q112_, QN112 - Sym3uN112 + 2*u1*Mu12);
    prim.set(i, Q113_, QN113 - Sym3uN113 + 2*u1*Mu13);
    prim.set(i, Q122_, QN122 - Sym3uN122 + 2*u1*Mu22);
    prim.set(i, Q123_, QN123 - Sym3uN123 + 2*u1*Mu23);
    prim.set(i, Q133_, QN133 - Sym3uN133 + 2*u1*Mu33);
    prim.set(i, Q222_, QN222 - Sym3uN222 + 2*u2*Mu22);
    prim.set(i, Q223_, QN223 - Sym3uN223 + 2*u2*Mu23);
    prim.set(i, Q233_, QN233 - Sym3uN233 + 2*u2*Mu33);
    prim.set(i, Q333_, QN333 - Sym3uN333 + 2*u3*Mu33);
  }
}

// prim and cons may reference the same array
void gas20_convert_prim_to_cons(int moffset,
  const dTensor1& prim, dTensor1& cons)
{
  using namespace TwentyMomentComponentID;
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
  const int Q111_ = moffset+_Q111;
  const int Q112_ = moffset+_Q112;
  const int Q113_ = moffset+_Q113;
  const int Q122_ = moffset+_Q122;
  const int Q123_ = moffset+_Q123;
  const int Q133_ = moffset+_Q133;
  const int Q222_ = moffset+_Q222;
  const int Q223_ = moffset+_Q223;
  const int Q233_ = moffset+_Q233;
  const int Q333_ = moffset+_Q333;

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
  const double Q111 = prim.get(Q111_);
  const double Q112 = prim.get(Q112_);
  const double Q113 = prim.get(Q113_);
  const double Q122 = prim.get(Q122_);
  const double Q123 = prim.get(Q123_);
  const double Q133 = prim.get(Q133_);
  const double Q222 = prim.get(Q222_);
  const double Q223 = prim.get(Q223_);
  const double Q233 = prim.get(Q233_);
  const double Q333 = prim.get(Q333_);
  const double Mu11 = M1*u1;
  const double Mu12 = M1*u2;
  const double Mu13 = M1*u3;
  const double Mu22 = M2*u2;
  const double Mu23 = M2*u3;
  const double Mu33 = M3*u3;
  const double Sym3uP111 = u1*P11*3;
  const double Sym3uP112 = u1*P12*2 + u2*P11;
  const double Sym3uP113 = u1*P13*2 + u3*P11;
  const double Sym3uP122 = u1*P22 + 2*u2*P12;
  const double Sym3uP123 = u1*P23 + u2*P13 + u3*P12;
  const double Sym3uP133 = u1*P33 + 2*u3*P13;
  const double Sym3uP222 = u2*P22*3;
  const double Sym3uP223 = u2*P23*2 + u3*P22;
  const double Sym3uP233 = u2*P33 + 2*u3*P23;
  const double Sym3uP333 = u3*P33*3;
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
  cons.set(Q111_, Q111 + Sym3uP111 + u1*Mu11);
  cons.set(Q112_, Q112 + Sym3uP112 + u1*Mu12);
  cons.set(Q113_, Q113 + Sym3uP113 + u1*Mu13);
  cons.set(Q122_, Q122 + Sym3uP122 + u1*Mu22);
  cons.set(Q123_, Q123 + Sym3uP123 + u1*Mu23);
  cons.set(Q133_, Q133 + Sym3uP133 + u1*Mu33);
  cons.set(Q222_, Q222 + Sym3uP222 + u2*Mu22);
  cons.set(Q223_, Q223 + Sym3uP223 + u2*Mu23);
  cons.set(Q233_, Q233 + Sym3uP233 + u2*Mu33);
  cons.set(Q333_, Q333 + Sym3uP333 + u3*Mu33);
}

// prim and cons may reference the same array
void gas20_convert_prim_to_cons(int moffset,
  const dTensor2& prim, dTensor2& cons)
{
  using namespace TwentyMomentComponentID;
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
  const int Q111_ = moffset+_Q111;
  const int Q112_ = moffset+_Q112;
  const int Q113_ = moffset+_Q113;
  const int Q122_ = moffset+_Q122;
  const int Q123_ = moffset+_Q123;
  const int Q133_ = moffset+_Q133;
  const int Q222_ = moffset+_Q222;
  const int Q223_ = moffset+_Q223;
  const int Q233_ = moffset+_Q233;
  const int Q333_ = moffset+_Q333;

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
    const double Q111 = prim.get(i, Q111_);
    const double Q112 = prim.get(i, Q112_);
    const double Q113 = prim.get(i, Q113_);
    const double Q122 = prim.get(i, Q122_);
    const double Q123 = prim.get(i, Q123_);
    const double Q133 = prim.get(i, Q133_);
    const double Q222 = prim.get(i, Q222_);
    const double Q223 = prim.get(i, Q223_);
    const double Q233 = prim.get(i, Q233_);
    const double Q333 = prim.get(i, Q333_);
    const double Mu11 = M1*u1;
    const double Mu12 = M1*u2;
    const double Mu13 = M1*u3;
    const double Mu22 = M2*u2;
    const double Mu23 = M2*u3;
    const double Mu33 = M3*u3;
    const double Sym3uP111 = u1*P11*3;
    const double Sym3uP112 = u1*P12*2 + u2*P11;
    const double Sym3uP113 = u1*P13*2 + u3*P11;
    const double Sym3uP122 = u1*P22 + 2*u2*P12;
    const double Sym3uP123 = u1*P23 + u2*P13 + u3*P12;
    const double Sym3uP133 = u1*P33 + 2*u3*P13;
    const double Sym3uP222 = u2*P22*3;
    const double Sym3uP223 = u2*P23*2 + u3*P22;
    const double Sym3uP233 = u2*P33 + 2*u3*P23;
    const double Sym3uP333 = u3*P33*3;
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
    cons.set(i, Q111_, Q111 + Sym3uP111 + u1*Mu11);
    cons.set(i, Q112_, Q112 + Sym3uP112 + u1*Mu12);
    cons.set(i, Q113_, Q113 + Sym3uP113 + u1*Mu13);
    cons.set(i, Q122_, Q122 + Sym3uP122 + u1*Mu22);
    cons.set(i, Q123_, Q123 + Sym3uP123 + u1*Mu23);
    cons.set(i, Q133_, Q133 + Sym3uP133 + u1*Mu33);
    cons.set(i, Q222_, Q222 + Sym3uP222 + u2*Mu22);
    cons.set(i, Q223_, Q223 + Sym3uP223 + u2*Mu23);
    cons.set(i, Q233_, Q233 + Sym3uP233 + u2*Mu33);
    cons.set(i, Q333_, Q333 + Sym3uP333 + u3*Mu33);
  }
}

// returns zero on success
int gas20_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *gas_min_speed,
  double *gas_max_speed)
{
  double minspd;
  double maxspd;
  gas10_get_extreme_wave_spds(nhat, Ql, Qr, moffset,
    &minspd, &maxspd);

  // 2.85 comes from the fastest wave speed for a Gaussian
  // distribution in the 35-moment model.
  // For the 20-moment model we should be able to replace
  // 2.85 with something between 2.85 and sq3=1.73,
  // but I don't know exactly what to choose. -eaj
  const double fastest_speed_over_slow_sound_speed = 2.85;
  const double add_speed
    = (.5*fastest_speed_over_slow_sound_speed/sq3)*(maxspd-minspd);
  const double avg_speed = (minspd+maxspd)*.5;
  minspd = avg_speed - add_speed;
  maxspd = avg_speed + add_speed;
  (*gas_min_speed) = minspd;
  (*gas_max_speed) = maxspd;
}

