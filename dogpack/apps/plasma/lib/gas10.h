#ifndef _gas10_h_
#define _gas10_h_
#include "gas05.h"
class dTensor1;
class dTensor2;

// See also FiveMomentComponentID
namespace TenMomentComponentID
{
  enum Enum{
    _rho=1,  // rho
    _M1 =2,  // 1-momentum
    _M2 =3,  // 2-momentum
    _M3 =4,  // 3-momentum
    _N11=5,  // energy_11
    _N12=6,  // energy_12
    _N13=7,  // energy_13
    _N22=8,  // energy_22
    _N23=9,  // energy_23
    _N33=10, // energy_33

    _u1  = _M1,
    _u2  = _M2,
    _u3  = _M3,
    _P11 = _N11,
    _P12 = _N12,
    _P13 = _N13,
    _P22 = _N22,
    _P23 = _N23,
    _P33 = _N33,
  };
}

void gas10_get_pressure_evals(double* evals,
  double P11, double P12, double P13, double P22, double P23, double P33);
void gas10_print_prim(int moffset, const dTensor2& prim);
void gas10_convert_cons_to_prim(int moffset,
  const dTensor2& cons, dTensor2& prim);
void gas10_convert_prim_to_cons(int moffset,
  const dTensor2& prim, dTensor2& cons);
void gas10_convert_cons_to_prim(int moffset,
  const dTensor1& cons, dTensor1& prim);
void gas10_convert_prim_to_cons(int moffset,
  const dTensor1& prim, dTensor1& cons);
void gas10_rescale_pressure_anisotropy(dTensor1& prim,
  double p, double rescaling_factor);

// inline makes sense for short functions
// with many arguments
inline double det3x3symm(
  double T11,
  double T12,
  double T13,
  double T22,
  double T23,
  double T33)
{
  // find components of a row of the adjugate
  double TA11 = T22*T33 - T23*T23;
  double TA12 = T23*T13 - T12*T33;
  double TA13 = T12*T23 - T22*T13;
  // find the determinant
  return TA11*T11 + TA12*T12 + TA13*T13;
}

int gas10_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *gas_min_speed,
  double *gas_max_speed);

double get_entropy_per_vol10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33);
void set_entropy_per_vol10(dTensor2& q_new, const dTensor2& qvals,
    int firstIdx, int _entropy_s);
double get_entropy_per_mass10(double rho,double M1,double M2,double M3,
  double N11,double N12,double N13,double N22,double N23,double N33);

#endif
