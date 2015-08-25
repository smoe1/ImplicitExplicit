#ifndef _gas20_h_
#define _gas20_h_
#include "gas10.h"
class dTensor1;
class dTensor2;

// See also TenMomentComponentID
namespace TwentyMomentComponentID
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
    _Q111=11,
    _Q112=12,
    _Q113=13,
    _Q122=14,
    _Q123=15,
    _Q133=16,
    _Q222=17,
    _Q223=18,
    _Q233=19,
    _Q333=20,

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

void gas20_convert_cons_to_prim(int moffset,
  const dTensor2& cons, dTensor2& prim);
void gas20_convert_prim_to_cons(int moffset,
  const dTensor2& prim, dTensor2& cons);
void gas20_convert_cons_to_prim(int moffset,
  const dTensor1& cons, dTensor1& prim);
void gas20_convert_prim_to_cons(int moffset,
  const dTensor1& prim, dTensor1& cons);

int gas20_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *gas_min_speed,
  double *gas_max_speed);

#endif
