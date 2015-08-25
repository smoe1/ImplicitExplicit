#ifndef gas_h
#define gas_h

class dTensor1;
void advection_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *min_speed,
  double *max_speed);

#endif
