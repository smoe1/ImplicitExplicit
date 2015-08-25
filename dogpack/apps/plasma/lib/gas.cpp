//#include "dog_math.h"
#include "PlasmaParams.h"
#include "debug.h"
#include "DogParams.h"
#include "DogState.h"
#include "tensors1d.h"
#include "gas05.h"
#include "dog_math.h"

class dTensor1;
// void gas_get_extreme_wave_spds(
//   const dTensor1& nhat,
//   const dTensor1& Ql,
//   const dTensor1& Qr,
//   int subsystem_firstindex,
//   int subsystem_size,
//   double *gas_min_speed,
//   double *gas_max_speed)
// {
//   const int moffset = subsystem_firstindex-1;
//   switch(subsystem_size)
//   {
//     case 10:
//       void gas10_get_extreme_wave_spds(
//         const dTensor1& nhat,
//         const dTensor1& Ql,
//         const dTensor1& Qr,
//         int moffset,
//         double *gas_min_speed,
//         double *gas_max_speed);
//       gas10_get_extreme_wave_spds(nhat, Ql, Qr, moffset,
//         gas_min_speed, gas_max_speed);
//       break;
//     case 5:
//       void gas05_get_extreme_wave_spds(
//         const dTensor1& nhat,
//         const dTensor1& Ql,
//         const dTensor1& Qr,
//         int moffset,
//         double *gas_min_speed,
//         double *gas_max_speed);
//       gas05_get_extreme_wave_spds(nhat, Ql, Qr, moffset,
//         gas_min_speed, gas_max_speed);
//       break;
//     default:
//       invalid_value_error(subsystem_size);
//   }
// }

void advection_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *min_speed,
  double *max_speed)
{
  const double& n1 = nhat.get(1);
  const double& n2 = nhat.get(2);

  const int _rho = 1;
  const int _M1 = 2;
  const int _M2 = 3;

  // Left states
  //
  const double rhol    = Ql.get(moffset+_rho);
  const double u1l     = Ql.get(moffset+_M1)/rhol;
  const double u2l     = Ql.get(moffset+_M2)/rhol;
  
  // Right states
  //
  const double rhor    = Qr.get(moffset+_rho);
  const double u1r     = Qr.get(moffset+_M1)/rhor;
  const double u2r     = Qr.get(moffset+_M2)/rhor;
      
  // Average states
  //
  const double u1      = 0.5e0*(u1l+u1r);
  const double u2      = 0.5e0*(u2l+u2r);
  
  // normal velocities
  //
  const double un  = n1*u1  + n2*u2;
  const double unl = n1*u1l + n2*u2l;
  const double unr = n1*u1r + n2*u2r;

  // Minimum speed
  (*min_speed) = Min(unl, un);
  
  // Maximum speed
  (*max_speed) = Max(unr, un);
}

void report_extreme_speeds(double *speeds1, double *speeds2)
{
  int first_maxwell_subsystem_idx=0;
  using namespace PlasmaModelType;
  switch(plasmaParams.get_model())
  {
    case g20:
    case g10:
    case g05:
    case i10e5:
      first_maxwell_subsystem_idx = 3;
      break;
    case p20:
    case p10:
    case p05:
      first_maxwell_subsystem_idx = 2;
      break;
    case mhd:
      break;
    default:
      unsupported_value_error(plasmaParams.get_model());
  }
  for(int i=1;i<first_maxwell_subsystem_idx;i++)
  {
    const double cs_light = maxwellParams.get_cs_light();
    if(speeds1[i]<-cs_light||speeds2[i]>cs_light)
    {
      switch(plasmaParams.get_allow_fast_sound())
      {
        case -1:
          // should change this to throw an exception
          // so that the current time and cell can be printed out, etc.
          errmsg_printf("in subsystem %d speed of light=%f is less than "
            "magnitude of gas speed: \n\t"
            "gas_min_speed=%f, "
            "gas_max_speed=%f",
            i, cs_light,
            speeds1[i],
            speeds2[i]);
        case 0:
          Wprintf("in subsystem %d speed of light=%f is less than "
            "magnitude of gas speed: \n\t"
            "gas_min_speed=%f, "
            "gas_max_speed=%f",
            i, cs_light,
            speeds1[i],
            speeds2[i]);
      }
    }
  }
}

