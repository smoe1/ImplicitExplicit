#include "dog_math.h"
#include "PlasmaParams.h"
#include "debug.h"
#include "assert.h"
#include "gas.h"
#include "gas20.h"
class dTensor1;
using namespace PlasmaModelType;

void SetWaveSpds(
  const dTensor1& nhat,
  const dTensor1& xedge,
  const dTensor1& Ql,
  const dTensor1& Qr,
  const dTensor1& Auxl,
  const dTensor1& Auxr,
  double* speeds1,double* speeds2)
{
  // for(int i=1;i<=plasmaParams.get_num_nonlinear_subsystems();i++)
  // {
  //   gas_get_extreme_wave_spds(nhat, Ql, Qr, 
  //     plasmaParams.get_subsystem_firstindex(i),
  //     plasmaParams.get_subsystem_size(i),
  //     &speeds1[i], &speeds2[i]);
  // }
  switch(plasmaParams.get_model())
  {
    case g05:
      gas05_get_extreme_wave_spds(nhat,Ql,Qr, 5, &speeds1[2],&speeds2[2]);
    case p05:
      gas05_get_extreme_wave_spds(nhat,Ql,Qr, 0, &speeds1[1],&speeds2[1]);
      break;
    case g10:
      gas10_get_extreme_wave_spds(nhat,Ql,Qr,10, &speeds1[2],&speeds2[2]);
    case p10:
      gas10_get_extreme_wave_spds(nhat,Ql,Qr, 0, &speeds1[1],&speeds2[1]);
      break;
    case g20:
      gas20_get_extreme_wave_spds(nhat,Ql,Qr,20, &speeds1[2],&speeds2[2]);
    case p20:
      gas20_get_extreme_wave_spds(nhat,Ql,Qr, 0, &speeds1[1],&speeds2[1]);
      break;
    case i10e5:
      gas10_get_extreme_wave_spds(nhat,Ql,Qr, 0, &speeds1[1],&speeds2[1]);
      gas05_get_extreme_wave_spds(nhat,Ql,Qr,10, &speeds1[2],&speeds2[2]);
      break;
    case mhd:
    default:
      unsupported_value_error(plasmaParams.get_model());
  }
  //
  // set wave speeds for maxwell system
  switch(plasmaParams.get_model())
  {
    case g20:
    case g10:
    case g05:
    case i10e5:
      switch(plasmaParams.get_num_riemann_subsystems())
      {
       default:
        invalid_value_error(plasmaParams.get_num_riemann_subsystems());
       case 7:
        advection_get_extreme_wave_spds(nhat, Ql, Qr, 
          plasmaParams.get_riemann_subsystem_firstindex(2)-1,
          &speeds1[7], &speeds2[7]);
        advection_get_extreme_wave_spds(nhat, Ql, Qr, 0,
          &speeds1[6], &speeds2[6]);
       case 5:
        speeds1[5] = plasmaParams.get_min_speed();
        speeds2[5] = plasmaParams.get_max_speed();
       case 4:
        speeds1[4] = plasmaParams.get_min_speed();
        speeds2[4] = plasmaParams.get_max_speed();
       case 3:
        speeds1[3] = -maxwellParams.get_cs_light();
        speeds2[3] = maxwellParams.get_cs_light();
      }
      break;
    case p20:
    case p10:
    case p05:
      switch(plasmaParams.get_num_riemann_subsystems())
      {
       default:
        invalid_value_error(plasmaParams.get_num_riemann_subsystems());
       case 4:
        // correct
        advection_get_extreme_wave_spds(nhat, Ql, Qr, 0,
          &speeds1[4], &speeds2[4]);
        // for debugging
        //speeds1[4] = plasmaParams.get_min_speed();
        //speeds2[4] = plasmaParams.get_max_speed();
       case 3:
        speeds1[3] = plasmaParams.get_min_speed();
        speeds2[3] = plasmaParams.get_max_speed();
       case 2:
        speeds1[2] = -maxwellParams.get_cs_light();
        speeds2[2] = maxwellParams.get_cs_light();
      }
      break;
    case mhd:
    default:
      unsupported_value_error(plasmaParams.get_model());
  }
  void report_extreme_speeds(double *s1, double *s2);
  report_extreme_speeds(speeds1,speeds2);
}

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
// This method is used if
// lib/2d/RiemannSolve.cpp is used but not if
// plasma/lib/2d/RiemannSolve.cpp is used.
//
void SetWaveSpd(
  const dTensor1& nhat,
  const dTensor1& xedge,
  const dTensor1& Ql,
  const dTensor1& Qr,
  const dTensor1& Auxl,
  const dTensor1& Auxr,
  double& s1,double& s2)
{
  // extreme speeds of the subsystems
  double speeds1[MAX_NUM_SUBSYSTEMS],speeds2[MAX_NUM_SUBSYSTEMS];
  assert_le(plasmaParams.get_num_riemann_subsystems(),MAX_NUM_SUBSYSTEMS);
  // initialize global extreme speeds
  speeds1[0] = 1.;
  speeds2[0] = -1.;
  SetWaveSpds(nhat, xedge, Ql, Qr, Auxl, Auxr, speeds1, speeds2);
  s1 = speeds1[0];
  s2 = speeds2[0];
  // if SetWaveSpds did not set global speeds then
  // set to be the extreme speeds of the subsystems
  if(speeds1[0] > speeds2[0])
  {
    s1 = speeds1[1];
    s2 = speeds2[1];
    for(int i=2; i<= plasmaParams.get_num_riemann_subsystems(); i++)
    {
      s1 = Min(s1,speeds1[i]);
      s2 = Max(s2,speeds2[i]);
    }
  }
}

