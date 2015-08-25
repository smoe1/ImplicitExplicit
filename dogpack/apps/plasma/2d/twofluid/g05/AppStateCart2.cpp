#include "debug.h"
#include "dog_math.h"
#include "tensors.h"
#include "AppStateCart2.h"
#include "AppSolverCart2.h"
#include "DogSolver.h"
#include "Positivity.h"
#include "Limiters.h"
#include "PlasmaParams.h"
#include "Components.h"
#include "float.h"
#include "stdlib.h"
#include "assert.h" // for debugging

double enforce_positivity_of_cell_averages(dTensorBC4& q)
{
  double pos0i[4];
  bool okay_i = true; //=gas10_enforce_positivity_of_cell_averages(_rho_i-1,q, pos0i);
  if(!okay_i && debug2)
  {
    println("ion positivity was:");
    printvn(pos0i[0])
    printvn(pos0i[1])
    printvn(pos0i[2])
    printvn(pos0i[3])
  }
  double pos0e[4];
  bool okay_e = true; //gas10_enforce_positivity_of_cell_averages(_rho_e-1,q, pos0e);
  if(!okay_e && debug2)
  {
    println("elc positivity was:");
    printvn(pos0e[0])
    printvn(pos0e[1])
    printvn(pos0e[2])
    printvn(pos0e[3])
  }

  // rescale positivity indicators to make them comparable
  pos0e[0]*=plasmaParams.get_mass_ratio();
  // should also rescale the electron pressure indicators
  // using the temperature ratio.

  // return minimum overall positivity indicator
  double retval = DBL_MAX;
  for(int k=0;k<4;k++) retval = Min(retval,pos0i[k]);
  for(int k=0;k<4;k++) retval = Min(retval,pos0e[k]);
  return retval;
}

double verify_positivity_of_cell_averages(const dTensorBC4& q)
{
  double pos0i[4];
  bool okay_i = true; //gas10_verify_positivity_of_cell_averages(_rho_i-1,q, pos0i);
  if(!okay_i && debug2)
  {
    println("ion positivity failed:");
    printvn(pos0i[0])
    printvn(pos0i[1])
    printvn(pos0i[2])
    printvn(pos0i[3])
  }
  double pos0e[4];
  bool okay_e = true; //gas10_verify_positivity_of_cell_averages(_rho_e-1,q, pos0e);
  if(!okay_e && debug2)
  {
    println("elc positivity failed:");
    printvn(pos0e[0])
    printvn(pos0e[1])
    printvn(pos0e[2])
    printvn(pos0e[3])
  }

  // rescale positivity indicators to make them comparable
  pos0e[0]*=plasmaParams.get_mass_ratio();
  // should also rescale the electron pressure indicators
  // using the temperature ratio.

  // return minimum overall positivity indicator
  double retval = DBL_MAX;
  for(int k=0;k<4;k++) retval = Min(retval,pos0i[k]);
  for(int k=0;k<4;k++) retval = Min(retval,pos0e[k]);
  return retval;
}

// virtual function that is called after each time stage
//
bool AppStateCart2::AfterAdvanceStage()
{
  // check that positivity of the cell average still holds.
  // if not throw/register an exception that
  // causes the time step to be redone with stricter
  // positivity enforcement (as indicated by some
  // sort of positivity flags).
  double dt = get_solver().get_dt();
  //double minval = verify_positivity_of_cell_averages(get_q());
  double minval = enforce_positivity_of_cell_averages(fetch_q());
  // this updates the positivity flags
  //DogSolver::fetch_solver().fetch_positivityState().update(minval,dt);
  if(minval >= 0.) return true;
  return false;
}

void AppStateCart2::BeforeStage(double dt)
{
}

