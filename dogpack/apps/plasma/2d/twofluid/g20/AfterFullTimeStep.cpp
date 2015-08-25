#include "debug.h"
#include "tensors.h"
#include "dog_math.h"
#include "gas20.h"
#include "Components.h"
#include "gasCart2.h"
#include "AppSolverCart2.h"
#include "Relax.h"
#include "Source2Fluid.h"
#include "Positivity.h"
#include "Limiters.h"
#include "L2Project.h"

static void relax_plasma(dTensor2& prim, const double d_time)
{
  SourceSolver::relax_gas(prim, d_time, GasModelType::gas20,
    SpeciesType::ION, _rho_i-1, _B1,_B2,_B3);
  SourceSolver::relax_gas(prim, d_time, GasModelType::gas20,
    SpeciesType::ELC, _rho_e-1, _B1,_B2,_B3);
}

// cons and prim may reference the same array
static void convert_cons_to_prim(const dTensor2& cons, dTensor2& prim)
{
  prim.copyfrom(cons); // no-op if they reference the same thing
  gas20_convert_cons_to_prim(0, cons, prim);
  gas20_convert_cons_to_prim(20, cons, prim);
}

// cons and prim may reference the same array
static void convert_prim_to_cons(const dTensor2& prim, dTensor2& cons)
{
  cons.copyfrom(prim); // no-op if they reference the same thing
  gas20_convert_prim_to_cons(0, prim, cons);
  gas20_convert_prim_to_cons(20, prim, cons);
}

static int check_positivity_cons(const dTensor1& Q)
{
  int gas10_check_positivity_cons(int moffset, const dTensor1& Q);
  int ret = gas10_check_positivity_cons(0, Q);
  ret = iMin(ret,gas10_check_positivity_cons(20, Q));
  return ret;
}

// assume that source term positivity limiters have already been applied
static void solve_source_term(
   const dTensor2& xpts,
   const dTensor2& qvals,
   const dTensor2& auxvals,
   dTensor2& q_new, void* data)
{
   q_new.copyfrom(qvals);
   dTensor2& prim = q_new;
   //double d_time = DogSolver:get_solver().get_dt(); //dt*beta;
   double d_time = *(double*)data;
   convert_cons_to_prim(q_new, prim);
   // the following three operations commute
   assert(plasmaParams.get_source_type_Emom() != SourceType::SPLIT);
   // solve_electro_momentum(prim, d_time); // one
   assert(plasmaParams.get_source_type_rotP() != SourceType::SPLIT);
   // rotate_pressure_tensors(prim, d_time); // three
   if(plasmaParams.get_source_type_isoP() == SourceType::SPLIT)
   {
     relax_plasma(prim, d_time); // two
   }
   // this conversion will automatically update the kinetic
   // energies based on the new velocities and pressures
   // (so electromomentum evolution effectively evolves KE)
   convert_prim_to_cons(prim, q_new);
}

// to satisfy the compiler
void AfterFullTimeStep(DogSolverCart2& solver)
{
  //eprintf("you forgot to override me!");
}

// advance the operator that has been split off from dogpack framework
//
// nonzero return value would indicate error. The exact solution
// has no stability constraint, so this always returns 0.
//
int AppStateCart2::advanceSplitTimeStep(double dt)
{
  bool using_split_source = 
     plasmaParams.get_source_type_Emom() == SourceType::SPLIT
  || plasmaParams.get_source_type_isoP() == SourceType::SPLIT
  || plasmaParams.get_source_type_rotP() == SourceType::SPLIT;

  if(!using_split_source) return 0;

  dTensorBC4& q = fetch_q();
  dTensorBC4& aux = fetch_aux();
  const int mxelems = q.getsize(1);
  const int myelems = q.getsize(2);
  const int    mbc = q.getmbc();
  // apply positivity limiters to source term quadrature points
  PositivityLimiter limiter(q,PositivityPointType::SOURCE);
  limiter.apply();
  //dprintf("solving source term for dt = %e",dt);
  L2Project_extra(1,mxelems,1,myelems,
      q,aux,q,&solve_source_term,&dt);
  return 0;
}
