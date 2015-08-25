#include "debug.h"
#include "tensors.h"
#include "dog_math.h"
#include "gas10.h"
#include "Components.h"
#include "gasCart2.h"
#include "AppSolverCart2.h"
#include "Relax.h"
#include "Source2Fluid.h"
#include "Positivity.h"
#include "Limiters.h"
#include "L2Project.h"

static void reset_entropy(
   const dTensor2& xpts,
   const dTensor2& qvals,
   const dTensor2& auxvals,
   dTensor2& q_new)
{
    // copy q_vals to q_new
    // q_new.copyall(q_vals);
    for(int i=1; i<=qvals.getsize(1); i++)
    for(int k=1;k<=qvals.getsize(2);k++)
    {
        q_new.set(i,k, qvals.get(i,k));
    }
    set_entropy_per_vol10(q_new, qvals, _rho_i, _entropy_i);
    set_entropy_per_vol10(q_new, qvals, _rho_e, _entropy_e);
}

#if 0
static void isotropize(
   const dTensor2& qvals,
   dTensor2& q_new,
   void* data)
{
  const int numpts=qvals.getsize(1);
  double d_time = *(double*)data;

  bool did_ion_relax = Relax::isotropize_spc(
      q_new,
      qvals,
      numpts,
      d_time,
      plasmaParams.get_iso_period_type(),
      plasmaParams.get_ion_base_iso_period(),
      plasmaParams.get_ion_mass(),
      _rho_i,
      _M1_i,
      _M2_i,
      _M3_i,
      _N11_i,
      _N12_i,
      _N13_i,
      _N22_i,
      _N23_i,
      _N33_i,
      _B1,
      _B2,
      _B3);
  bool did_elc_relax = Relax::isotropize_spc(
      q_new,
      qvals,
      numpts,
      d_time,
      plasmaParams.get_iso_period_type(),
      plasmaParams.get_elc_base_iso_period(),
      plasmaParams.get_elc_mass(),
      _rho_e,
      _M1_e,
      _M2_e,
      _M3_e,
      _N11_e,
      _N12_e,
      _N13_e,
      _N22_e,
      _N23_e,
      _N33_e,
      _B1,
      _B2,
      _B3);
}
#endif

static void isotropize_pressure_tensors(dTensor2& prim, const double d_time)
{
  SourceSolver::isotropize_pressure(prim, d_time,
    SpeciesType::ION, _rho_i-1, _B1,_B2,_B3);
  SourceSolver::isotropize_pressure(prim, d_time,
    SpeciesType::ELC, _rho_e-1, _B1,_B2,_B3);
}

static void rotate_pressure_tensors(dTensor2& prim, const double d_time)
{
  const double ion_q_over_m = plasmaParams.get_ion_q_over_m();
  const double elc_q_over_m = plasmaParams.get_elc_q_over_m();
  SourceSolver::rotate_pressure_tensor(prim, d_time, ion_q_over_m,
    _rho_i-1, _B1,_B2,_B3);
  SourceSolver::rotate_pressure_tensor(prim, d_time, elc_q_over_m,
    _rho_e-1, _B1,_B2,_B3);
}

static void solve_electro_momentum(dTensor2& prim, const double d_time)
{
  // pass all the electromagnetic components so that
  // we can use zero index to signal that a component is zero.
  SourceSolver::solve_electro_momentum(
    prim, d_time,
    _rho_i, _rho_e,
    _B1, _B2, _B3,
    _E1, _E2, _E3);
}

// cons and prim may reference the same array
static void convert_cons_to_prim(const dTensor2& cons, dTensor2& prim)
{
  prim.copyfrom(cons); // no-op if they reference the same thing
  gas10_convert_cons_to_prim(0, cons, prim);
  gas10_convert_cons_to_prim(10, cons, prim);
}

static void print_prim(const dTensor2& prim)
{
  printf("\n ion pressures:\n");
  gas10_print_prim(0, prim);
  printf(" elc pressures:\n");
  gas10_print_prim(10, prim);
}

// cons and prim may reference the same array
static void convert_prim_to_cons(const dTensor2& prim, dTensor2& cons)
{
  cons.copyfrom(prim); // no-op if they reference the same thing
  gas10_convert_prim_to_cons(0, prim, cons);
  gas10_convert_prim_to_cons(10, prim, cons);
}

static int check_positivity_cons(const dTensor1& Q)
{
  int gas10_check_positivity_cons(int moffset, const dTensor1& Q);
  int ret = gas10_check_positivity_cons(0, Q);
  ret = iMin(ret,gas10_check_positivity_cons(10, Q));
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
   if(plasmaParams.get_source_type_Emom() == SourceType::SPLIT)
     solve_electro_momentum(prim, d_time); // one
   if(plasmaParams.get_source_type_isoP() == SourceType::SPLIT)
     isotropize_pressure_tensors(prim, d_time); // two
   if(plasmaParams.get_source_type_rotP() == SourceType::SPLIT)
     rotate_pressure_tensors(prim, d_time); // three
   // this conversion will automatically update the kinetic
   // energies based on the new velocities and pressures
   // (so electromomentum evolution effectively evolves KE)
   if(debug3)
   {
     print_prim(prim);
   }
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
