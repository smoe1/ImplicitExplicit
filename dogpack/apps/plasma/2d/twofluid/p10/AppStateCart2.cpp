#include "debug.h"
#include "assert.h"
#include "dog_math.h"
#include "gas10.h" // for ini_doc
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Source2Fluid.h"
#include "Components.h"
#include "App.h"
#include "AppStateCart2.h"
#include "Limiters.h"
#include "L2Project.h"
#include "Params.h"
#include "float.h"
#include <string>

// ******* section: AppStateCart2 *******

static void solve_electro_momentum(dTensor2& prim, const double d_time)
{
  SourceSolver::solve_electro_momentum_sympair(
    prim, d_time,
    _rho_i,
    _B1, _B2, _E3);
}

static void isotropize_pressure_tensors(dTensor2& prim, const double d_time)
{
  SourceSolver::isotropize_pressure(prim, d_time,
    SpeciesType::ION, _rho_i-1, _B1,_B2,0);
}

static void rotate_pressure_tensors(dTensor2& prim, const double d_time)
{
  const double ion_q_over_m = plasmaParams.get_ion_q_over_m();
  SourceSolver::rotate_pressure_tensor(prim, d_time, ion_q_over_m,
    _rho_i-1, _B1,_B2,0);
}

// cons and prim may reference the same array
static void convert_cons_to_prim(const dTensor2& cons, dTensor2& prim)
{
  prim.copyfrom(cons); // no-op if they reference the same thing
  gas10_convert_cons_to_prim(0, cons, prim);
}

static void print_prim(const dTensor2& prim)
{
  printf("\n ion pressures:\n");
  gas10_print_prim(0, prim);
}

// cons and prim may reference the same array
static void convert_prim_to_cons(const dTensor2& prim, dTensor2& cons)
{
  cons.copyfrom(prim); // no-op if they reference the same thing
  gas10_convert_prim_to_cons(0, prim, cons);
}

static int check_positivity_cons(const dTensor1& Q)
{
  int gas10_check_positivity_cons(int moffset, const dTensor1& Q);
  int ret = gas10_check_positivity_cons(0, Q);
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

void AppStateCart2::ReportAfterStep()const
{
    const dTensorBC4& aux = get_aux();
    const dTensorBC4& q = get_q();
    const double t = get_time();

    assert(dogParams.get_mcapa()<1); // without capacity function
    void WriteConservation(const dTensorBC4& q, double t);
    WriteConservation(q, t);

    if(Params::get_problem()==Params::ODE)
        return;

    double OutputMagneticFluxes(
        int n_B1,
        int n_B2,
        const dTensorBC4& q,
        double t);
    const double bottom_to_rght_flux = OutputMagneticFluxes(_B1, _B2, q, t);

    void Output_center_E3(const dTensorBC4& q, double t, int _E3);
    Output_center_E3(q, t, _E3);
    void Output_center_gas_states(
      const dTensorBC4& q, double t, const char* outputfile,
      int species_offset, bool fiveMoment);
    const char* get_outputdir();
    string outputfile_i = string(get_outputdir())+"/xpoint_gas_i.dat";
    Output_center_gas_states(q, t, outputfile_i.c_str(), _rho_i-1,false);
    void Output_center_cell_state(const dTensorBC4& q, double t);
    Output_center_cell_state(q, t);

    // output the rate of production of entropy at the center
    //
    const int maux = aux.getsize(3);
    if(maux>0)
    {
      assert_eq(maux,2);
      void output_component_at_center(
        const dTensorBC4& q, double t, const char* outputfile,
        int idx);
      outputfile_i = string(get_outputdir())+"/xpoint_ent_prod_i.dat";
      output_component_at_center(aux, t, outputfile_i.c_str(),1);
    }

    void output_local_divergence_error( int method_order, double dx, double dy,
        int n_B1, int n_B2, const dTensorBC4& q, double t);
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    output_local_divergence_error(
      dogParams.get_space_order(), dx, dy, _B1, _B2, q, t);

    static bool unit_flux_triggered = bottom_to_rght_flux >=1;
    if(!unit_flux_triggered && bottom_to_rght_flux >= 1.)
    {
      unit_flux_triggered = true;
      dprintf("writing restart file q8000.dat at time %f;"
        "\n\tbottom_to_rght_flux = %f", t, bottom_to_rght_flux);
      get_solver().write_restart(8000);
    }
}

static double enforce_positivity_of_cell_averages(dTensorBC4& q)
{
  double pos0i[4];
  bool okay_i = gas10_enforce_positivity_of_cell_averages(_rho_i-1,q, pos0i);
  if(!okay_i && debug2)
  {
    println("ion positivity was:");
    printvn(pos0i[0])
    printvn(pos0i[1])
    printvn(pos0i[2])
    printvn(pos0i[3])
  }

  // return minimum overall positivity indicator
  double retval = DBL_MAX;
  for(int k=0;k<4;k++) retval = Min(retval,pos0i[k]);
  return retval;
}

static double verify_positivity_of_cell_averages(const dTensorBC4& q)
{
  double pos0i[4];
  bool okay_i = gas10_verify_positivity_of_cell_averages(_rho_i-1,q, pos0i);
  if(!okay_i && debug2)
  {
    println("ion positivity failed:");
    printvn(pos0i[0])
    printvn(pos0i[1])
    printvn(pos0i[2])
    printvn(pos0i[3])
  }

  // return minimum overall positivity indicator
  double retval = DBL_MAX;
  for(int k=0;k<4;k++) retval = Min(retval,pos0i[k]);
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

