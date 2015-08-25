#include "AppSolver.h"
#include "Source2Fluid.h"
#include "Components.h"
#include "gas05.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Limiters.h"
#include "L2Project.h"
#include "IniDocument.h" // for ini_doc
#include "float.h" // for DBL_MAX
#include "dog_math.h" // for Min

// ******* section: AppSolver *******

void AppSolverCart2::init()
{
  if(!&get_state())
  {
    set_state(new AppStateCart2);
    fetch_state().init();
  }

  DogSolverCart2::init();
}

void AppSolverCart2::initParams()
{
    DogSolverCart2::initParams();
    InitApp(); // get application parameters
}

// ******* section: AppState *******

void solve_electro_momentum(dTensor2& prim, const double d_time)
{
  SourceSolver::solve_electro_momentum_sympair(
    prim, d_time,
    _rho_i,
    _B1, _B2, _E3);
}

// cons and prim may reference the same array
void convert_cons_to_prim(const dTensor2& cons, dTensor2& prim)
{
  prim.copyfrom(cons); // no-op if they reference the same thing
  gas05_convert_cons_to_prim(0, cons, prim);
}

// cons and prim may reference the same array
void convert_prim_to_cons(const dTensor2& prim, dTensor2& cons)
{
  cons.copyfrom(prim); // no-op if they reference the same thing
  gas05_convert_prim_to_cons(0, prim, cons);
}

int check_positivity_cons(const dTensor1& Q)
{
  int gas05_check_positivity_cons(int moffset, const dTensor1& Q);
  int ret = gas05_check_positivity_cons(0, Q);
  return ret;
}

// assume that source term positivity limiters have already been applied
void solve_source_term(
   const dTensor2& xpts,
   const dTensor2& qvals,
   const dTensor2& auxvals,
   dTensor2& q_new, void* data)
{
   q_new.copyfrom(qvals);
   dTensor2& prim = q_new;
   //double d_time = DogSolver:get_solver().get_dt(); //dt*beta;
   double d_time = *(double*)data;
   // it would not be necessary to convert variables
   // if I had written solve_electro_momentum
   // to work with the momentum instead of the velocities
   convert_cons_to_prim(q_new, prim);
   if(plasmaParams.get_source_type_Emom() == SourceType::SPLIT)
     solve_electro_momentum(prim, d_time); // one
   convert_prim_to_cons(prim, q_new);
}

// to satisfy the compiler
//void AfterFullTimeStep(DogSolverCart2& solver)
//{
//  //eprintf("you forgot to override me!");
//}

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
    Output_center_gas_states(q, t, outputfile_i.c_str(), _rho_i-1,true);
    void Output_center_cell_state(const dTensorBC4& q, double t);
    Output_center_cell_state(q, t);

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
  // need to create version for five-moment gas
  bool okay_i = true; // = gas10_enforce_positivity_of_cell_averages(_rho_i-1,q, pos0i);
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

