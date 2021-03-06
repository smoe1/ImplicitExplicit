#include <assert.h>
#include "debug.h"
#include "tensors.h"
#include "constants.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Components.h"
#include "AppSolverCart2.h"
#include <string>
#include "GEM.h"

void ConSoln(
    const dTensorBC4& aux,
    const dTensorBC4& q,
    double t)
{
  eprintf("you forgot to override me by defining ReportAfterStep!");
}

void AppSolverCart2::truncate_logfiles(double time)const
{
  // the default will actually work
  DogSolver::truncate_logfiles(time);
  //truncate_logfile(time, "conservation.dat");
  //truncate_logfile(time, "recon_flux.dat");
  //truncate_logfile(time, "xpoint_gas_i.dat");
  //truncate_logfile(time, "xpoint_gas_e.dat");
  //truncate_logfile(time, "xpoint_cell_state.dat");
  //truncate_logfile(time, "xpoint_E3.dat");
}

void AppStateCart2::ReportAfterStep()const
{
    const dTensorBC4& aux = get_aux();
    const dTensorBC4& q = get_q();
    const double t = get_time();

    assert(dogParams.get_mcapa()<1); // without capacity function
    void WriteConservation(const dTensorBC4& q, double t);
    WriteConservation(q, t);

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
    string outputfile_e = string(get_outputdir())+"/xpoint_gas_e.dat";
    Output_center_gas_states(q, t, outputfile_i.c_str(), _rho_i-1,false);
    Output_center_gas_states(q, t, outputfile_e.c_str(), _rho_e-1,true);
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
      outputfile_e = string(get_outputdir())+"/xpoint_ent_prod_e.dat";
      output_component_at_center(aux, t, outputfile_i.c_str(),1);
      output_component_at_center(aux, t, outputfile_e.c_str(),2);
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
