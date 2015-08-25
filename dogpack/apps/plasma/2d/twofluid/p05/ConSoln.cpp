#include <assert.h>
#include "tensors.h"
#include <string>
#include "Components.h"
#include<assert.h>
#include "DogParams.h"
#include "constants.h"
#include "GEM.h"
using namespace std;

void ConSoln(
    const dTensorBC4& aux,
    const dTensorBC4& q,
    double t)
{
    increment_steps_that_have_been_accepted();
    // restrict reporting time-steps to
    // DogParams::num_subintervals many per frame
    //
    // initialized upon first call
    //
    static double next_time = t;
    static double next_frame_time = t;
    //
    if(t-next_frame_time >= -EPSILON)
    {
      next_time = next_frame_time + dogParams.get_frame_subinterval();
      next_frame_time += dogParams.get_frame_interval();
    }
    else if(t - next_time >= -EPSILON)
    {
      next_time += dogParams.get_frame_subinterval();
    }
    else
    {
      return;
    }

    void OutputMagneticFluxes(
        int n_B1,
        int n_B2,
        const dTensorBC4& q,
        double t);
    OutputMagneticFluxes(
        _B1, _B2, q, t);

    assert(dogParams.get_mcapa()<1); // without capacity function
    void WriteConservation(const dTensorBC4& q, double t);
    WriteConservation(q, t);

    void Output_center_E3(const dTensorBC4& q, double t, int _E3);
    Output_center_E3(q, t, _E3);
    void Output_center_gas_states(
      const dTensorBC4& q, double t, const char* outputfile,
      int species_offset, bool fiveMoment=false);
    const char* get_outputdir();
    string outputfile = string(get_outputdir())+"/xpoint_gas_i.dat";
    Output_center_gas_states(q, t, outputfile.c_str(),0,true);
    void Output_center_cell_state(const dTensorBC4& q, double t);
    Output_center_cell_state(q, t);
    //
    // output the rate of production of entropy at the center
    //
    const int maux = dogParams.get_maux();
    if(maux>0)
    {
      assert_eq(maux,1);
      void output_component_at_center(
        const dTensorBC4& q, double t, const char* outputfile,
        int idx);
      outputfile = string(get_outputdir())+"/xpoint_ent_prod_i.dat";
      output_component_at_center(aux, t, outputfile.c_str(),1);
    }
}

