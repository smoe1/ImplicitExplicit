#if 0
#include "tensors.h"
#include "gas10.h"
#include "Components.h"
#include "gasCart2.h"
#include "DogSolverCart2.h"

void reset_entropy(
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
}

// Function that is called after a full time step
// (i.e., after all stages are complete)
void AfterFullTimeStep(DogSolverCart2& solver)
{
  UpdateEntropies(solver.get_dogState().get_dt(),
    solver.fetch_aux(),solver.fetch_q(), _entropy_i,1,&reset_entropy);
}
#endif
