#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"
#include "SLState.h"

void AfterSLFrame(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2&
    solver, SL_state sl_state )
{

    // local variables
    const int mx      = q.getsize(1);
    const int my      = q.getsize(2);
    const int meqn    = q.getsize(3);
    const int kmax    = q.getsize(4);
    const int mbc     = q.getmbc();
    const int maux    = aux.getsize(3);

    // Sync up the electric field with the current solution.  
    // TODO - I don't remember why a Poisson solve is called at the end of a
    // time step ... ?  Is is for outputting solution to file?
    void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, 
        DogSolverCart2& solver );
    BeforeStep(dt, aux, q, solver );

    return;
}
