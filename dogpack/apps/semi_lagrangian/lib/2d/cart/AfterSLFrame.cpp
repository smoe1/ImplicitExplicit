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
    const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d  = int(sqrt(mpoints));

    // Sync up the electric field with the current solution.  
//  void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, void* data );
//  BeforeStep(dt, aux, q, NULL );

    return;
}
