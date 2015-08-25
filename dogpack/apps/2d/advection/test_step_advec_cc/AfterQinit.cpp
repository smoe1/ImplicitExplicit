#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"
class DogSolverCart2;
void AfterQinit(DogSolverCart2& solver)
{

    DogStateCart2* dogStateCart2 = &solver.fetch_state();
    dTensorBC4& qnew             = dogStateCart2->fetch_q();
    dTensorBC4& aux              = dogStateCart2->fetch_aux();

    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
    if( dogParams.using_moment_limiter() )
    { 
        ApplyPosLimiter(aux, qnew); 
    }

}
