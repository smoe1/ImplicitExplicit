#include <cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include <iostream>
#include "DogSolverCart2.h"

void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver)
{

    // local variables
    double x,y;
    const int mx      = q.getsize(1);
    const int my      = q.getsize(2);
    const int meqn    = q.getsize(3);
    const int kmax    = q.getsize(4);
    const int mbc     = q.getmbc();
    const int maux    = aux.getsize(3);
    const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d  = int(sqrt(mpoints));

    // redo the velocities.  these are set in AuxFunc //
    void L2Project(int istart, int iend, int jstart, int jend,
               const dTensorBC4* qin,
               const dTensorBC4* auxin, dTensorBC4* Fout,
               void (*Func)(const dTensor2&,const dTensor2&,
                            const dTensor2&,dTensor2&));
    void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
             const dTensor2& NOT_USED_2, dTensor2& auxvals);
    L2Project(1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, aux, &AuxFunc );

}


void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q,  DogSolverCart2& solver, void* data)
{

    double t = *( (double*)data );
    BeforeStep(dt,aux,q,solver);

}
