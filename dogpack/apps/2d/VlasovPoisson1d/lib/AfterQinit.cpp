#include "dogdefs.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"
#include "DogParams.h"
#include "VlasovParams.h"
#include "DogStateCart2.h"

void AfterQinit(DogSolverCart2& solver)
{
    DogStateCart2* dogStateCart2 = &solver.fetch_state();
    dTensorBC4& qnew = solver.fetch_state().fetch_q();
    dTensorBC4& aux  = solver.fetch_state().fetch_aux();

    const int mx     = qnew.getsize(1);
    const int my     = qnew.getsize(2);
    const int maux   = aux.getsize(3);
    const int meqn   = qnew.getsize(3);
    const int kmax   = qnew.getsize(4);
    const int mbc    = qnew.getmbc();
    const int kmax1d = dogParams.get_space_order();

    // limit at all of the quadrature points
    void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
    if( dogParams.using_moment_limiter() )
    { ApplyPosLimiter(aux, qnew); }

    // make sure all Q1's are positive:
    for (int i=1-mbc; i<=mx+mbc; i++)
        for (int j=1-mbc; j<=my+mbc; j++)
            for (int me=1; me<=meqn; me++)
            {
                if( qnew.get(i,j,me,1) < -1e-12 )
                {
                    printf("  found a bad cell in AfterQinit! \n");
                    printf("  cell is located at (i,j) = (%d,%d)\n", i,j);
                    for(int k=1; k<=kmax; k++)
                    {
                        qnew.set(i,j,me,k, 0.0);
                    }
                }
            }

    // Compute the integration constant (ion density); this will be used for all future calls
    // to ComputeElectricField

    double intConst = 0.0;
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        { intConst += qnew.get(i,j,1,1)*dogParamsCart2.get_prim_vol(); }
    intConst = intConst/(dogParamsCart2.get_xhigh()-dogParamsCart2.get_xlow());
    vlasovParams.set_integration_constant(intConst);

    dTensorBC3 rho_u(mx,meqn,kmax1d,mbc,1);
    void IntegrateQ1dMoment1(const dTensorBC4& q2d, dTensorBC3& q1d);
    IntegrateQ1dMoment1(qnew, rho_u);  // first moment

    // Compute average momentum:
    intConst = 0.0;
    for(int i=1; i<=mx; i++)
    { intConst += rho_u.get(i,1,1); }
    intConst *= dogParamsCart2.get_dx()/(dogParamsCart2.get_xhigh()-dogParamsCart2.get_xlow());
    vlasovParams.set_avg_momemtum(intConst);

    printf("*** In AfterQinit, we are saving rho_u0 = %2.15e\n", vlasovParams.get_avg_momentum() );
}
