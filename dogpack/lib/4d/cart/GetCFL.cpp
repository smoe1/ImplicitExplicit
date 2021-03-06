#include "tensors.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart4.h"
#include "DogSolverCart4.h"

double DogSolverCart4::GetCFL(double dt) const
{

    const dTensorBC5& smax = get_smax();
    const int mx = smax.getsize(1);
    const int my = smax.getsize(2);
    const int mz = smax.getsize(3);
    const int mw = smax.getsize(4);

    const double prim_vol = dogParamsCart4.get_prim_vol();
    const double dt_over_prim_vol = dt/prim_vol;
    dTensor1 cflx(mx);          // parallelization requires memory allocation
    dTensor2 cflxy(mx,my);      // parallelization requires memory allocation
    dTensor3 cflxyz(mx,my,mz);  // parallelization requires memory allocation
    cflx.setall(0.);
    cflxy.setall(0.);
    cflxyz.setall(0.);

/*
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    for (int k=1; k<=mz; k++)
    {
        double cfl_ijk = 0.;
        for (int ell=1; ell<=mz; ell++)
        {
            const double tmp = Max(Max(smax.get(i, j, k, 1), smax.get(i, j, k, 2)), smax.get(i, j, k, 3));
            const double local_cfl = tmp*dt_over_prim_vol;
            cfl_ij = Max(local_cfl, cfl_ij);
        }
        cflxy.set(i,j, cfl_ij);
    }

#pragma omp parallel for
    for (int i=1; i<=mx; i++)
    {
        double cfl_i=0.;
        for (int j=1; j<=my; j++)
        {
            cfl_i = Max(cflxy.get(i,j), cfl_i);
        }
        cflx.set(i, cfl_i);
    }

    double cfl=0.;
    for (int i=1; i<=mx; i++)
    {
        cfl = Max(cflx.get(i), cfl);
    }

    return cfl;

*/

return 0.;
}
