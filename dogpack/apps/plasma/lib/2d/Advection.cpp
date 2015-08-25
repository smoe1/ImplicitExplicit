#include "dogdefs.h"

// advect the scalar S using the velocity field
void AdvectionFluxFunc(
    const dTensor2& Q,
    dTensor3& flux,
    int advIdx,
    int rhoIdx)
{
    if(Q.getsize(2)<advIdx) return;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        const double S = Q.get(j,advIdx);
        const double rho_inv = 1./Q.get(j,rhoIdx);
        const double u1 = Q.get(j,rhoIdx+1)*rho_inv;
        const double u2 = Q.get(j,rhoIdx+2)*rho_inv;
        flux.set(j,advIdx,1, u1*S);
        flux.set(j,advIdx,2, u2*S);
    }
}

// advect the scalar S using the velocity field
void AdvectionFluxFunc1(
    const dTensor2& Q,
    dTensor2& flux,
    int advIdx,
    int rhoIdx)
{
    if(Q.getsize(2)<advIdx) return;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        const double S = Q.get(j,advIdx);
        const double rho = Q.get(j,rhoIdx);
        const double u1 = Q.get(j,rhoIdx+1)/rho;
        flux.set(j,advIdx, u1*S);
    }
}

// advect the scalar S using the velocity field
void AdvectionFluxFunc2(
    const dTensor2& Q,
    dTensor2& flux,
    int advIdx,
    int rhoIdx)
{
    if(Q.getsize(2)<advIdx) return;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        const double S = Q.get(j,advIdx);
        const double rho = Q.get(j,rhoIdx);
        const double u2 = Q.get(j,rhoIdx+2)/rho;
        flux.set(j,advIdx, u2*S);
    }
}
