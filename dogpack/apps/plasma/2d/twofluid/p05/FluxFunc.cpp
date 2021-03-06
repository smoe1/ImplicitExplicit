#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
// 5-moment two-fluid plasma equations
//
class dTensor2;
class dTensor3;
void FluxFunc(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor3& flux)
{
    void FiveMomentFluxFunc(
        int n_offset,
        const dTensor2& Q,
        dTensor3& flux);
    FiveMomentFluxFunc(0, Q, flux);

    void SymPair_MaxwellFluxFunc(
        int n_offset,
        const dTensor2& Q,
        dTensor3& flux);
    SymPair_MaxwellFluxFunc(5, Q, flux);

    void AdvectionFluxFunc(
        const dTensor2& Q,
        dTensor3& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
        AdvectionFluxFunc(Q, flux, _entropy_i, _rho_i);
}

void FluxFunc1(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor2& flux)
{
    void FiveMomentFluxFunc1(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    FiveMomentFluxFunc1(0, Q, flux);

    void SymPair_MaxwellFluxFunc1(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    SymPair_MaxwellFluxFunc1(5, Q, flux);

    void AdvectionFluxFunc1(
        const dTensor2& Q,
        dTensor2& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc1(Q, flux, _entropy_i, _rho_i);
}

void FluxFunc2(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor2& flux)
{
    void FiveMomentFluxFunc2(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    FiveMomentFluxFunc2(0, Q, flux);

    void SymPair_MaxwellFluxFunc2(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    SymPair_MaxwellFluxFunc2(5, Q, flux);

    void AdvectionFluxFunc2(
        const dTensor2& Q,
        dTensor2& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc2(Q, flux, _entropy_i, _rho_i);
}
