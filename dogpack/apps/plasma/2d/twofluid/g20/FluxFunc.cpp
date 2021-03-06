#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
// 10-moment two-fluid plasma equations
//
class dTensor2;
class dTensor3;
void FluxFunc(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor3& flux)
{
    void TwentyMomentFluxFunc(
        int n_offset,
        const dTensor2& Q,
        dTensor3& flux);
    TwentyMomentFluxFunc(0, Q, flux);
    TwentyMomentFluxFunc(20, Q, flux);

    void MaxwellFluxFunc(
        int n_offset,
        const dTensor2& Q,
        dTensor3& flux);
    MaxwellFluxFunc(40, Q, flux);

    void AdvectionFluxFunc(
        const dTensor2& Q,
        dTensor3& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc(Q, flux, _entropy_i, _rho_i);
    if(Q.getsize(2)<_entropy_e) return;
    AdvectionFluxFunc(Q, flux, _entropy_e, _rho_e);
}

void FluxFunc1(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor2& flux)
{
    void TwentyMomentFluxFunc1(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    TwentyMomentFluxFunc1(0, Q, flux);
    TwentyMomentFluxFunc1(20, Q, flux);

    void MaxwellFluxFunc1(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    MaxwellFluxFunc1(40, Q, flux);

    void AdvectionFluxFunc1(
        const dTensor2& Q,
        dTensor2& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc1(Q, flux, _entropy_i, _rho_i);
    if(Q.getsize(2)<_entropy_e) return;
    AdvectionFluxFunc1(Q, flux, _entropy_e, _rho_e);
}

void FluxFunc2(
    const dTensor2& xpts,
    const dTensor2& Q,
    const dTensor2& Aux,
    dTensor2& flux)
{
    void TwentyMomentFluxFunc2(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    TwentyMomentFluxFunc2(0, Q, flux);
    TwentyMomentFluxFunc2(20, Q, flux);

    void MaxwellFluxFunc2(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    MaxwellFluxFunc2(40, Q, flux);

    void AdvectionFluxFunc2(
        const dTensor2& Q,
        dTensor2& flux,
        int advIdx,
        int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc2(Q, flux, _entropy_i, _rho_i);
    if(Q.getsize(2)<_entropy_e) return;
    AdvectionFluxFunc2(Q, flux, _entropy_e, _rho_e);
}
