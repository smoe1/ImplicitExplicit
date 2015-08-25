#include "Components.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
// 20-moment two-fluid plasma equations
//
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

    void SymPair_MaxwellFluxFunc(
        int n_offset,
        const dTensor2& Q,
        dTensor3& flux);
    SymPair_MaxwellFluxFunc(20, Q, flux);
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

    void SymPair_MaxwellFluxFunc1(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    SymPair_MaxwellFluxFunc1(20, Q, flux);
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

    void SymPair_MaxwellFluxFunc2(
        int n_offset,
        const dTensor2& Q,
        dTensor2& flux);
    SymPair_MaxwellFluxFunc2(20, Q, flux);
}
