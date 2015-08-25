#ifndef _POISSON_SOLVE_RADIAL_H__
#define _POISSON_SOLVE_RADIAL_H__


void setGaussLegendrePoints1d(dTensor1& x1d, dTensor1& w1d);

void L2Project(int mopt, int istart, int iend,
           const dTensor2& node,
           const dTensorBC3& qin, 
           const dTensorBC3& auxin,  
           dTensorBC3& Fout,
           void (*Func)(const dTensor1&, 
            const dTensor2&, const dTensor2&, dTensor2&));

// Source term function (this is used to compute r * f(r) )
void SourceTermFunc_Extra(const dTensor1& xpts, 
        const dTensor2& qvals, 
        const dTensor2& auxvals,
        dTensor2& fvals);

// Actual source term function in -Laplacian( phi ) = f(r)
void SourceTermFunc(const dTensor1& xpts, 
        const dTensor2& qvals, 
        const dTensor2& auxvals,
        dTensor2& fvals);

///////////////////////////////////////////////////////////////////////////////
//
// Routine to compute quadrature rules on the following integrals:
//
//    M_{k,l} := 1/dr \int_{-1}^{1} phi^{(k)} phi^{(l)} / ({\xi}+{\xi}_i) d\xi
//
// It's possible to perform these integrals exactly.  For now, I'm going to
// use high-order quadrature rules.  (-DS)
//
// Parameters:
//
//    xi_i = 2 r_i / dr  ( the local transfromation )
//
// TODO - this should be an inline function to reduce function call overhead
//
///////////////////////////////////////////////////////////////////////////////
void ComputeLocalM( 
    double xi_i, const dTensor1& w1d, const dTensor1& spts, 
    const dTensor2& phi, dTensor2& M);

// See $(DOGPACK)/lib/GaussElimMatrixInv.
//
// This routine is called on each local matrix for inversion.
void GaussElimMatrixInv(dTensor2 inMat, dTensor2& outMat);

// see dog_math.cpp:
void TensorMultiply( const dTensor2& A, const dTensor2& B, dTensor2& C );

#endif
