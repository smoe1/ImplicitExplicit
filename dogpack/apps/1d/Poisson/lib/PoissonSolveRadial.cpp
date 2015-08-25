#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "assert.h"
#include "DogParamsCart1.h"
#include "PoissonSolveRadial.h"

///////////////////////////////////////////////////////////////////////////////
//     
// Function to solve the radially symmetric 2D Poisson equation in 1D:
//
//                -1/r( r u' )' = f(x) on [a,b]
//                -u'(a) = gamma    u(b) = beta
//
// This equation is solved via a discontinuous Galerkin (DG) method
// of order 1,2,3,4, or 5. In the DG method, we actually solve the
// Poisson equation by writing it as a system of two first order
// differential equations:
//
//                (r E )' = r f(r) on [a,b]
//                -u'     =   E    on [a,b]
//
//                E(a) = gamma    u(b) = beta
//
// Some combinations of other boundary conditions u(a) = alpha; 
// u'(b) = delta Can be specified by first calling /PoissonSolver/ConvertBC.cpp
//
// In the case of a = 0, because the solution should be radially
// symmetric, we enforce E(a) = 0.  (TODO - not implemented)
//
// Direct integration of ( r E )' = r rho produces:
//
//      E(r) = a/r E(a) + (1/r) \int_a^r s {\rho}(s) ds
//
// Multiplying this by a test function produces the necessary
// coefficients for the electric field.  This is performed locally at each
// cell.  That is, we set E(a) = E( r_{i-1/2}^{-} ) to evaluate the electric 
// field.
//
// Parameters:
//
//     Fvals(mx, 1, morder, mbc) - Legendre weights for right-hand     
//                                 side function f(x).
//                    
//     q(mx, meqn, morder, mbc) - legendre weights for solution u.
//     q(:,mstart,:)   = 'electric field'  .... the gradient, -u_x = E
//     q(:,mstart+1,:) = 'potential'       .... the solution u(x)
//
///////////////////////////////////////////////////////////////////////////////

void PoissonSolveRadial(const int mstart, 
    const dTensor2& node,
    const double gamma, const double beta, 
    const dTensorBC3& aux, const dTensorBC3& qin,
    const dTensorBC3& F, dTensorBC3& q)
{

    //parameters
    const int mx = q.getsize(1);  assert_eq( mx, node.getsize(1)-1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);

    const int mbc = q.getmbc();

    const int morder = kmax;
    dTensorBC3 Fs( mx, meqn, kmax, mbc );

    const double dx = dogParamsCart1.get_dx();
    const double dr = dx;

    // quick error check
    if ( (mstart+1)>meqn )
    {
        cout << " ERROR in PoissonSolve_Radial: mstart is incorrectly set" << endl;
        cout << "    mstart = " << mstart << endl;
        cout << "      meqn = " << meqn << endl;
        exit(1);
    }

    // Set quadrature weights and points used for the entire integration
    // process:

    // number of points used (we can increase this if desired)
    const int mpoints = morder;

    // quadrature weights and points (in canonical variables):
    dTensor1  wgt(mpoints);             // quadrature weights
    dTensor1 spts(mpoints);             // points in the canonical variables
    setGaussPoints1d( spts, wgt );

    // Placeholder for Legendre polynomials:
    const int KMAX = 6;
    dTensor2 phi   ( mpoints, KMAX );

    // ------------------------------------------------------------------
    //
    //     FORWARD SUBSTITUTION TO SOLVE FOR  ``ELECTRIC FIELD''
    //
    // ------------------------------------------------------------------
    q.setall(0.);

    // special case for first element:
    double a   = node.get(1,1);  assert_gt( a, 0. );

    double rimh = a;
    double ri   = a+0.5*dx;
    double riph = a+1.0*dx;

    // Set up matrices for the current element:
    dTensor2 Bm( morder, morder );
    dTensor2  A( morder, morder );
    for( int i1=1; i1 <= morder; i1++ )
    for( int i2=1; i2 <= morder; i2++ )
    {
        Bm.set(i1,i2, rimh * pow(-1.0, i2+1) * sqrt( 2.0*double(i1)-1.0 )*sqrt( 2.0*double(i2)-1.0 ));
    }
    A.set(1,1, ri + 0.5*dr );
    A.set(1,2, 0.5*sq3*( 2.0*ri + dr ) );
    A.set(2,1, 0.5*sq3*( 2.0*ri*dr + dr*dr - 4.0*ri ) / dr );
    A.set(2,2, 3.0*ri + 1.5*dr - 1.0 );
    dTensor2 Ainv( morder, morder );
    GaussElimMatrixInv(A, Ainv );

    // Matrix multiplying right hand side:
    dTensor2 R(morder, morder );
    R.set(1,1, ri );
    R.set(1,2, sq3/6.0*dr );
    R.set(2,2, ri         );
    R.set(2,1, C.get(1,2) );

    TensorMultiply( Ainv, Bm, Bm );
    TensorMultiply( Ainv,  R,  R );

//  printf("checking the matrix\n");
//  for( int i=1; i <= morder ; i++ )
//  for( int j=1; j <= morder; j++ )
//  {
//      if( i == j )
//      { assert_almost_eq( R.get(i,j), 1.0 ); }
//      else
//      { assert_almost_eq( R.get(i,j), 0.0 ); }
//  }
//  printf("passed inspection\n");

    // First order method:
    q.set(1,1,1, ( gamma*(rimh/riph) + dx*(ri/riph)*F.get(1,1,1) ) );

    // TODO - this isn't correct:
    q.set(1,1,2, 0.0 );

    for( int i=2; i <= mx; i++ )
    {

        // Evaluate E at the left hand side of r_{i-1/2}:

        double rimh  = node.get(i,1   );
        double ri    = rimh + 0.5*dx;
        double riph  = ri   + 0.5*dx;

        double Ea = q.get(i-1,1,1);
        q.set(i,1,1, ( Ea*(rimh/riph) + dx*(ri/riph)*F.get(i,1,1) ) );
    }

    // ------------------------------------------------------------------
    //
    //     BACKWARD SUBSTITUTION TO SOLVE FOR  ``POTENTIAL''
    //
    // ------------------------------------------------------------------

    // special case for first element:
    q.set(mx,2,1, beta + dx*q.get(mx,1,1) );

    // every other case:
    for( int i=mx-1; i >= 1; i-- )
    for( int k=1; k <= kmax; k++ )
    {
        q.set(i,2,k, q.get(i+1,2,k) + dx*q.get(i,1,k) );
    }

}// end of function PoissonSolve.cpp


// Wrapper function for performing the following projection:
//
//     Fs := s * f(s) 
// 
// onto the Legendre coefficients
// This is a user-supplied routine that sets the source term
void SourceTermFunc_Extra(const dTensor1& xpts, 
        const dTensor2& qvals, 
        const dTensor2& auxvals,
        dTensor2& fvals)
{

    // The `real' source term function: f(r):
    SourceTermFunc( xpts, qvals, auxvals, fvals );

    // Multiply the actual source term by the extra factor of s:
    const int numpts =  xpts.getsize();
    const int   meqn = qvals.getsize(2);
    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn;   m++)
    {
        double x = xpts.get(i);
        fvals.set(i, m, x*fvals.get(i,m) );
    }

}

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
///////////////////////////////////////////////////////////////////////////////
void ComputeLocalM( 
    double xi_i, const dTensor1& w1d, const dTensor1& spts, 
    const dTensor2& phi, dTensor2& M)
{

    const int kmax    = M.getsize(1);    assert_eq( kmax,       M.getsize(2) );
    const int mpoints = phi.getsize(1);  assert_eq( mpoints, spts.getsize( ) );

    for( int k=1; k <= kmax; k++ )
    for( int l=1; l <= kmax; l++ )
    {
        double tmp = 0.;
        for( int mp=1; mp <= mpoints; mp++ )
            tmp += phi.get(mp, k) * phi.get(mp, l) / ( spts.get(mp) + xi_i );
        M.set(k, l, tmp );
    }
}
