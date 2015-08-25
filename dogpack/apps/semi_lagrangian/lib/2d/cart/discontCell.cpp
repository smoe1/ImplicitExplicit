#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "discontCell.h"
#include "stdlib.h"
#include "DogParams.h"

// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (discontCell.cpp)
//
//  Fill in a description here...
//    This class provides a method for projecting this cell onto 
//     the legendre basis.
//
// --------------------------------------------------------------------------
//        THIS PART CURRENTLY IS UNOPERATIONAL
//        if discont_direction == 0: the integration is broken up into 4 steps.
//            cell1 = lower left cell
//            cell2 = lower right cell
//            cell3 = upper left cell
//            cell4 = upper right cell
// --------------------------------------------------------------------------
//
//
//        discont_direction == 1: the integration is broken up into 2 pieces.
//
//            cell1 = left cell
//            cell2 = right cell
//
//        discont_direction == 2: the integration is broken up into 2 pieces
//
//            cell1 = bottom cell
//            cell2 = upper cell 
//
//     Gaussian quadrature is used to compute the integration.
//
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//        Private members of class:
//
//    int mpoints, mpoints1d, kmax
//    dTesnor1 *width1, *width2;
//   dTensor1 *scut;            //location of the cuts for each cell
//    
//  int discont_direction;        //type of cell constructed.  
//
//    dTensor1 *wgt;        //gauss quad weights
//    dTensor2 *q1, *q2;    //Legendre weights on left and right half of cell
//    
//   q's values on left and right half of cell evaluated at each quadrature
//   point
//    dTensor1 *qvals1, *qvals2;
//    
//    legendre polynomial values evaluated on gauss quad grid points
//    size is (mpoints, kmax) for each of these
//    dTensor2 *phi1, *phi2, *centeredPhi2, *centeredPhi1;
//
//   // gaussian quadrature points necessary for integration
//    dTensor2 *spts(2*mpoints1d^2, 2)
//            
//   Points are labeled according to the following format:
//        1:mpoints1d - are the points in the 'first row'
//        mpoints1d+1 : 2*mpoints1d - are the hpoints in the 'second row'.
//        2mpoints1d+1 : 3*mpoints1d - are the hpoints in the 'third row'.
//
// --------------------------------------------------------------------------

discontCell::discontCell(int kmax, int discont_direction)
    // Constructor
    // POST: Create a discontinuous cell 
{

    //type of discontinuity that will occur
    this->discont_direction = discont_direction;    
    this->kmax = kmax;    //number of polynomials used

    mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    mpoints1d = int(sqrt(mpoints));
    scut      = new dTensor1(mpoints1d);    //location of cut (located in [-1,1])

    // complete list of all quadrature points //
    spts   = new dTensor2(2*mpoints1d*mpoints1d,2);

    // width of left/right cells inside [-1,1]
    width1 = new dTensor1(mpoints1d);
    width2 = new dTensor1(mpoints1d);
    wgt    = new dTensor1( 2*mpoints1d*mpoints1d );  // this is a function of widths

    ///////////////////////////////////////////////////////////////////////////
    // set up gaussian quadrature points on left and right half of cell
    //  find the 1d gaussian quad wgts and points
    w1d = new dTensor1(mpoints1d);
    x1d = new dTensor1(mpoints1d);
    set1DGaussPoints(); 
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // the integrated legendre polynomials (integrated in "rows")
    intPhi = new dTensor4(mpoints1d, kmax, kmax, 2);
    ///////////////////////////////////////////////////////////////////////////

    // this is the number of quadrature points used for the limiting process for
    // each 1D integral.
    if( dogParams.using_moment_limiter() )
    {

        // TODO - if there is a source term, how do we deal with kpts ? //
        this->Kpts         = (int) ceil( double (mpoints1d / 2.0) );
        if( dogParams.get_source_term() > 0 )
        { this->Kpts = mpoints1d; }

        this->limited_spts = new dTensor2( 2 * this->Kpts * mpoints1d, 2 );
    }

}

discontCell::~discontCell()
    // Destructor
    // POST: discontCell no longer exists
{
    delete scut;  
    delete width1; 
    delete width2;
    delete w1d; 
    delete x1d; 
    delete wgt;
    delete intPhi;
    delete spts;

    if( dogParams.using_moment_limiter() )
    { delete limited_spts; }

}

void discontCell::setScut(const dTensor1& initscut )
    // POST: Function to set the cut for the location of the discontinuity
{

    for(int n=1; n<= mpoints1d; n++)
    {

        scut->set(n, initscut.get(n) );
        // Set up areas for left and right half cells
        if ( fabs(1.0-scut->get(n) ) < 1e-17)
        {
            scut->set(n,1.0);
        } else if( fabs(-1.0-scut->get(n)) < 1e-17)
        {
            scut->set(n,-1.0);
        }
        width1->set(n, (scut->get(n)+1.0) ); 
        width2->set(n, (1.0-scut->get(n)) );
    }

    // set quadrature weights (for use with source term ... )
    // "Left" half of domain //
    int n = 0;
    for( int i=1; i <= mpoints1d; i++ )
    for( int j=1; j <= mpoints1d; j++ )
    {
        n++;
        wgt->set(n, w1d->get(i)*width1->get(j)*w1d->get(j)/2.0 );
    }

    // "Right" half of domain //
    for( int i=1; i <= mpoints1d; i++ )
    for( int j=1; j <= mpoints1d; j++ )
    {
        n++;
        wgt->set(n, w1d->get(i)*width2->get(j)*w1d->get(j)/2.0 );
    }

    // phi needs to know what scut is in order to evaluate itself
    setPhi();

}

void discontCell::set1DGaussPoints()
{

    // ---------------------------------
    // Set 1D quadrature weights and points
    // ---------------------------------
    switch ( mpoints1d )
    {
        case 1:
            w1d->set(1, 2.0 );

            x1d->set(1, 0.0 );
            break;

        case 2:
            w1d->set(1,  1.0 );
            w1d->set(2,  1.0 );

            x1d->set(1, -1.0/sq3 );
            x1d->set(2,  1.0/sq3 );
            break;

        case 3:
            w1d->set(1,  5.0/9.0 );
            w1d->set(2,  8.0/9.0 );
            w1d->set(3,  5.0/9.0 );

            x1d->set(1,  -sq3/sq5 );
            x1d->set(2,  0.0 );
            x1d->set(3, sq3/sq5 );
            break;

        case 4:
            w1d->set(1, (18.0 - sq3*sq10)/36.0 );
            w1d->set(2, (18.0 + sq3*sq10)/36.0 );
            w1d->set(3, w1d->get(2) );
            w1d->set(4, w1d->get(1) );

            x1d->set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d->set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d->set(3, -x1d->get(2) );
            x1d->set(4, -x1d->get(1) );          
            break;

        case 5:         
            w1d->set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d->set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d->set(3, 128.0/225.0 );
            w1d->set(4, w1d->get(2) );
            w1d->set(5, w1d->get(1) );

            x1d->set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d->set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d->set(3,  0.0 );
            x1d->set(4, -x1d->get(2) );
            x1d->set(5, -x1d->get(1) );
            break;
    }


}// end of setGausspts

///////////////////////////////////////////////////////////////////////////////
// 
//    Routine for computing the projection onto the legendre basis.
//
// Parameters:
// 
//   dTensor2(mpoints1d, kmax) q1new = coefficient weights on legendre polys for left
//   half of grid cell
//
//   dTensor2(mpoints1d, kmax) q2new = coefficient weights on legendre polys for right
//   half of grid cell
//
//  Returns:
//
//   dTensor1(kmax) qnew = coefficient weights on legendre polys over entire
//   cell.
//
///////////////////////////////////////////////////////////////////////////////
void discontCell::project(const dTensor2& q1new, const dTensor2& q2new, 
        dTensor1& qnew)
{

    for(int k = 1; k <= kmax; k++)
    {
        double tmp = 0.0;
        for(int j = 1; j<= mpoints1d; j++ )
            for(int n = 1; n <= kmax; n++ )
            {
                tmp += intPhi->get(j,k,n,1)*q1new.get(j,n) +
                       intPhi->get(j,k,n,2)*q2new.get(j,n);
            }

        qnew.set( k, tmp / 4.0 );

    }//end of setting weights;

    ///////////////// TODO ////////////////////////////////////////////////
    // This is here simply for debugging                                 //
////if( !dogParams.get_source_term()  & 
////    qnew.get(1) < -1e-10 & dogParams.using_moment_limiter() )
////{
////    printf("****   saving a negative Q1 in discontCell.cpp: ****\n");

////    for( int k=1; k <= kmax; k++ )
////    {
////        printf("     qnew.get( %2d ) = %+2.8e\n", k, qnew.get(k));
////    }

////    for(int j = 1; j<= mpoints1d; j++ )
////    for( int k=1; k <= kmax; k++ )
////    {
////        printf("     q1.get( %2d, %2d ) = %+2.8e\n", j, k, q1new.get(j,k));
////        printf("     q2.get( %2d, %2d ) = %+2.8e\n", j, k, q2new.get(j,k));
////    }

////    printf("*****  widths are given by: ****\n");
////    for( int j=1; j <= mpoints1d; j++ )
////    {
////        printf("  width1->get( %d ) = %2.3e\n", j, width1->get(j) );
////        printf("  width2->get( %d ) = %2.3e\n", j, width2->get(j) );
////        assert( fabs( width1->get(j) + width2->get(j) - 2.0 ) < 1e-13 );
////    }

////    dTensor2* limited_spts = this->get_limited_spts();
////    printf("  discont_direction = %d\n", discont_direction );
////    for( int i=1; i <= limited_spts->getsize(1); i++ )
////    {
////        double xi  = limited_spts->get(i,1);
////        double eta = limited_spts->get(i,2);
////        printf("i=%2d.   xi = %+2.8e, eta = %+2.8e \n", i, xi,eta);
////    }

////    // force a quit
////    exit(1);
////}
///////////////////// TODO ////////////////////////////////////////////////

}// end of project

//////////////////////////////////////////////////////////////////////////////
//  This function is essentially a duplicate of what's in ApplyPosLimiter and
//  should probably be combined.
//////////////////////////////////////////////////////////////////////////////
void discontCell::evalPhi( const dTensor2& loc_spts, dTensor2& phi)
{
    const int mpoints     = loc_spts.getsize(1);
    const int space_order = mpoints1d;

    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi  = loc_spts.get(m,1);
        const double eta = loc_spts.get(m,2);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1].
        switch( space_order )
        {
            case 5:  // fifth order                                 
                phi.set( 15, m, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( 14, m, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi.set( 13, m, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi.set( 12, m, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set( 11, m, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 4:  // fourth order
                phi.set( 10, m,  sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( 9,  m,   sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( 8,  m,   sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( 7,  m,   sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 3:  // third order
                phi.set( 6, m,   sq5*(1.5*eta2 - 0.5) );
                phi.set( 5, m,   sq5*(1.5*xi2 - 0.5) );
                phi.set( 4, m,   3.0*xi*eta );                  

            case 2:  // second order                
                phi.set( 3, m,  sq3*eta );
                phi.set( 2, m,  sq3*xi  );

            case 1:  // first order
                phi.set( 1, m, 1.0 );

                break;                

            default:
                unsupported_value_error(space_order);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Evalutate the function, phi(k)(xi,eta)
// not pretty because of hte intensive function calls, but this works 
///////////////////////////////////////////////////////////////////////////////
inline double discontCell::evalPhi(const int& k, const double &xi, 
        const double &eta)
{

    const double xi2 = xi*xi;
    const double xi3 = xi*xi2;
    const double xi4 = xi*xi3;
    const double eta2 = eta*eta;
    const double eta3 = eta*eta2;
    const double eta4 = eta*eta3;

    // Legendre basis functions evaluated at (xi,eta) in the
    // interval [-1,1]x[-1,1].
    switch( k )
    {
        case 15: 
            return ( 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );

        case 14: 
            return(  105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );

        case 13: 
            return(  5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );

        case 12: 
            return(  sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );

        case 11: 
            return(  sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

        case 10:
            return(  sq7*(2.5*eta3 - 1.5*eta) );

        case 9:
            return(  sq7*(2.5*xi3 - 1.5*xi) );

        case 8:
            return(   sq3*sq5*xi*(1.5*eta2 - 0.5) );

        case 7: 
            return(   sq3*sq5*eta*(1.5*xi2 - 0.5) );

        case 6:
            return(   sq5*(1.5*eta2 - 0.5) );

        case 5:
            return(   sq5*(1.5*xi2 - 0.5) );

        case 4:
            return(   3.0*xi*eta );            

        case 3:
            return(  sq3*eta );

        case 2:
            return(  sq3*xi  );

        case 1: 
            return(  1.0 );

    }

    return 0.0;

}// end of function evalPhi
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Function for evaluating the integral of q(xi,eta) in one direction.
//
//  For discontinuity in the 'x'direction, we have:
//
//  F^(k) ( eta_j ) = \sum_n ( Q1^(n) * intPhi(j,k,1,n) + 
//                             Q2^(n) * intPhi(j,k,2,n).)
//
//  intPhi(j,k,1,n) = \int_{-1}^s \phil^(n) \phi^k
//  intPhi(j,k,2,n) = \int_{s}^1  \phir^(n) \phi^k
//
//  Where j is the index of the 'row'
///////////////////////////////////////////////////////////////////////////////
void discontCell::setPhi()
{

    int k, n, i, j;
    double tmp1, tmp2;

    //x1( i, j, 1 ) = values for phi^(k)  (centered poly)
    //x1( i, j, 2 ) = values for phi^(n)  (off centered poly)
    dTensor3 x1( mpoints1d, mpoints1d, 2 );  // quadrature points for 'left'
    dTensor3 x2( mpoints1d, mpoints1d, 2 );  // quadrature points on 'right' 

    // NOTE: The "points" where we need to evaluate phi are not the same
    // as the points where we evaluate the test function.
    // For example: if scut = 0.99, then we actually evaluate the old phi from
    // [-1.0, -0.99] rather than [0.99,1.0]
    for( i=1; i<= mpoints1d; i++ )
        for( j=1; j<= mpoints1d; j++ )
        {
            tmp1 =  (scut->get(j) - (-1.0) ) / 2.0 * x1d->get(i) +
                (scut->get(j) + (-1.0) ) / 2.0;
            tmp2 =  ( 1.0 - scut->get(j) ) / 2.0 * x1d->get(i) +
                (scut->get(j) + 1.0 ) / 2.0;

            x1.set(i, j, 1, tmp1);
            x2.set(i, j, 1, tmp2);

            tmp1 -= (scut->get(j) - 1.0 );
            tmp2 -= (scut->get(j) + 1.0 );

            x1.set(i, j, 2, tmp1);
            x2.set(i, j, 2, tmp2);

        }

    for( j=1; j<=mpoints1d; j++)
        for( k=1; k<=kmax; k++)
            for( n=1; n<=kmax; n++)
            {
                tmp1 = 0.0;
                tmp2 = 0.0;
                for( i=1; i<= mpoints1d; i++)
                {

                    switch(discont_direction)
                    {
                        case 1:  //advect in x-direction
                            tmp1 += w1d->get(i)*width1->get(j)/2.0 *
                                evalPhi(k, x1.get(i, j, 1), x1d->get(j) ) *
                                evalPhi(n, x1.get(i, j, 2), x1d->get(j) );
                            tmp2 += w1d->get(i)*width2->get(j)/2.0 *
                                evalPhi(k, x2.get(i, j, 1), x1d->get(j) ) *
                                evalPhi(n, x2.get(i, j, 2), x1d->get(j) );
                            break;

                        case 2: //advect in y-direction
                            tmp1 += w1d->get(i)*width1->get(j)/2.0 *
                                evalPhi(k, x1d->get(j), x1.get(i, j, 1) ) *
                                evalPhi(n, x1d->get(j), x1.get(i, j, 2) );
                            tmp2 += w1d->get(i)*width2->get(j)/2.0 *
                                evalPhi(k, x1d->get(j), x2.get(i, j, 1) ) *
                                evalPhi(n, x1d->get(j), x2.get(i, j, 2) );
                            break;
                    }// end of switch statement

                }

                intPhi->set( j, k, n, 1, w1d->get(j) * tmp1 );
                intPhi->set( j, k, n, 2, w1d->get(j) * tmp2 );

            }

    double xi, eta;

    // "Left" half of domain //
    k = 0;
    for( i=1; i<= mpoints1d; i++ )
        for( j=1; j<= mpoints1d; j++ )
        {

            k++;

            xi = (scut->get(j) - (-1.0) ) / 2.0 * x1d->get(i) +
                 (scut->get(j) + (-1.0) ) / 2.0;

            eta = x1d->get(j);
            spts->set(k, 1, xi);
            spts->set(k, 2, eta);

        }

    // "Right" half of domain //
    for( i=1; i<= mpoints1d; i++ )
        for( j=1; j<= mpoints1d; j++ )
        {

            k++;

            xi = (1.0 - scut->get(j) ) / 2.0 * x1d->get(i) +
                 (1.0 + scut->get(j) ) / 2.0;
            eta = x1d->get(j);

            spts->set(k, 1, xi);
            spts->set(k, 2, eta);

        }

}

///////////////////////////////////////////////////////////////////////////////
// This routine simply translates all the quadrature points 
// from (xi, eta) -> (x,y).
//
// One needs to supply an xc, yc and dx, dy in order to do this.
//
// This function is necessary for adding a source term.
//
///////////////////////////////////////////////////////////////////////////////
void discontCell::set_xpts(dTensor1* xc_vec, dTensor1* dx_vec, dTensor2* xpts)
{

    // cell centers and widths //
    const double xc = xc_vec->get(1);
    const double yc = xc_vec->get(2);
    const double dx = dx_vec->get(1);
    const double dy = dx_vec->get(2);

    const int mpts  = xpts->getsize(1);

    // "Left" and "Right" half of domain //
    for( int m=1; m <= mpts; m++ )
    {
            xpts->set(m, 1, 0.5*dx*spts->get(m,1) + xc );
            xpts->set(m, 2, 0.5*dy*spts->get(m,2) + yc );
    }

}

///////////////////////////////////////////////////////////////////////////////
// This routine simply translates all the quadrature points 
// from (xi, eta) -> (x,y).
//
// The quadrature points here are the ones used in the limiting process.
//
// One needs to supply an xc, yc and dx, dy in order to do this.
//
// This function is necessary for adding a source term.  get_limited_spts needs
// to be called before this function can be called
//
///////////////////////////////////////////////////////////////////////////////
void discontCell::set_limited_xpts(dTensor1* xc_vec, dTensor1* dx_vec, dTensor2* xpts)
{

    // cell centers and widths //
    const double xc = xc_vec->get(1);
    const double yc = xc_vec->get(2);
    const double dx = dx_vec->get(1);
    const double dy = dx_vec->get(2);

    const int mpts  = xpts->getsize(1);

    // "Left" and "Right" half of domain //
    for( int m=1; m <= limited_spts->getsize(1); m++ )
    {
        xpts->set(m, 1, 0.5*dx*limited_spts->get(m,1) + xc );
        xpts->set(m, 2, 0.5*dy*limited_spts->get(m,2) + yc );
    }

}

// integrate out the time //
void discontCell::integrate_psi_in_time(const dTensor3& psi, dTensor2& psi_new)
{

    const int mpts    = psi.getsize(1);
    const int meqn    = psi.getsize(2);

    for( int n=1; n <= mpts; n++ )
    for( int me=1; me <= meqn; me++ )
    {
        double tmp = 0.0;
        for( int t=1; t <= mpoints1d; t++ )
        {
            tmp += psi.get(n, me, t) * w1d->get(t);
        }
        psi_new.set(n, me, tmp );
    }

}

///////////////////////////////////////////////////////////////////////////////
// this routine accomodates time integration.  The format of psi is:
//
//   psi( 1:mpts, 1:meqn, 1:mpoints1d )
//
// The third component is psi evaluated at the time quadrature points
///////////////////////////////////////////////////////////////////////////////
void discontCell::project_source_w_time(const dTensor3& psi, dTensor2& psi_new)
{

    const int mpts = psi.getsize(1);
    const int meqn = psi.getsize(2);

    // perform the projection now that we have psi at each quadrature point //
    dTensor2 phi( kmax,  mpts);
    evalPhi     ( *spts, phi );

    for( int me=1; me <= meqn; me++ )
    for( int k=1; k <= kmax; k++ )
    {
        double tmp = 0.0;
        for( int n=1; n <= mpts; n++ )
        { // spatial integration
            for( int t=1; t <= mpoints1d; t++ )
            { // time integration
                
                tmp += w1d->get(t) * wgt->get(n) * psi.get(n,me,t) * phi.get(k, n);
            }
        }
        psi_new.set(me, k, tmp / 4.0 );
    }
 
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// The difference between this routine and the other one is that psi has already
// been integrated in time
///////////////////////////////////////////////////////////////////////////////
void discontCell::project_source_w_time(const dTensor2& psi, dTensor2& psi_new)
{

    const int mpts = psi.getsize(1);
    const int meqn = psi.getsize(2);

    // perform the projection now that we have psi at each quadrature point //
    dTensor2 phi( kmax,  mpts);
    evalPhi     ( *spts, phi );

    for( int me=1; me <= meqn; me++ )
    for( int k=1; k <= kmax; k++ )
    {
        double tmp = 0.0;
        for( int n=1; n <= mpts; n++ )
        { // spatial integration
            tmp += wgt->get(n) * psi.get(n,me) * phi.get(k, n);
        }
        psi_new.set(me, k, tmp / 4.0 );
    }
 
}
///////////////////////////////////////////////////////////////////////////////

dTensor1* discontCell::get_x1d()
{ return x1d; }

dTensor1* discontCell::get_w1d()
{ return w1d; }

dTensor2* discontCell::get_spts()
{ return spts; }

dTensor2* discontCell::get_limited_spts()
{

    dTensor1* x1dloc = new dTensor1(Kpts);
    switch ( Kpts )
    {
        case 1:
            x1dloc->set(1, 0.0 );
            break;

        case 2:
            x1dloc->set(1, -1.0/sq3 );
            x1dloc->set(2,  1.0/sq3 );
            break;

        case 3:
            x1dloc->set(1,  -sq3/sq5 );
            x1dloc->set(2,  0.0 );
            x1dloc->set(3, sq3/sq5 );
            break;

        case 4:
            x1dloc->set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1dloc->set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1dloc->set(3, -x1dloc->get(2) );
            x1dloc->set(4, -x1dloc->get(1) );          
            break;

        case 5:         
            x1dloc->set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1dloc->set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1dloc->set(3,  0.0 );
            x1dloc->set(4, -x1dloc->get(2) );
            x1dloc->set(5, -x1dloc->get(1) );
            break;

        default:
            printf("  bad Kpts chosen\n");
            exit(1);

    }

    dTensor1* scut_old = new dTensor1( scut->getsize() );
    for( int m=1; m <= mpoints1d; m++ )
    {
        // note: the points where the legendre polys are evaluated at time t^n
        // are different than the points where these are evaluated at time
        // t^{n+1}.
        //
        // I.e. if scut = 0.1, then [0.1, 1.0] gets mapped to [-1.0, -0.1]
        scut_old->set(m, -1.0 * scut->get(m) );
    }

    int m = 0;

    // "left" half of cell //
    for( int i = 1; i <= mpoints1d; i++ )
    for( int j = 1; j <= Kpts; j++ )
    {

        m++;
        double eta = x1d->get(i);
        double xi  = 0.5 * ( x1dloc->get(j) * (1.0 + scut_old->get(i)) + scut_old->get(i)-1.0);

        switch(discont_direction)
        {
            case 1:  //advect in x-direction
                
                limited_spts->set(m, 1, xi);
                limited_spts->set(m, 2, eta);
                break;

            case 2: //advect in y-direction

                limited_spts->set(m, 1, eta);
                limited_spts->set(m, 2, xi);
                break;

        }// end of switch statement

    }

    // "right" half of cell //
    for( int i = 1; i <= mpoints1d; i++ )
    for( int j = 1; j <= Kpts; j++ )
    {
        m++;
        double eta = x1d->get(i);
        double xi  = 0.5 * ( x1dloc->get(j) * (1.0 - scut_old->get(i)) + scut_old->get(i) + 1.0);

        switch(discont_direction)
        {
            case 1:  //advect in x-direction
                
                limited_spts->set(m, 1, xi);
                limited_spts->set(m, 2, eta);
                break;

            case 2: //advect in y-direction

                limited_spts->set(m, 1, eta);
                limited_spts->set(m, 2, xi);
                break;

        }// end of switch statement

    }

    delete x1dloc;
    delete scut_old;

    return limited_spts;

}
