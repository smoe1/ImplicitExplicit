#include <stdio.h>
#include <stdlib.h>
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"
#include "DogParams.h" // for reading limiter

// Function that traces a quadrature point back in time to its old value, and
// computes the index of its shift.  It is expected that xi \in [-1,1]
void translateXi( double large_eta, const dTensor1& xi, dTensor1 &xinew, iTensor1& ishift )
{
    
    for( int m=1; m <= xi.getsize(); m++ )
    {
        double xi_old = xi.get(m) - 2.0 * large_eta;
        ishift.set( m, (int) floor( 0.5 * (xi_old+1.0) ) );
        xinew.set(  m, xi_old - 2.0 * ( ishift.get(m) ) );
    }

}

void StepAdvec(const double& dt, const dTensor2& node,
  dTensorBC3& auxvals, dTensorBC1& smax,
  dTensorBC3& qold, // limiter will modify qold if turned on
  dTensorBC3& qnew)
{
    //////////////////////////////////////////////////////////////////////////////
    // Function for stepping advection equation forward in time.
    //
    // The advection equation: q_t + u q_x = 0 is solved by projecting the exact
    // solution, qex(x,t) = qn( x-v*t ) against all the basis functions.  Care is
    // taken to do the integration exactly.
    //
    // It is assumed that aux(1,1,1) = u is constant.
    //
    //////////////////////////////////////////////////////////////////////////////


    //-local parameters -----------------------------------------------------//
    const int melems  = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

    //save the center of each grid cell
    const double dx = node.get(2,1)-node.get(1,1);
    const double speed = auxvals.get(1,1,1);

    for(int j=1-mbc; j<=melems+mbc; j++)
    { smax.set(j, Max(smax.get(j), fabs(speed) ) ); }

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points for the unit interval
    //////////////////////////////////////////////////////////////////////////

    const int numpts = 2*kmax;

    // compute quadrature points for where Q needs to be sampled
    const double s_area = 2.0;
    const double large_eta = speed*dt/dx;
    const double frac_eta  = (double) ( large_eta - floor( large_eta ) );

    // integration points and weights
    dTensor1 wgt(numpts), spts(numpts);
    void setIntegrationPoints(double frac_eta, dTensor1& spts, dTensor1& wgt);
    setIntegrationPoints(frac_eta, spts, wgt);

    dTensor1 spts_old ( numpts );
    iTensor1 ishift   ( numpts );
    void translateXi( 
        double eta, const dTensor1& xi, dTensor1 &xinew, iTensor1& ishift );
    translateXi( large_eta, spts, spts_old, ishift );

    // Legendre Basis functions, their derivatives, evaluated at all the
    // necessary poitns

    const int MAX_KMAX = 6;  // maximum possible supported order of accuracy ...

    dTensor2 phi(numpts, MAX_KMAX), phi_x(numpts,MAX_KMAX);
    dTensor2 phi_old(numpts, MAX_KMAX), phi_x_old(numpts,MAX_KMAX);
    void evaluateLegendrePolys( 
            const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x);
    evaluateLegendrePolys( dx, spts,     phi,     phi_x);
    evaluateLegendrePolys( dx, spts_old, phi_old, phi_x_old);


    // Check and apply the limiter if necessary --------------- //
    if(dogParams.using_moment_limiter())
    {

        for(int i=1-mbc; i<=melems+mbc; i++)
        for(int me=1; me <= meqn; me++ )
        {

            // -------------------------------------------------------------- //
            // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
            // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
            // -------------------------------------------------------------- //
            double m = qold.get(i,me,1);
            for(int mp=1; mp <= phi_old.getsize(1); mp++)
            {
                // evaluate q at spts(mp) //
                double qnow = 0.0;
                for( int k=1; k <= kmax; k++ )
                {
                    qnow += qold.get(i,me,k) * phi_old.get(mp,k);
                }
                m = Min(m, qnow );
            }

            double Q1 = qold.get(i,me,1);

            double theta = 0.0;
            if( fabs( Q1 - m ) < EPSILON ) { theta = 0.0; }
            else{ theta = Min( 1.0, (Q1 - EPSILON) / (Q1 - m) ); }

            // limit q //
            for( int k=2; k <= kmax; k++ )
            {
                qold.set(i,me,k, qold.get(i,me,k) * theta );
            }

        }

    }
    // End of limiter //
    // //////////////////////////////////////////// //

    #pragma omp parallel for
    for(int i=1; i<=melems; i++)
    {

        //loop over each equation
        for(int me=1; me<= meqn; me++)
        {

            for( int k=1; k <= kmax; k++)
            {
                double sum = 0.0;
                for( int m=1; m <= numpts; m++ )
                {

                    // hard coded periodicity is enforced here!
                    int io = (int)( i + ishift.get(m) );
                    io = iMod((io-1), melems)+1;

                    // evaluate q at the old point:
                    double qval = 0.0;
                    for( int kq=1; kq <= kmax; kq++ )
                    {
                        qval += qold.get(io, me, kq) * phi_old.get(m, kq);
                    }
                    sum += wgt.get(m) * qval * phi.get(m, k);
                }
                qnew.set(i, me, k, sum / s_area );
            }

        }//end of loop over each cell
    }//end of loop over each eqn  

    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    SetBndValues(node, auxvals, qnew);

}

void setIntegrationPoints(double frac_eta, dTensor1& spts, dTensor1& wgt)
{

    const int mpts   = wgt.getsize();
    const int morder = (int) (mpts / 2);

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points (standard before a transformation)
    //////////////////////////////////////////////////////////////////////////
    dTensor1 w1d(morder), x1d(morder);
    switch ( morder )
    {

        case 1:
            w1d.set(1,  2.0e0 );
            x1d.set(1, 0.0e0 );

            break;

        case 2:
            w1d.set(1,   1.0 );
            w1d.set(2,   1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );

            break;

        case 3:
            w1d.set(1, 5.0e0/9.0e0 );
            w1d.set(2, 8.0e0/9.0e0 );
            w1d.set(3, 5.0e0/9.0e0 );

            x1d.set(1, -sq3/sq5 );
            x1d.set(2,  0.0e0 );
            x1d.set(3,  sq3/sq5 );

            break;

        case 4:
            w1d.set(1, (18.0 - sqrt(30.0))/36.0 );
            w1d.set(2, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(3, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(4, (18.0 - sqrt(30.0))/36.0 );

            x1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            x1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

            break;

        case 5:
            w1d.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
            w1d.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(3,  128.0/225.0 );
            w1d.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

            x1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            x1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(3,  0.0 );
            x1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

            break;

        case 6:

            w1d.set(1,  0.171324492379170);
            w1d.set(2,  0.360761573048139);
            w1d.set(3,  0.467913934572691);
            w1d.set(4,  0.467913934572691);
            w1d.set(5,  0.360761573048139);
            w1d.set(6,  0.171324492379170);


            x1d.set(1,  0.932469514203152);
            x1d.set(2,  0.661209386466265);
            x1d.set(3,  0.238619186083197);
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );


    }

    // location of discontinuity and width of left and right hand side
    const double scut      = -1.0 + 2.0*frac_eta;
    const double left_len  = scut - (-1.0);
    const double right_len = 1.0  - scut;

    // left half of the integral
    const double bl = scut;
    const double al = -1.0;

    // right half of the integral
    const double br = 1.0;
    const double ar = scut;

    // modify weights and find modified quadrature points to evaluate q
    for( int m=1; m <= morder; m++ )
    {

        // weights
        wgt.set(m,           left_len / 2.0 * w1d.get(m) );
        wgt.set(m + morder, right_len / 2.0 * w1d.get(m) );

        // points
        double x = x1d.get(m);
        spts.set(m,         0.5*( (bl-al)*x + (al+bl) ) );
        spts.set(morder+m,  0.5*( (br-ar)*x + (ar+br) ) );
    }

}
