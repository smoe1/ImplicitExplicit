#include <stdio.h>
#include <stdlib.h>
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"
//////////////////////////////////////////////////////////////////////////////
// Function for stepping 1D advection equation forward in time.
//
// The advection equation: q_t + u q_x = 0 is solved by projecting the exact
// solution, qex(x,t) = qn( x-v*t ) against all the basis functions.  Care is
// taken to do the integration exactly.
//
//////////////////////////////////////////////////////////////////////////////
void StepAdvec1D(
    double dt, double xlow, double dx,
    const dTensor1& speeds,
    const dTensorBC3& qold, dTensorBC3& qnew)
{

    //-local parameters -----------------------------------------------------//
    const int melems  = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

    // Twice as many points are required to perform the projections
    const int numpts = 2*kmax;

    // compute quadrature points for where Q needs to be sampled
    const double s_area = 2.0;
    for(int me=1; me<= meqn; me++)
    {

        ////////////////// Set up speeds for current equation /////////////
        const double large_eta = speeds.get(me)*dt/dx;
        const double frac_eta  = (double) ( large_eta - floor( large_eta ) );

        // integration points and weights
        dTensor1 wgt(numpts), spts(numpts);
        void setIntegrationPoints1D(double frac_eta, dTensor1& spts, dTensor1& wgt);
        setIntegrationPoints1D(frac_eta, spts, wgt);

        dTensor1 spts_old ( numpts );
        iTensor1 ishift   ( numpts );
        void translateXi( 
            double eta, const dTensor1& xi, dTensor1 &xinew, iTensor1& ishift );
        translateXi( large_eta, spts, spts_old, ishift );

        // Legendre Basis functions, their derivatives, evaluated at all the
        // necessary poitns
        dTensor2 phi(numpts, 5);
        dTensor2 phi_old(numpts, 5);
        void evaluateLegendrePolys1D( const dTensor1& spts, dTensor2& phi );
        evaluateLegendrePolys1D( spts,     phi     );
        evaluateLegendrePolys1D( spts_old, phi_old );



        // This part isn't parallelized, because this routine is assumed to be
        // `fast' - as in it's called many times from a higher dimensional code
        // ...
        // #pragma omp parallel for
        for(int i=1; i<=melems; i++)
        {

            for( int k=1; k <= kmax; k++)
            {
                double sum = 0.0;
                for( int m=1; m <= numpts; m++ )
                {

                    // periodicity is enforced here!
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

        }
    }

    //void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    //SetBndValues(node, auxvals, qnew);

}

void setIntegrationPoints1D(double frac_eta, dTensor1& spts, dTensor1& wgt)
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

// function that evaluates legendre polynomials and all their derivatives
// at the list of points given by spts
void evaluateLegendrePolys1D( 
    const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x)
{

    const int numpts = spts.getsize();

    for (int m=1; m<=numpts; m++)
    {
        // Legendre basis functions at each grid point
        phi.set( m,1, 1.0 );
        phi.set( m,2, sq3*spts.get(m) );
        phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
        phi.set( m,4, 0.5*sq7*spts.get(m)
                *(5.0*pow(spts.get(m),2) - 3.0) );
        phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );

        // 1st derivative of Legendre basis functions at each grid point
        phi_x.set( m,1, 0.0 );
        phi_x.set( m,2, 2.0*sq3/dx );
        phi_x.set( m,3, 6.0*sq5*spts.get(m)/dx );
        phi_x.set( m,4, 3.0*sq7*(5.0*pow(spts.get(m),2)-1.0)/dx );
        phi_x.set( m,5, 15.0*spts.get(m)*
                (7.0*pow(spts.get(m),2)-3.0)/dx );
    }
}

// function that evaluates legendre polynomials and all their derivatives
// at the list of points given by spts
void evaluateLegendrePolys1D(const dTensor1& spts, dTensor2& phi )
{

    const int numpts = spts.getsize();

    for (int m=1; m<=numpts; m++)
    {
        // Legendre basis functions at each grid point
        phi.set( m,1, 1.0 );
        phi.set( m,2, sq3*spts.get(m) );
        phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
        phi.set( m,4, 0.5*sq7*spts.get(m)
                *(5.0*pow(spts.get(m),2) - 3.0) );
        phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );

    }
}

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
