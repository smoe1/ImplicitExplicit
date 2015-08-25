#include "stdlib.h"
#include "stdio.h"
#include "dog_math.h"
#include "constants.h"
#include "dg_cell.h"

dg_cell::dg_cell( )
{
    // Constructor
    // POST: Create a dG cell that is able to handle integration of any
    // polynomial

}

dg_cell::~dg_cell()
{
    // Destructor
    // POST: Cell no longer exists

} // end of destructor


// function that sets numpts quadrature points into w1d and x1d //
void dg_cell::SetW1dX1d(dTensor1 &w1d, dTensor1& x1d)
{

    switch ( w1d.getsize() )
    {
        case 1:
            w1d.set(1, 2.0e0 );
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
}

// evaluate phi(1) - phi(5) at all the points provided by spts //
void dg_cell::evaluatePolys( const dTensor1& spts, dTensor2& phi )
{

    const int numpts = spts.getsize();
    const int kmax   = phi.getsize(2);

    for (int m=1; m<=numpts; m++)
    {
        // Legendre basis functions at each grid point

        switch( kmax )
        {
            case 5:
                phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                        - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );

            case 4:
                phi.set( m,4, 0.5*sq7*spts.get(m)
                        *(5.0*pow(spts.get(m),2) - 3.0) );

            case 3:
                phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );

            case 2:
                phi.set( m,2, sq3*spts.get(m) );

            case 1:
                phi.set( m,1, 1.0 );
                break;

            default:
                printf("  shame on you!\n");
        }

    }
}

void dg_cell::integratePhi(double a, double b, dTensor1& phik )
{

    for( int k=1; k <= phik.getsize(); k++ )
    { phik.set(k, integratePhi(a,b,k) ); }

}


double dg_cell::integratePhi(double a, double b, int k )
{

    // minimum number of points required for evaluating the integral:
//  const int mpts = 1 + (int) ( k / 2 );

    // TODO -- DON'T NEED THIS MANY POINTS !!! //
    const int mpts = 3;

    // set quadrature weights and points (in interval [-1,1] //
    dTensor1 w1d(mpts), x1d(mpts);
    this->SetW1dX1d(w1d, x1d);

    // translated points, and modified 'weights'
    dTensor1 wgt(mpts), spts(mpts); 
    for( int m=1; m <= mpts; m++ )
    {
        wgt.set   (m, 0.5 * (b-a) * w1d.get(m) );
        spts.set  (m, 0.5 * ( (b+a) + (b-a) * x1d.get(m) ) );
    }

    // evaluate the polynomials at all the translated points
    dTensor2 phi(mpts, 5);
    this->evaluatePolys( spts, phi );

    // perform the integration //
    double sum = 0.0;
    for( int m=1; m <= mpts; m++ )
    { sum += wgt.get(m) * phi.get(m, k ); }

    return sum;

}

/*
void EvaluateIntegrationPoints(
    int kmax, double eta, dTensor1& ishift, const dTensor1& spts,
    dTensor3& phi, dTensor2& wgt)
{

    const int mpts = spts.getsize();

    // npts is chosen to be able to integrate a polynomial of degree 
    // 'kmax' exactly.
    int npts = 0;
    if( kmax == 1 )
    { npts = 1; }
    else if( kmax <= 3 )
    { npts = 2; }
    else
    { npts = 3; }

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points (standard before a transformation)
    //////////////////////////////////////////////////////////////////////////
    dTensor1 w1d(npts), x1d(npts);
    void SetW1dX1d(dTensor1 &w1d, dTensor1& x1d, int numpts );
    SetW1dX1d(w1d, x1d, npts );

    // compute shift index and set all the correct points
    for( int m=1; m <= mpts; m++ )
    {

        const double xm  = spts.get(m);  // in the interval [-1,1]
        const double xms = xm - 2.0*eta; // in some box
        ishift.set(m, (int) ( floor( 0.5 * (xms + 1.0 ) ) ) );

        const double b1 = xm;
        const double a1 = -1.0;
        const double b2 = xms - 2.0*ishift.get(m);
        const double a2 = -1.0;

        // local quadrature points:
        dTensor1 sptsloc1(npts), sptsloc2(npts);
        dTensor2 philoc(npts, 5);
        for( int n=1; n <= npts; n++ )
        {
            sptsloc2.set(n, 0.5*(b2-a2)*x1d.get(n) + 0.5*(b2+a2) );
            sptsloc1.set(n, 0.5*(b1-a1)*x1d.get(n) + 0.5*(b1+a1) );
        }

        // evaluate phi at these local quadrature points:
        void evaluateLegendrePolys(const dTensor1& spts, dTensor2& phi );
        evaluateLegendrePolys(sptsloc1, philoc1 );
        evaluateLegendrePolys(sptsloc2, philoc2 );

        for( int n=1; n <= npts; n++ )
        {

            wgt_tnp1.set (m, n, w1d.get(n) * (a1+1.0) / 2.0 );
            wgt_tn.set   (m, n, w1d.get(n) * (a2+1.0) / 2.0 );
            for( int k=1; k <= kmax; k++ )
            {
                phin.set( m,n,k, philoc1.get(n,k) );
                phis.set( m,n,k, philoc2.get(n,k) );
            }
        }

    }

    // modify weights and find modified quadrature points to evaluate q
    const int mpts = ishift.getsize();
    for( int m=1; m <= mpts; m++ )
    for( int n=1; n <= npts; n++ )
    {

        // weights
        wgt.set(1, m,  left_len.get(m) / 2.0 * w1d.get(m) );
        wgt.set(2, m, right_len.get(m) / 2.0 * w1d.get(m) );

        // points
        double x = x1d.get(m);
        spts.set(m,         0.5*( (bl-al)*x + (al+bl) ) );
        spts.set(morder+m,  0.5*( (br-ar)*x + (ar+br) ) );
    }

}
*/
