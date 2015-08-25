#include "StepAdvec.h"

void set1DGaussPoints(dTensor1* w1d, dTensor1* x1d)
{

    assert( w1d->getsize() == x1d->getsize() );

    // ---------------------------------
    // Set 1D quadrature weights and points
    // ---------------------------------
    switch ( w1d->getsize() )
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

void evaluate_phi2d( int space_order, const dTensor2& spts, dTensor2& phi)
{
    const int mpoints = spts.getsize(1);

    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi  = spts.get(m,1);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;

        const double eta = spts.get(m,2);
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
