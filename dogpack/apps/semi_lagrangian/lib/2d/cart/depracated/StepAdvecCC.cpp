#include "StepAdvec.h"

///////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation:
//
//            q_t + u q_x + v q_y = 0 
//
// This routine was never finished, and is inoperable.
//
// The idea was to copy the single step routine placed in the 1D library, but it
// will take a bit more work to make it fully operable. (-DS)
//
///////////////////////////////////////////////////////////////////////////////
void StepAdvecCC(double dt, const dTensorBC4& qold, const dTensor2& speeds,
    dTensorBC4& qnew )
{

    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int space_order = dogParams.get_space_order();

    // number of points used for integration on each cell
    const int numpts = 2*space_order*space_order;

    // cell widths are uniform.  need only look at cell (1,1)
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double x1c = dogParamsCart2.get_xc(1); //center of grid cell (1,1)
    const double y1c = dogParamsCart2.get_yc(1);
    double scutx = 0.0;
    double scuty = 0.0;
    double xcut = 0.0;
    double ycut = 0.0;
    //-----------------------------------------------------------------------//

    const double s_area = 4.0;

    // Fractional eta values - 
    //      frac_eta( 1:meqn, 1 ) = x-direction, 
    //      frac_eta( 1:meqn, 2 ) = y direction
    dTensor2 frac_eta( meqn, 2 );
    for( int m=1; m <= meqn; m++ )
    {
        double large_eta = speeds.get(m,1) * dt / dx;
        frac_eta.set( m, 1, (double)( large_eta - floor( large_eta ) ) );

        large_eta = speeds.get(m,2) * dt / dy;
        frac_eta.set( m, 2, (double)( large_eta - floor( large_eta ) ) );
    }


    // integration points and weights
    dTensor2  wgt(meqn, numpts);
    dTensor3 spts(meqn, numpts, 2);
    void setIntegrationPoints(dTensor2& frac_eta, dTensor3& spts, dTensor2& wgt);
    setIntegrationPoints(frac_eta, spts, wgt);

    // 'old' integration points:
    dTensor2 spts_old ( numpts, 2 );
    iTensor1 ishift   ( numpts );
    iTensor1 jshift   ( numpts );
    void translateXi( 
        double eta, const dTensor1& xi, dTensor1 &xinew, iTensor1& ishift );
    translateXi( large_eta, spts, spts_old, ishift );

    // Legendre Basis functions, their derivatives, evaluated at all the
    // necessary poitns

    dTensor2 phi(meqn, numpts, kmax);
    dTensor2 phi_old(numpts, kmax);
    void evaluateLegendrePolys( const dTensor2& spts, dTensor3& phi );
    evaluateLegendrePolys( spts,     phi      );
    evaluateLegendrePolys( spts_old, phi_old  );

    // Check and apply the limiter if necessary --------------- //
    if(dogParams.using_moment_limiter())
    {

        for(int i=1; i<=mx; i++)
        for(int j=1; j<=my; j++)
        for(int me=1; me <= meqn; me++ )
        {

            // -------------------------------------------------------------- //
            // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
            // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
            // -------------------------------------------------------------- //
            double m = qold.get(i,j,me,1);
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

            double Q1 = qold.get(i,j,me,1);

            double theta = 0.0;
            if( fabs( Q1 - m ) < EPSILON ) { theta = 0.0; }
            else{ theta = Min( 1.0, (Q1 - EPSILON) / (Q1 - m) ); }

            // limit q //
            for( int k=2; k <= kmax; k++ )
            {
                qold.set(i,j,me,k, qold.get(i,j,me,k) * theta );
            }

        }

    }
    // End of limiter //
    // //////////////////////////////////////////// //

    for(int i=1; i<=mx; i++)
    for(int j=1; j<=my; j++)
    for(int me=1; me <= meqn; me++ )
    {

        for( int k=1; k <= kmax; k++)
        {
            double sum = 0.0;
            for( int m=1; m <= numpts; m++ )
            {

                // hard coded periodicity is enforced here!
                int io = (int)( i + ishift.get(m) );
                io = iMod((io-1), mx)+1;

                int jo = (int)( j + jshift.get(m) );
                jo = iMod((jo-1), my)+1;

                // evaluate q at the old point:
                double qval = 0.0;
                for( int kq=1; kq <= kmax; kq++ )
                {
                    qval += qold.get(io, jo, me, kq) * phi_old.get(m, kq);
                }
                sum += wgt.get(m) * qval * phi.get(m, k);
            }
            qnew.set(i, j, me, k, sum / s_area );
        }

    }

    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    SetBndValues(node, auxvals, qnew);

}

void setIntegrationPoints(dTensor2& frac_eta, dTensor2& spts, dTensor2& wgt)
{

    const int morder = dogParams.get_space_order();

    const int meqn   = wgt.getsize(1); // == frac_eta.getsize(1)
    const int numpts = wgt.getsize(2);

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points (standard before a transformation)
    //////////////////////////////////////////////////////////////////////////

    dTensor1 w1d(morder), x1d(morder);
    switch ( morder )
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

    for( int me=1; me <= meqn; me++ )
    {

        // location of discontinuity and width of left and right hand side
        const double scut      = -1.0 + 2.0*frac_eta.get(me,1);
        const double left_len  = scut - (-1.0);
        const double right_len = 1.0  - scut;

        // left half of the integral
        const double bl = scut;
        const double al = -1.0;

        // right half of the integral
        const double br = 1.0;
        const double ar = scut;

        dTensor1 wgt_x( 2*morder );

        // modify weights and find modified quadrature points to evaluate q
        for( int m=1; m <= morder; m++ )
        {

            // weights
            wgt_x.set(m,           left_len / 2.0 * w1d.get(m) );
            wgt_x.set(m + morder, right_len / 2.0 * w1d.get(m) );

            // points
            double x = x1d.get(m);
            spts.set(m,         0.5*( (bl-al)*x + (al+bl) ) );
            spts.set(morder+m,  0.5*( (br-ar)*x + (ar+br) ) );
        }

    }

}// end of function StepAdvecConstCoeff.cpp
