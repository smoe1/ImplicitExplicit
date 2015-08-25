#include "StepAdvecHybrid.h"

///////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation:
//
//            q_t + u(y) q_x + v(x) q_y = 0 
//
// forward in time using an operator split method.  Time integration is
// only performed in one direction - x or y.
//
// Parameters
//
// dt - time step taken in each equation
// qold(mx,my,meqn,kmax) - q before translation
// qnew - q after being translated and projected onto the legendre basis
// advec_vel(mx, my, mpoints1d, 2) - advection velocity.
// direction - direction for velocity to take place
//			== 1: travel in x-direction (i.e. v = 0)
//			== 2: travel in y-direction (i.e. u = 0)
//
//  u1(my, mpoints1d) = u(y) - the 1d gauss point values for speeds
//  u2(mx, mpoints1d) = v(x) - the 1d gauss point values for speeds
//
///////////////////////////////////////////////////////////////////////////////

// Wrapper function to split StepAdvec into chunks
void StepAdvecHybrid(const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
        const dTensorBC2& u1, const dTensorBC2& u2, int direction,
        SL_state& sl_state)
{
    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();

    switch( direction )
    {
        case 1:

            StepAdvec1lines(1, my, dt, qold, qnew, u1, u2, sl_state);
            break;

        case 2:

            StepAdvec2lines(1, mx, dt, qold, qnew, u1, u2);

            break;
    }

    void SetBndValues(dTensorBC4& q, dTensorBC4& aux);
    SetBndValues(aux, qnew);

}

// Experimental methods where Q is projected onto `lines', followed by many 1D
// solves
void StepAdvec1lines(int j_start, int j_end, const double& dt, 
        dTensorBC4& qold,  // limiter will modify qold
        dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2, SL_state& sl_state)
{

    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();

    // cell widths are uniform.  need only look at cell (1,1)
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double x1c = dogParamsCart2.get_xc(1); // center of grid cell (1,1)
    const double y1c = dogParamsCart2.get_yc(1);

    const double ylow = dogParamsCart2.get_ylow();
    const double xlow = dogParamsCart2.get_xlow();

    const int space_order = dogParams.get_space_order();
    const int mpoints1d   = space_order;
    const int kmax1d = space_order;
    //-----------------------------------------------------------------------//

    dTensor3 a(kmax1d,kmax1d,kmax);
    void init_projection_tensor1( dTensor3& a );
    init_projection_tensor1( a );

    dTensor1 w1d( kmax1d );
    dTensor1 x1d( kmax1d );
    void set1DGaussPoints(dTensor1* w1d, dTensor1* x1d);
    set1DGaussPoints(&w1d, &x1d);


    //loop over every cell
#pragma omp parallel for
    for(int j=j_start; j<=j_end; j++)
    {

        dTensor1 speeds(mpoints1d);
        for( int me=1; me <= mpoints1d; me++ )
        { speeds.set(me, u1.get(j,me) ); }

        // Project onto 1D coefficients
        dTensorBC3 q1d   ( mx, kmax1d, kmax1d, mbc );
        dTensorBC3 q1dnew( mx, kmax1d, kmax1d, mbc );

        // Project current slice onto 1D problems:
        for(int i=1; i<=mx; i++)
        {
            for( int m=1; m <= kmax1d; m++ )
            {
                for( int k1=1; k1 <= kmax1d; k1++ )
                {
                    double tmp = 0.;
                    for( int k2=1; k2 <= kmax; k2++ )
                    {
                        tmp += qold.get(i,j,1,k2) * a.get(m,k1,k2);
                    }
                    q1d.set(i,m,k1, tmp);
                }
            }

        }

        // Solve all of the 1D problems for this collumn
        RKSolve1D( sl_state.t, dt, xlow, dx, speeds, q1d, q1dnew );

        // Integrate up to the 2D weights using each slice
        // q1d( m, k1 ) = sum_{k2=1:kmax2d} 0.5 * w1d(m) * q2d(:,k2) * a( m, k1, k2 ) 
        for( int i=1; i <= mx; i++ )
        {
            for( int k2=1; k2 <= kmax; k2++ )
            {
                double tmp = 0.;
                for( int m=1; m <= kmax1d; m++ )
                {
                    for( int k1=1; k1 <= kmax1d; k1++ )
                    {
                        tmp += 0.5 * w1d.get(m) * a.get(m,k1,k2) * q1dnew.get(i,m,k1);
                    }
                }
                qnew.set(i, j, 1, k2, tmp );
            }

        }


    } 

}


// Experimental methods where Q is projected onto `lines', followed by many 1D
// solves
void StepAdvec2lines(int i_start, int i_end, const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2)
{

    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();

    // cell widths are uniform.  need only look at cell (1,1)
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double x1c = dogParamsCart2.get_xc(1); // center of grid cell (1,1)
    const double y1c = dogParamsCart2.get_yc(1);

    const double ylow = dogParamsCart2.get_ylow();

    const int space_order = dogParams.get_space_order();
    const int mpoints1d   = space_order;
    const int kmax1d = space_order;
    //-----------------------------------------------------------------------//

    dTensor3 a(kmax1d,kmax1d,kmax);
    void init_projection_tensor2( dTensor3& a );
    init_projection_tensor2( a );

    dTensor1 w1d( kmax1d );
    dTensor1 x1d( kmax1d );
    void set1DGaussPoints(dTensor1* w1d, dTensor1* x1d);
    set1DGaussPoints(&w1d, &x1d);


    //loop over every cell
#pragma omp parallel for
    for(int i=i_start; i<=i_end; i++)
    {

        dTensor1 speeds(mpoints1d);
        for( int me=1; me <= mpoints1d; me++ )
        { speeds.set(me, u2.get(i,me) ); }

        // Project onto 1D coefficients
        dTensorBC3 q1d   ( my, kmax1d, kmax1d, mbc );
        dTensorBC3 q1dnew( my, kmax1d, kmax1d, mbc );

        // Project current slice onto 1D problems:
        for( int j=1; j <= my; j++ )
        {
            for( int m=1; m <= kmax1d; m++ )
            {
                for( int k1=1; k1 <= kmax1d; k1++ )
                {
                    double tmp = 0.;
                    for( int k2=1; k2 <= kmax; k2++ )
                    {
                        tmp += qold.get(i,j,1,k2) * a.get(m,k1,k2);
                    }
                    q1d.set(j,m,k1, tmp);
                }
            }

        }

        // Solve all of the 1D problems for this collumn
        StepAdvec1D( dt, ylow, dy, speeds, q1d, q1dnew );

        // Integrate up to the 2D weights using each slice
        // q1d( m, k1 ) = sum_{k2=1:kmax2d} 0.5 * w1d(m) * q2d(:,k2) * a( m, k1, k2 ) 
        for( int j=1; j <= my; j++ )
        {
            for( int k2=1; k2 <= kmax; k2++ )
            {
                double tmp = 0.;
                for( int m=1; m <= kmax1d; m++ )
                {
                    for( int k1=1; k1 <= kmax1d; k1++ )
                    {
                        tmp += 0.5 * w1d.get(m) * a.get(m,k1,k2) * q1dnew.get(j,m,k1);
                    }
                }
                qnew.set(i,j,1,k2, tmp );
            }

        }


    } 

}

// Construct projection tensor a.  This tensor satisfies the relation:
//
// q1d( m, k1 ) = sum_{k2=1:kmax2d} q2d(:,k2) * a( m, k1, k2 ) 
//
// a can is computed via:
//
//      a(m,k1,k2) = 
//
//          1 / 2 \int_{-1}^1 \phi1d^(k1)(xi) \phi2d^(k2)(xi,eta_m)\ dxi
//
//      m == m^th `row'
//     k1 == 1D quadrature weight
//     k2 == 2D quadrature weight
void init_projection_tensor1( dTensor3& a )
{

    const int kmax1d  = a.getsize(1);
    const int numpts  = a.getsize(2);
    const int kmax2d  = a.getsize(3);

    assert( kmax1d == numpts );

    dTensor1 w1d(kmax1d);
    dTensor1 spts(kmax1d);
    void set1DGaussPoints(dTensor1* w1d, dTensor1* x1d);
    set1DGaussPoints( &w1d, &spts);

    dTensor2 phi(numpts, 5);
    void evaluateLegendrePolys1D( const dTensor1& spts, dTensor2& phi );
    evaluateLegendrePolys1D( spts, phi );

    dTensor2 spts2d(numpts, 2);
    dTensor2  phi2d(kmax2d, numpts );
    for( int m=1; m <= numpts; m++ )
    {

        // 2D points for current row:
        for( int n=1; n <= numpts; n++ )
        {
            spts2d.set(n,1, spts.get( n ) );
            spts2d.set(n,2, spts.get( m ) );
        }
        
        // sample 2D basis function at each point:
        void evaluate_phi2d( int space_order, const dTensor2& spts, dTensor2& phi);
        evaluate_phi2d( kmax1d, spts2d, phi2d );

        // compute \int phi2d( xi_m, eta ) * phi1d( eta )\ deta
        for( int k1=1; k1 <= kmax1d; k1++)
        for( int k2=1; k2 <= kmax2d; k2++)
        {
            double tmp = 0.;
            for( int n=1; n <= numpts; n++ )
            {
                tmp += w1d.get(n) * phi2d.get(k2,n) * phi.get(n,k1);
            }

            // q1d( m, k1 ) = sum_{k2=1:kmax2d} q2d(:,k2) * a( m, k1, k2 ) 
            a.set(m,k1,k2, 0.5 * tmp );

        }

 
    }

}

// Construct projection tensor a.  This tensor satisfies the relation:
//
// q1d( m, k1 ) = sum_{k2=1:kmax2d} q2d(:,k2) * a( m, k1, k2 ) 
//
// a can is computed via:
//
//      a(m,k1,k2) = 
//
//          1 / 2 \int_{-1}^1 \phi1d^(k1)(eta) \phi2d^(k2)(xi_m,eta)\ deta
//
//      m == m^th `row'
//     k1 == 1D quadrature weight
//     k2 == 2D quadrature weight
void init_projection_tensor2( dTensor3& a )
{

    const int kmax1d  = a.getsize(1);
    const int numpts  = a.getsize(2);
    const int kmax2d  = a.getsize(3);

    assert( kmax1d == numpts );

    dTensor1 w1d(kmax1d);
    dTensor1 spts(kmax1d);
    void set1DGaussPoints(dTensor1* w1d, dTensor1* x1d);
    set1DGaussPoints( &w1d, &spts);

    dTensor2 phi(numpts, 5);
    void evaluateLegendrePolys1D( const dTensor1& spts, dTensor2& phi );
    evaluateLegendrePolys1D( spts, phi );

    dTensor2 spts2d(numpts, 2);
    dTensor2  phi2d(kmax2d, numpts );
    for( int m=1; m <= numpts; m++ )
    {

        // 2D points for current row:
        for( int n=1; n <= numpts; n++ )
        {
            spts2d.set(n,1, spts.get( m ) );
            spts2d.set(n,2, spts.get( n ) );
        }
        
        // sample 2D basis function at each point:
        void evaluate_phi2d( int space_order, const dTensor2& spts, dTensor2& phi);
        evaluate_phi2d( kmax1d, spts2d, phi2d );

        // compute \int phi2d( xi_m, eta ) * phi1d( eta )\ deta
        for( int k1=1; k1 <= kmax1d; k1++)
        for( int k2=1; k2 <= kmax2d; k2++)
        {
            double tmp = 0.;
            for( int n=1; n <= numpts; n++ )
            {
                tmp += w1d.get(n) * phi2d.get(k2,n) * phi.get(n,k1);
            }

            // q1d( m, k1 ) = sum_{k2=1:kmax2d} q2d(:,k2) * a( m, k1, k2 ) 
            a.set(m,k1,k2, 0.5 * tmp );

        }

 
    }

}
