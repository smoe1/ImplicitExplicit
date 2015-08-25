#include "StepAdvec.h"

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
void StepAdvec(const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
        const dTensorBC2& u1, const dTensorBC2& u2, int direction,
        SL_state& sl_state)
{
    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);

    void SetBndValues(dTensorBC4& q, dTensorBC4& aux);
    SetBndValues(qold, aux );

    switch( direction )
    {
        case 1:

            StepAdvec1(1, my, dt, qold, qnew, u1, u2, sl_state);
            break;

        case 2:

            StepAdvec2(1, mx, dt, qold, qnew, u1, u2);
            break;
    }

    SetBndValues(qnew, aux );

}

// this function can be run concurrently
void StepAdvec1(int j_start, int j_end, const double& dt,
        dTensorBC4& qold, dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2, SL_state &sl_state)
{
    //
    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();

    // cell widths are uniform.  need only look at cell (1,1) for dx, dy
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // centers of grid cell (1,1)
    const double x1c = dogParamsCart2.get_xc(1);
    const double y1c = dogParamsCart2.get_yc(1);

    const int space_order = dogParams.get_space_order();
    const int mpoints1d   = space_order;
    //-----------------------------------------------------------------------//

#pragma omp parallel for
    for(int j=j_start; j<=j_end; j++)
    {

        /////////////// necessary terms for when using a source term //////////////
        dTensor1 xc_vec(2);
        dTensor1 dx_vec(2);
        dx_vec.set(1, dx);
        dx_vec.set(2, dy);
        dTensor1 tpts ( mpoints1d );
        dTensor2 xpts ( 2*mpoints1d*mpoints1d, 2 );

        // Can use this smaller psi if you want to do integration exactly //
        //dTensor2 psi  ( 2*mpoints1d*mpoints1d, meqn );
        dTensor3 psi  ( 2*mpoints1d*mpoints1d, meqn, mpoints1d );
        if( dogParams.get_source_term() > 0 )
        {

            discontCell* cell = new discontCell(kmax, 1);

            // find the quadrature time points //
            dTensor1* x1d = cell->get_x1d();
            for( int n=1; n <= mpoints1d; n++ )
            { tpts.set(n, 0.5*dt*x1d->get(n) + 0.5*(2.0*sl_state.t + dt) ); }
            delete cell;
        }
        //-----------------------------------------------------------------------//

        //legendre weight for "left/bottom" half of cell
        dTensor2 qleft(mpoints1d, kmax);		

        //legendre weights for "right/top" half of cell
        dTensor2 qright(mpoints1d, kmax);

        dTensor1 qcell(kmax);		//legendre weights for current cell
        dTensor1 scut(mpoints1d);   //location of discontinuity
        discontCell* cell = new discontCell(kmax, 1);

        for(int me=1; me<= meqn; me++)
        {
            for(int n=1; n<= mpoints1d; n++)
            {

                double lg_cfl       = u1.get(j,n)*dt/dx;
                int num_cell_shift  = (int)( floor(lg_cfl) );

                //location of discontinuity (in [-1,1])
                double scutx = 2.0 * ( lg_cfl - floor( lg_cfl ) ) - 1.0;
                if( fabs(scutx - 1.0) < 1.0e-15 )
                { scutx = 1.0; }
                scut.set(n,scutx);
            }

            cell->setScut(scut);

            // limit this entire 'row' if using the limiter:
            if(dogParams.using_moment_limiter())
            {

                dTensor2* spts = cell->get_limited_spts();

                // ---------------------------------------------------------- //
                // sample basis at all points where we want to check solution //
                // ---------------------------------------------------------- //
                void SamplePhiAtPositivePoints(const int& space_order, 
                        const dTensor2& spts, dTensor2& phi);
                dTensor2 phi(spts->getsize(1), kmax);
                SamplePhiAtPositivePoints(space_order, *spts, phi);

                for(int i=1; i<=mx; i++)
                {

                    // -------------------------------------------------------------- //
                    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
                    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
                    // -------------------------------------------------------------- //
                    double m = 0.0;
                    for(int mp=1; mp <= phi.getsize(1); mp++)
                    {
                        // evaluate q at spts(mp) //
                        double qnow = 0.0;
                        for( int k=1; k <= kmax; k++ )
                        {
                            qnow += qold.get(i,j,me,k) * phi.get(mp,k);
                        }

                        // note: for most computations, int_psi = 0 //
                        m = Min(m, qnow );
                    }

                    const double eps = 1e-12;
                    double Q1 = qold.get(i,j,me,1);

                    double theta = 0.0;
                    if( fabs( Q1 - m ) < 1.0e-14 ) { theta = 0.0; }
                    else{ theta = Min( 1.0, (Q1 - eps) / (Q1 - m) ); }

                    // limit q //
                    for( int k=2; k <= kmax; k++ )
                    {
                        qold.set(i,j,me,k, qold.get(i,j,me,k) * theta );
                    }

                }
            }

            for(int i=1; i<=mx; i++)
            {

                for( int n=1; n<= mpoints1d; n++)
                {

                    double speedx = u1.get(j,n);
                    int num_cell_shift = (int)( floor(speedx*dt/dx) );

                    /////////////////////////////////////////////////////
                    // Set ql and qr
                    /////////////////////////////////////////////////////
                    //io = 'old' cell that provides data for current cell i
                    int io = (int)(i - num_cell_shift);
                    assert_ge( io, 1-mbc  );
                    assert_le( io, mx+mbc );

                    // periodic boundary condtions would be enforced if we
                    // used iMod here:
                    //io = iMod((io-1),mx) + 1;	

                    // use ghost cells for the data!
                    for(int k= 1; k<=kmax; k++)
                    {
                        qright.set(n, k, qold.get(io,   j, me, k) );
                        qleft.set (n, k, qold.get(io-1, j, me, k) );
                    }//end of loop over each polynomial

                }

                //create a discontinuous cell and project onto legendre basis
                cell->project(qleft, qright, qcell);

                //save legendre weights into qnew
                for(int k=1; k<=kmax; k++)
                {
                    qnew.set(i,j, me, k, qcell.get(k));
                }
                // if a source term exists //
                if( dogParams.get_source_term() > 0 )
                {

                    xc_vec.set(1, dogParamsCart2.get_xc(i) );  
                    xc_vec.set(2, dogParamsCart2.get_yc(j) );

                    cell->set_xpts( &xc_vec, &dx_vec, &xpts );

                    /////// This section uses quadrature in time ////
                    void IntegratePsiBackward( double tnew, const dTensor2 &xpts, 
                                        const dTensor1 &tpts, dTensor3 &psi);
                    IntegratePsiBackward( sl_state.t + dt, xpts, tpts, psi );
                    //////////////////////////////////////////////////

                    dTensor2 qextra(meqn, kmax);
                    cell->project_source_w_time( psi, qextra );

                    //save legendre weights into qnew
                    for(int k=1; k<=kmax; k++)
                    {
                        //////////////////////////////////////////////////
                        /////// This section uses quadrature in time ////
                        qnew.set(i,j,me,k, qnew.get(i,j,me,k) + 
                                            0.5*dt*qextra.get(me, k));
                    }
                }
            }
        }

        delete cell;

    }

}

void StepAdvec2(int i_start, int i_end, const double& dt, 
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

    const int space_order = dogParams.get_space_order();
    const int mpoints1d   = space_order;
    //-----------------------------------------------------------------------//

    //loop over every cell
#pragma omp parallel for
    for(int i=i_start; i<=i_end; i++)
    {

        // Placed here because they need to be private to each thread -- //
        discontCell* cell = new discontCell(kmax, 2);

        dTensor1 qcell(kmax);		//legendre weights for current cell
        dTensor1 scut(mpoints1d);   //location of discontinuity

        //legendre weight for "left/bottom" half of cell
        dTensor2 qleft(mpoints1d, kmax);		

        //legendre weights for "right/top" half of cell
        dTensor2 qright(mpoints1d, kmax);
        ///////////////////////////////////////////////////////////////////

        for(int me=1; me<= meqn; me++)
        {
            for( int n=1; n<= mpoints1d; n++)
            {

                double lg_cfl       = u2.get(i,n)*dt/dy;
                int num_cell_shift  = (int)( floor(lg_cfl) );

                //location of discontinuity (in [-1,1])
                double scuty = 2.0 * ( lg_cfl - floor( lg_cfl ) ) - 1.0;
                if( fabs(scuty - 1.0) < 1.0e-15 )
                { scuty = 1.0; }

                scut.set(n,scuty);

            }

            // this is the same for every collumn, and is an expensive
            // call
            cell->setScut(scut);

            // limit this entire 'collumn' if using the limiter:
            if(dogParams.using_moment_limiter())
            {

                dTensor2* spts = cell->get_limited_spts();
                // ---------------------------------------------------------- //
                // sample basis at all points where we want to check solution //
                // ---------------------------------------------------------- //
                void SamplePhiAtPositivePoints(const int& space_order, 
                        const dTensor2& spts, dTensor2& phi);
                dTensor2 phi(spts->getsize(1), kmax);
                SamplePhiAtPositivePoints(space_order, *spts, phi);

                for(int j=1; j<=my; j++)
                {
                    // -------------------------------------------------------------- //
                    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
                    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
                    // -------------------------------------------------------------- //
                    double m = 0.0;
                    for(int mp=1; mp <= phi.getsize(1); mp++)
                    {
                        // evaluate q at spts(mp) //
                        double qnow = 0.0;
                        for( int k=1; k <= kmax; k++ )
                        {
                            qnow += qold.get(i,j,me,k) * phi.get(mp,k);
                        }
                        m = Min(m, qnow);
                    }

                    const double eps = 1e-12;
                    double theta = 0.0;
                    double Q1 = qold.get(i,j,me,1);

                    if( fabs( Q1 - m ) < 1.0e-14 ) { theta = 0.0; }
                    else{ theta = Min( 1.0, (Q1 - eps) / (Q1 - m) ); }

                    // limit q //
                    for( int k=2; k <= kmax; k++ )
                    {
                        qold.set(i,j,me,k, qold.get(i,j,me,k) * theta );
                    }

                }
            }

            for(int j=1; j<=my; j++)
            {
                for( int n=1; n<= mpoints1d; n++)
                {

                    /////////////////////////////////////////////////////
                    // Set ql and qr
                    /////////////////////////////////////////////////////					
                    //jo = 'old' cell that provides data for current cell j
                    double speedy = u2.get(i, n);
                    int num_cell_shift = (int) ( floor(speedy * dt/ dy) );
                    int jo = (int)(j - num_cell_shift);

                    //save legendre weight from qold using
                    //using periodic boundary conditions
                    jo = iMod((jo-1),my) + 1;

                    for(int k= 1; k<=kmax; k++)
                    {

                        qright.set(n, k, qold.get(i, jo, me, k) );
                        if( jo - 1 == 0)
                        {
                            qleft.set(n, k, qold.get(i, my, me, k) );
                        }else
                        {
                            qleft.set(n, k, qold.get(i, jo-1, me, k) );
                        }
                    }//end of loop over each polynomial
                }

                //create a discontinuous cell and project onto legendre basis
                cell->project(qleft, qright, qcell);

                //save legendre weights into qnew
                for(int k=1; k<=kmax; k++)
                {
                    qnew.set(i,j, me, k, qcell.get(k));
                }

            }
        }

        delete cell;

    }

}

/*
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
    void init_projection_tensor( dTensor3& a );
    init_projection_tensor( a );

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
*/

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
void init_projection_tensor( dTensor3& a )
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
