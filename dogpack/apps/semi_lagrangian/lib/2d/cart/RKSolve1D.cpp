#include "RKSolve1D.h"

//////////////////////////////////////////////////////////////////////////////
//
// Function for stepping 1D advection equation forward in time.
//
//////////////////////////////////////////////////////////////////////////////
void RKSolve1D(
    double tn,
    double dt, double xlow, double dx,
    const dTensor1& speeds,
    dTensorBC3& qold,
    dTensorBC3& qnew)
{

    // TODO - this isn't the case for generic velocities, and is currently
    // poorly coded!!!
    //
    // Note: for VP, speeds(1:meqn) = velocity values ...
    //
    //

    //-local parameters -----------------------------------------------------//
    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

    double smax = 0.;
    for( int m=1; m <= meqn; m++ )
    { smax = Max( smax, fabs( speeds.get(m) ) ); }

    // TODO - this should be a user supplied parameter, as opposed to a hard
    // coded constant here!
    const double cfl = 0.1;

    // h = micro substep. TODO - check that these things actually add up to the
    // correct total time ...
    const int N = iMax( 1, (int)( ceil( ( fabs(dt) * smax ) / ( dx * cfl ) ) ) );
    double h = double( dt / N );

    double t = RKStep1D( tn, h, xlow, dx, speeds, qold, qnew); 
    for( int n=2; n <= N; n++ )
    {
        // CopyQ( qnew, qold );  WHY DOESN'T COPYQ WORK HERE? -DS
        // Construct Full Update ...
        for( int i=1-mbc; i <= mx+mbc ; i++ )
        for( int m=1; m <= meqn; m++ )
        for( int k=1; k <= kmax; k++ )
        { qold.set(i,m,k, qnew.get(i,m,k) ); }

        t = RKStep1D( t, h, xlow, dx, speeds, qold, qnew); 
    }

    // TODO - rip this line out ...
//  if( fabs( t - (tn+dt) ) > EPSILON )
//  {
//      double err = fabs( t - (tn+dt) );
//      printf( "Warning; t-tnp1 != 0;  the error = %2.15e\n", err );
//  }

}

// the only purpose of this function is for subcycling.  the dt passed into here
// represents a small dt, whereas the dt passed into RKSolve1D represents a
// macro time step
double RKStep1D( double tn,
    double dt, double xlow, double dx,
    const dTensor1& speeds,
    dTensorBC3& qold,  // Set Boundary conditions modifies qold ...
    dTensorBC3& qnew)
{

    //-local parameters -----------------------------------------------------//
    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

   
    // Memory for each of the sub stages TODO 
    dTensorBC3 k1(mx,meqn,kmax,mbc);
    dTensorBC3 k2(mx,meqn,kmax,mbc);
    dTensorBC3 k3(mx,meqn,kmax,mbc);
    dTensorBC3 k4(mx,meqn,kmax,mbc);

    double t = tn;
    const double sgn_dt = sgn( dt );

    // Stage 1:
    ConstructL1D( t, sgn_dt, xlow, dx, speeds, qold, k1 );
    t = EulerStep( tn, 0.5*dt, qold, k1, qnew );

    // Stage 2:
    ConstructL1D( t, sgn_dt, xlow, dx, speeds, qnew, k2 );
    t = EulerStep( tn, 0.5*dt, qold, k2, qnew );

    // Stage 3:
    ConstructL1D( t, sgn_dt, xlow, dx, speeds, qnew, k3 );
    t = EulerStep( tn, dt, qold, k3, qnew );

    // Stage 4:
    ConstructL1D( t, sgn_dt, xlow, dx, speeds, qnew, k4 );

    // Construct Full Update ...
    for( int i=1; i <= mx ; i++ )
    for( int m=1; m <= meqn; m++ )
    for( int k=1; k <= kmax; k++ )
    {
        double tmp = 0.;

        tmp +=     k1.get(i,m,k);
        tmp += 2.0*k2.get(i,m,k);
        tmp += 2.0*k3.get(i,m,k);
        tmp +=     k4.get(i,m,k);

        qnew.set(i,m,k, qold.get(i,m,k) + dt*tmp/6.0 );

    }
    /////////////////////////////////////////////////////


    // --------------------------------------------------------

    SetPeriodicBndy1D( qnew );
    assert( fabs(tn+dt - t ) < 1e-12 );
    return tn+dt;

}

//
// Ripped, and gutted from the 1D library ... -DS
//
// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
void ConstructL1D(
        double t,
        const double sgn_dt,
        const double xlow, const double dx,
        const dTensor1& speeds,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC3& Lstar )
{

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    //const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // Flux values, interior integrals and source term, respectively.
    dTensorBC2 Fm(melems,meqn,mbc);
    dTensorBC2 Fp(melems,meqn,mbc);
    dTensorBC3  N(melems,meqn,kmax,mbc);
    dTensorBC3 Psi(melems,meqn,kmax,mbc);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetPeriodicBndy1D( q );

    // Loop over interior edges and solve Riemann problems
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {
        dTensor1 Ql(meqn);
        dTensor1 Qr(meqn);
        dTensor1 Auxl(meqn);
        dTensor1 Auxr(meqn);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );

            for (int k=1; k<=kmax; k++)
            {
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *q.get(i-1,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            }
        }

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {

            // Advection equation: forced upwind flux chosen here:
            // TODO - accomodate negative velocities ...

            // TODO - what's the correct thing to place here!?

            double qv = 0.;
            if( sgn_dt * speeds.get(m) > 0 )
            { qv = speeds.get(m) * Ql.get(m); }
            else
            { qv = speeds.get(m) * Qr.get(m); }

            Fm.set(i,  m, qv );
            Fp.set(i-1,  m, qv );

        }
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( f(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    if ( kmax > 1 )
    {  L2ProjectAdvection(1,1-mbc,melems+mbc,dx,q,speeds,N);  }

    else
    { N.setall(0.); }

    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( dogParams.get_source_term() > 0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2ProjectSource(t, 1-mbc, melems+mbc, dx, 
            speeds,  Psi, &HybridSourceTermFunc1D);
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    if ( dogParams.get_source_term() == 0 )
    {

        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
        for (int m=1; m<=meqn; m++)
        for (int k=1; k<=kmax; k++)
        {
            double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;

            Lstar.set(i,m,k, tmp );	      
        }
    }
    else  // With Source Term
    {

        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)	    
                for (int k=1; k<=kmax; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
                        + Psi.get(i,m,k);

                    Lstar.set(i,m,k, tmp );
                }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
//  LstarExtra(aux,q,Lstar);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part VI: artificial viscosity limiter
    // ---------------------------------------------------------
//  if (dogParams.get_space_order()>1)
//  {
//      if (dogParams.using_viscosity_limiter())
//      {  ArtificialViscosity(aux,q,Lstar);  }
//  }
    // ---------------------------------------------------------
}

// Single Eulerian Step
double EulerStep( double t, double dt, 
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& qnew )
{

    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells

    for( int i=1; i <= mx ; i++ )
    for( int m=1; m <= meqn; m++ )
    for( int k=1; k <= kmax; k++ )
    {
        double tmp = qold.get(i,m,k) + dt*Lstar.get(i,m,k);
        qnew.set(i,m,k, tmp);
    }

    return t+dt;

}

// Forced periodic boundary conditions
void SetPeriodicBndy1D( dTensorBC3 &qnew )
{

    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

    for( int mb=1; mb <= mbc; mb++ )
    for( int m=1; m <= meqn; m++ )
    for( int k=1; k <= kmax; k++ )
    {
        qnew.set(1-mb,  m, k, qnew.get(mx-mb+1, m , k ) );
        qnew.set(mx+mb, m, k, qnew.get(mb,      m , k ) );
    }

}

// ---------------------------------------------------------------------
//
// Ripped, and gutted from the 1D library ... -DS
//
// All-purpose routine for computing the L2-projection
// of various functions onto:
//     mopt==0:   the Legendre basis
//     mopt==1:   the derivatives of Legendre basis
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC3 auxin(1-mbc:mnodes+mbc,maux,mpoints)
//           dTensorBC3   qin(1-mbc:mnodes+mbc,meqn,mpoints)
//           dTensorBC3  Fout(1-mbc:mnodes+mbc,mlength,mpoints)
// ---------------------------------------------------------------------

void L2ProjectAdvection(int mopt, int istart, int iend,
        double dx,
        const dTensorBC3& qin,
        const dTensor1& speeds,  
        dTensorBC3& Fout
    )
{    

    const int KMAX = 6;

    const int meqn    = qin.getsize(2);
    const int mlength = Fout.getsize(2); // i.e., number of equations ...
    const int mpoints = Fout.getsize(3);  

    dTensor1 wgt(mpoints), spts(mpoints);
    dTensor2 phi(mpoints,KMAX), phi_x(mpoints,KMAX);

    // ---------------------------------------------
    // Check for trivial case in the case of mopt==1
    // ---------------------------------------------
    if ( mpoints == mopt )
    { Fout.setall(0.); }
    else
    {
        // ---------------------------------
        // Set quadrature weights and points
        // ---------------------------------
        switch ( (mpoints-mopt) )
        {
            case 1:
                wgt.set(1,  2.0e0 );
                spts.set(1, 0.0e0 );

                break;

            case 2:
                wgt.set(1,   1.0 );
                wgt.set(2,   1.0 );

                spts.set(1, -1.0/sq3 );
                spts.set(2,  1.0/sq3 );

                break;

            case 3:
                wgt.set(1, 5.0e0/9.0e0 );
                wgt.set(2, 8.0e0/9.0e0 );
                wgt.set(3, 5.0e0/9.0e0 );

                spts.set(1, -sq3/sq5 );
                spts.set(2,  0.0e0 );
                spts.set(3,  sq3/sq5 );

                break;

            case 4:
                wgt.set(1, (18.0 - sq30)/36.0 );
                wgt.set(2, (18.0 + sq30)/36.0 );
                wgt.set(3, (18.0 + sq30)/36.0 );
                wgt.set(4, (18.0 - sq30)/36.0 );

                spts.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
                spts.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
                spts.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
                spts.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

                break;

            case 5:
                wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
                wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
                wgt.set(3,  128.0/225.0 );
                wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
                wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

                spts.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
                spts.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
                spts.set(3,  0.0 );
                spts.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
                spts.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

                break;

            case 6:

                wgt.set(1,  0.171324492379170);
                wgt.set(2,  0.360761573048139);
                wgt.set(3,  0.467913934572691);
                wgt.set(4,  0.467913934572691);
                wgt.set(5,  0.360761573048139);
                wgt.set(6,  0.171324492379170);


                spts.set(1, -0.932469514203152);
                spts.set(2, -0.661209386466265);
                spts.set(3, -0.238619186083197);
                spts.set(4, -spts.get(3) );
                spts.set(5, -spts.get(2) );
                spts.set(6, -spts.get(1) );

                break;

        }

        // Loop over each quadrature point to construct Legendre polys
        for (int m=1; m<=(mpoints-mopt); m++)
        {
            const double xi  = spts.get(m);
            const double xi2 = xi*xi;
            const double xi3 = xi2*xi;
            const double xi4 = xi3*xi;
            const double xi5 = xi4*xi;

            // Legendre basis functions at each grid point
            phi.set( m,1, 1.0 );
            phi.set( m,2, sq3*xi );
            phi.set( m,3, 0.5*sq5*( 3.0*(xi2) - 1.0 ) );
            phi.set( m,4, 0.5*sq7*xi
                    *(5.0*(xi2) - 3.0) );
            phi.set( m,5, (105.0/8.0)*(xi4) 
                    - (45.0/4.0)*(xi2) + (9.0/8.0) );
            phi.set( m, 6, (63.0/8.0)*sq11 * ( (xi5) - (10.0/9.0)*(xi3)+ (5.0/21.0)*xi ) );

            // 1st derivative of Legendre basis functions at each grid point
            phi_x.set( m, 1, 0.0 );
            phi_x.set( m, 2, 2.0*sq3/dx );
            phi_x.set( m, 3, 6.0*sq5*xi/dx );
            phi_x.set( m, 4, 3.0*sq7*(5.0*(xi2)-1.0)/dx );
            phi_x.set( m, 5, 15.0*xi* (7.0*(xi2)-3.0)/dx );
            phi_x.set( m, 6, (2.0/dx)*(63.0/8.0)*sq11*(5.0*(xi4)-(10.0/3.0)*(xi2)+5.0/21.0) );

        }

        // ----------------------------------
        // Loop over all elements of interest
        // ----------------------------------    
        const double s_area = 2.0;
        for (int i=istart; i<=iend; i++)
        {

            double xc = dogParamsCart2.get_xlow() + (double(i)-0.5)*dx;

            // each of these three items needs to be private to each thread ..
            dTensor1 xpts(mpoints);
            dTensor2 qvals(mpoints,meqn);
            dTensor2 fvals(mpoints,mlength);

            // Loop over each quadrature point
            for (int m=1; m<=(mpoints-mopt); m++)
            {
                // grid point x
                xpts.set( m, xc + 0.5*dx*spts.get(m) );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {
                    // Sample q at xpts(:)
                    qvals.set(m,me, 0.0 );
                    for (int k=1; k<=mpoints; k++)
                    {
                        qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin.get(i,me,k) );
                    }

                    // TODO - HARD CODED ADVECTION EQUATION HERE!
                    fvals.set(m, me, qvals.get(m,me)*speeds.get(me) );
                    //fvals.set(m, me, qvals.get(m,me) );

                }

                // Call user-supplied function to set fvals
                // Func(xpts,qvals,auxvals,fvals);

            }

            // Evaluate integrals
            if (mopt==0) // project onto Legendre basis
            {
                for (int m1=1; m1<=mlength; m1++)
                    for (int m2=1; m2<=mpoints; m2++)
                    {
                        double tmp = 0.0;
                        for (int k=1; k<=mpoints; k++)
                        {
                            tmp = tmp + wgt.get(k)*fvals.get(k,m1)
                                *phi.get(k,m2);
                        }
                        Fout.set(i,m1,m2, tmp/s_area );         
                    }

            }
            else // project onto derivatives of Legendre basis
            {
                for (int m1=1; m1<=mlength; m1++)             
                for (int m2=1; m2<=mpoints; m2++)
                {
                    double tmp = 0.0;
                    for (int k=1; k<=(mpoints-mopt); k++)
                    {
                        tmp = tmp + wgt.get(k)*fvals.get(k,m1)
                            *phi_x.get(k,m2);
                    }
                    Fout.set(i,m1,m2, tmp/s_area );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
//
// Ripped, and gutted from the 1D library ... -DS
//
// All-purpose routine for computing the L2-projection
// of various functions onto:
//     mopt==0:   the Legendre basis
//     mopt==1:   the derivatives of Legendre basis
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC3 auxin(1-mbc:mnodes+mbc,maux,mpoints)
//           dTensorBC3   qin(1-mbc:mnodes+mbc,meqn,mpoints)
//           dTensorBC3  Fout(1-mbc:mnodes+mbc,mlength,mpoints)
// ---------------------------------------------------------------------

void L2ProjectSource(double t, int istart, int iend,
        double dx, 
        const dTensor1& speeds,  
        dTensorBC3& Fout,
        void (*Func)(
            double t, 
            const dTensor1& speeds, 
            const dTensor1& xpts,
            dTensor2& fvals))
{    

    const int KMAX = 6;
    const int meqn = Fout.getsize(2);
    const int mlength = Fout.getsize(2);
    const int mpoints = Fout.getsize(3);  
    dTensor1 wgt(mpoints), spts(mpoints);
    dTensor2 phi(mpoints,KMAX), phi_x(mpoints,KMAX);

    // ---------------------------------
    // Set quadrature weights and points
    // ---------------------------------
    switch ( mpoints )
    {
        case 1:
            wgt.set(1,  2.0e0 );
            spts.set(1, 0.0e0 );

            break;

        case 2:
            wgt.set(1,   1.0 );
            wgt.set(2,   1.0 );

            spts.set(1, -1.0/sq3 );
            spts.set(2,  1.0/sq3 );

            break;

        case 3:
            wgt.set(1, 5.0e0/9.0e0 );
            wgt.set(2, 8.0e0/9.0e0 );
            wgt.set(3, 5.0e0/9.0e0 );

            spts.set(1, -sq3/sq5 );
            spts.set(2,  0.0e0 );
            spts.set(3,  sq3/sq5 );

            break;

        case 4:
            wgt.set(1, (18.0 - sq30)/36.0 );
            wgt.set(2, (18.0 + sq30)/36.0 );
            wgt.set(3, (18.0 + sq30)/36.0 );
            wgt.set(4, (18.0 - sq30)/36.0 );

            spts.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            spts.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            spts.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            spts.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

            break;

        case 5:
            wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
            wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
            wgt.set(3,  128.0/225.0 );
            wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
            wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

            spts.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            spts.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            spts.set(3,  0.0 );
            spts.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            spts.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

            break;

        case 6:

            wgt.set(1,  0.171324492379170);
            wgt.set(2,  0.360761573048139);
            wgt.set(3,  0.467913934572691);
            wgt.set(4,  0.467913934572691);
            wgt.set(5,  0.360761573048139);
            wgt.set(6,  0.171324492379170);


            spts.set(1, -0.932469514203152);
            spts.set(2, -0.661209386466265);
            spts.set(3, -0.238619186083197);
            spts.set(4, -spts.get(3) );
            spts.set(5, -spts.get(2) );
            spts.set(6, -spts.get(1) );

            break;

    }

    // Loop over each quadrature point to construct Legendre polys
    for (int m=1; m<=mpoints; m++)
    {
        const double xi  = spts.get(m);
        const double xi2 = xi*xi;
        const double xi3 = xi2*xi;
        const double xi4 = xi3*xi;
        const double xi5 = xi4*xi;

        // Legendre basis functions at each grid point
        phi.set( m,1, 1.0 );
        phi.set( m,2, sq3*xi );
        phi.set( m,3, 0.5*sq5*( 3.0*(xi2) - 1.0 ) );
        phi.set( m,4, 0.5*sq7*xi
                *(5.0*(xi2) - 3.0) );
        phi.set( m,5, (105.0/8.0)*(xi4) 
                - (45.0/4.0)*(xi2) + (9.0/8.0) );
        phi.set( m, 6, (63.0/8.0)*sq11 * ( (xi5) - (10.0/9.0)*(xi3)+ (5.0/21.0)*xi ) );

        // 1st derivative of Legendre basis functions at each grid point
//      phi_x.set( m, 1, 0.0 );
//      phi_x.set( m, 2, 2.0*sq3/dx );
//      phi_x.set( m, 3, 6.0*sq5*xi/dx );
//      phi_x.set( m, 4, 3.0*sq7*(5.0*(xi2)-1.0)/dx );
//      phi_x.set( m, 5, 15.0*xi* (7.0*(xi2)-3.0)/dx );
//      phi_x.set( m, 6, (2.0/dx)*(63.0/8.0)*sq11*(5.0*(xi4)-(10.0/3.0)*(xi2)+5.0/21.0) );

    }

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------    
    const double s_area = 2.0;
    for (int i=istart; i<=iend; i++)
    {

        double xc = dogParamsCart2.get_xlow() + (double(i)-0.5)*dx;

        // each of these three items needs to be private to each thread ..
        dTensor1 xpts(mpoints);
        dTensor2 qvals(mpoints,meqn);
        dTensor2 fvals(mpoints,mlength);

        // Loop over each quadrature point
        for (int m=1; m<=mpoints; m++)
        {
            // grid point x
            xpts.set( m, xc + 0.5*dx*spts.get(m) );


            // Call user-supplied function to set fvals
            Func( t, speeds, xpts, fvals);

        }

        // Evaluate integrals
        // project onto Legendre basis
        for (int m1=1; m1<=mlength; m1++)
        for (int m2=1; m2<=mpoints; m2++)
        {
            double tmp = 0.0;
            for (int k=1; k<=mpoints; k++)
            {
                tmp = tmp + wgt.get(k)*fvals.get(k,m1) *phi.get(k,m2);
            }
            Fout.set(i,m1,m2, tmp/s_area );         
        }

    }
}
