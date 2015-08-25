#include <cmath>
#include <sys/time.h>
#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data.h"
#include "DogSolveSL.h"   // function declarations 
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "Legendre2d.h"
#include "VPDataCart2.h"

//////////////////////////////////////////////////////////////////////////////
// solver written specifically for the advection equation.
// currently supports advection equations of the form:
//
//  q_t + u(y) q_x + v(x) q_y = 0
//
//////////////////////////////////////////////////////////////////////////////

// better time function:
// simple wrapper function for printing time spent running the routines
double GetTime()
{
    struct timeval t;
    int rc = gettimeofday(&t, NULL);
    assert(rc == 0);
    return t.tv_sec + t.tv_usec/1e6; //
}

// semi-lagrangian solver
void DogSolverCart2::DogSolveUser(double tstart, double tend)
{
    // this accomodates the maximum number of stages allowed for the
    // split method ... and is a work in progress ...
    const int MAX_STAGES = 13;

    void ComputeElecField(double t, 
            const dTensor2* node1d, 
            const dTensorBC4* qvals, 
            dTensorBC3* Evals);

    void ApplyPosLimiter(const dTensorBC4& aux, 
            dTensorBC4& q);

    const edge_data& EdgeData = Legendre2d::instance().get_edgeData();
    dTensorBC3& smax = fetch_smax();
    DogStateCart2* dogStateCart2 = &fetch_state();

    dTensorBC4& qnew = fetch_state().fetch_q();
    dTensorBC4& qold = fetch_state_old().fetch_q();
    dTensorBC4& aux  = fetch_state().fetch_aux();

    const int nv = dogParams.get_nv();
    const double* cflv = dogParams.get_cflv();
    // --------------------------------------------------------------
    // define local variables
    int m_accept;
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int maux = aux.getsize(3);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();

    const int ndims = 2;

    int n_step = 0;
    double t = tstart;
    double tn = t;
    double told = 0.0;
    double dt = get_dt();
    const double CFL_max = cflv[1];
    const double CFL_target = cflv[2];
    double cfl = 0.0;
    double dtmin = dt;
    double dtmax = dt;

    dTensorBC4   qstar(mx,my,meqn,kmax,mbc);    //temporary place holder
    dTensorBC4 aux_old(mx,my,maux,kmax,mbc);

    //////////////////////////////////////////////////////////////////////////
    // Setup for manipulating velocities.  Store 2 arrays, 1 for each velocity
    //////////////////////////////////////////////////////////////////////////
    // array for storing advection speeds in each cell
    const int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int mpoints1d = int(sqrt(mpoints));
    const int morder = mpoints1d;
    dTensorBC2 u1(my, mpoints1d,  mbc, ndims-1);    // speed, u(y) (ndims == 1 )
    dTensorBC2 u2(mx, mpoints1d,  mbc, ndims-1);    // speed, v(x) (ndims == 1 )

    ////////Set up any extra state variables associated with this problem /////
    SL_state sl_state;
    sl_state.split_time = new dTensor2(MAX_STAGES, 2);
    sl_state.aux1d      = new dTensorBC4( Max(mx,my), 2, 4, mpoints1d, mbc, 1);
    sl_state.node1d     = new dTensor2(mx+1,1);
    for( int i=1; i <=(mx+1); i++)
    { sl_state.node1d->set(i, 1, dogParamsCart2.get_xl(i) ); }

    sl_state.qold = &qold;
    sl_state.qnew = &qnew;

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    //legendre polys evaluated at quad points
    dTensor2 phi1d(mpoints1d,dogParams.get_space_order());
    dTensor1   x1d(mpoints1d);
    dTensor1   wgt(mpoints1d);

    // ---------------------------------
    // Set quadrature points
    // ---------------------------------
    switch ( mpoints1d )
    {
        case 1:

            x1d.set(1, 0.0 );
            wgt.set(1, 2.0 );
            break;

        case 2:

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            wgt.set(1,  1.0 );
            wgt.set(2,  1.0 );
            break;

        case 3:
            x1d.set(1,  -sq3/sq5 );
            x1d.set(2,  0.0 );
            x1d.set(3, sq3/sq5 );
            wgt.set(1, 5.0e0/9.0e0 );
            wgt.set(2, 8.0e0/9.0e0 );
            wgt.set(3, 5.0e0/9.0e0 );
            break;

        case 4:
            x1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d.set(3, -x1d.get(2) );
            x1d.set(4, -x1d.get(1) );          

            wgt.set(1, (18.0 - sqrt(30.0))/36.0 );
            wgt.set(2, (18.0 + sqrt(30.0))/36.0 );
            wgt.set(3, (18.0 + sqrt(30.0))/36.0 );
            wgt.set(4, (18.0 - sqrt(30.0))/36.0 );
            break;

        case 5:         

            x1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );

            wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
            wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
            wgt.set(3,  128.0/225.0 );
            wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
            wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );
            break;
    }

    //evaluate the 1d legendre polynomials at 1d sample points
    for(int m=1; m<=mpoints1d; m++)
    {
        const double xi  = x1d.get(m);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;

        // Legendre basis functions evaluated at (xi,eta) in the
        // interval [-1,1]x[-1,1].
        switch( mpoints1d )
        {
            case 5:  // fifth order
                phi1d.set( m,5, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );

            case 4:  // fourth order
                phi1d.set( m,4, sq7*(2.5*xi3 - 1.5*xi) );

            case 3:  // third order
                phi1d.set( m,3, sq5*(1.5*xi2 - 0.5) );

            case 2:  // second order            
                phi1d.set( m,2, sq3*xi  );

            case 1:  // first order
                phi1d.set( m,1, 1.0 );
        }

    }//end of evaluating legendre polys at sample grid points indexed by spts
    ///////////////////////////////////////////////////////////////////////////

    double time1, time2;   // running time values
    time1 = GetTime();     //running time for this program (can remove this..)

    // fourth order splitting coefficients
    dTensor1* dt_stages;
    if( dogParams.get_time_order() >= 2 )
    {
        dt_stages = new dTensor1(MAX_STAGES);
        sl_state.dt_stages = dt_stages;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Main Time Stepping Loop
    ///////////////////////////////////////////////////////////////////////////
    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            printf(" Error in DogSolveAdvec.cpp: "); 
            printf(" Exceeded allowed # of time steps \n");
            printf("    n_step = %i\n",n_step);
            printf("        nv = %i\n",nv);
            printf("\n");
            exit(1);
        }

        // copy qnew and aux in case we reject the step
        CopyQ(qnew,qold);
        CopyQ(qnew,qstar);
        CopyQ(aux,aux_old);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time          
            told = t;
            tn = t;
            if (told+dt > tend) { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero (this will be saved in
            // SetAdvecSpeed)
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)
                {
                    smax.set(i,j,1, 0.0e0 );
                    smax.set(i,j,2, 0.0e0 );
                }

            sl_state.dt = dt;
            sl_state.t  = sl_state.tn = tn;

            void InitSLState( const dTensorBC4& q, const dTensorBC4& aux, 
                    SL_state& sl_state );
            InitSLState( qnew, aux, sl_state );

            CopyQ(qold, qstar);
            /////////////////////////////////////
            // Perform a full time step
            /////////////////////////////////////
            switch ( dogParams.get_time_order() )
            {
                case 2:  // 2nd order in time (strang split method)

                    // Stage 1
                    SetAdvecSpeed_MOD(phi1d, qstar, vpDataCart2.Efield, 
                            vpDataCart2.v1d, smax, 2, u1, u2);  
                    StepAdvec(0.5*dt, qstar, qnew, aux, u1, u2, 2, sl_state ); 

                    // Stage 2              
                    SetAdvecSpeed_MOD(phi1d, qnew, vpDataCart2.Efield, 
                            vpDataCart2.v1d, smax, 1, u1, u2);
                    StepAdvec(dt, qnew, qstar, aux, u1, u2, 1, sl_state ); 

                    // Stage 3
                    for(int i=1; i<=mx; i++)
                        for(int k=1; k<=dogParams.get_space_order(); k++)
                        { vpDataCart2.Efield_old->set(i,1,k, vpDataCart2.Efield->get(i,1,k) ); }
                    ComputeElecField(tn, vpDataCart2.node1d,
                            &qstar, vpDataCart2.Efield);      
                    SetAdvecSpeed_MOD(phi1d, qstar, vpDataCart2.Efield, 
                            vpDataCart2.v1d, smax, 2, u1, u2);  
                    StepAdvec(0.5*dt, qstar, qnew, aux, u1, u2, 2, sl_state);

                    // More limiting to keep qnew positive at Gauss quadrature points          
                    if(dogParams.using_moment_limiter())
                    {  ApplyPosLimiter(aux, qnew);  }

                    // Correct electric field
                    void CorrectElectricField(const double dt,
                            const dTensorBC4* qold,
                            const dTensorBC4* qnew,
                            const dTensorBC3* Efield_old,
                            dTensorBC3* Efield);          
                    CorrectElectricField(dt,
                            &qold,
                            &qnew,
                            vpDataCart2.Efield_old,
                            vpDataCart2.Efield);
                    break; 

                default:
                    // still here?  too bad!  Pick a valid time stepping
                    // method
                    fprintf(stderr, 
                            "Bad time stepping method chosen in DogSolveAdvec\n");
                    fprintf(stderr, 
                            "dogParams.get_time_order() = %d\n",
                            dogParams.get_time_order());
                    exit(1);

            }// end of taking a full time step

            // compute cfl number
            cfl = GetCFL(dt);

            // output time step information
            if (dogParams.get_verbosity()>0) 
            {
                printf("DogSolve2D ... Step %5i",n_step);
                printf("   CFL = %6.3f",cfl);
                printf("   dt = %11.3e",dt);
                printf("   t = %11.3e\n",t);
            }

            if (cfl>0.0)
            {   dt = Min(dogParams.get_max_dt(),dt*CFL_target/cfl); }
            else
            { dt = dogParams.get_max_dt(); }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
            {   // accept 
                m_accept = 1;
                dogStateCart2->set_time(t);
                dogParams.set_time(t); // time hack
            }
            else 
            {   //reject
                t = told;
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step...");
                    printf("CFL number too large\n");

                    // find index of larger value...
                    int imax = 1; int jmax = 1; int dmax = 1;
                    for( int i = 1; i <= smax.getsize(1); i++ )
                        for( int j = 1; j <= smax.getsize(2); j++ )
                            for( int d = 1; d <= ndims; d++ )
                            {
                                if( smax.get(i,j,d) > smax.get(imax,jmax,dmax) )
                                { imax = i;   jmax = j; dmax = d;}
                            }
                    printf("   imax = %d, jmax = %d, dmax = %d\n", imax, jmax, dmax );
                    printf("   smax(imax,jmax,dmax) = %3.2e\n", smax.get(imax, jmax, dmax ) );

                }

                // copy qold into qnew
                CopyQ(qold,qnew);
                CopyQ(aux_old,aux);
            }
        }

        void AfterFullSLTimeStep( dTensorBC4& aux, dTensorBC4& qnew, double t );
        AfterFullSLTimeStep(aux, qnew, t );

        // apply the limiter - This way global integrals stay positive
        void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
        if(dogParams.using_moment_limiter())
        { ApplyPosLimiter(aux, qnew); }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    void AfterSLFrame(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2&
            solver, SL_state sl_state );
    AfterSLFrame(dt, aux, qnew, *this, sl_state );

    // set initial time step for next call to DogSolveAdvec
    set_dt(dt);

    // Clock the running time for this call to DogSolveAdvec
    time2 = GetTime();
    if(dogParams.get_verbosity()>0)
    {
        printf("         DogSolveAdvec Running Time = %3g seconds\n\n",time2-time1);
    }

    if( dogParams.get_time_order() >= 2 )
    {
        delete sl_state.split_time;
        delete sl_state.aux1d;
        delete sl_state.node1d;
        delete dt_stages;
    }


}
///////////////////////////////////////////////////////////////////////////////
