#include <cmath>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data.h"
#include "DogSolveSL.h"   // function declarations 
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "Legendre2d.h"

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

void Set_dt_stages(double dt, int morder, dTensor1* dt_stages );
void SetSplitTime (double dt, int morder, double tn,  SL_state& sl_state,
        dTensor1* dt_stages );

void StrangStep( double dt, double tn, dTensor1* dt_stages, const dTensor2& phi, 
   dTensorBC2& u1, dTensorBC2& u2, SL_state& sl_state,
   dTensorBC3& smax, 
   dTensorBC4& qold,   // non-const because of the pos. limiter
   dTensorBC4& aux, dTensorBC4& qstar, dTensorBC4& qnew,
   DogSolverCart2& solver );

// semi-lagrangian solver
void DogSolverCart2::DogSolveUser(double tstart, double tend)
{

    // this accomodates the maximum number of stages allowed for the
    // split method ... and is a work in progress ...
    const int MAX_STAGES = 13;

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

    dTensorBC4   qstar (mx,my,meqn,kmax,mbc);    //temporary place holder
    dTensorBC4   qstar2(mx,my,meqn,kmax,mbc);    //temporary place holder
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

    // sample grid points (gauss points)
    dTensor2* spts = new dTensor2(mpoints, 2);

    //legendre polys evaluated at spts
    dTensor2 phi(mpoints, kmax);
    dTensor1 x1d(mpoints1d);
    dTensor1 wgt(mpoints1d);

    // ---------------------------------
    // Set quadrature points
    // ---------------------------------
    switch ( mpoints1d )
    {
        case 1:

            x1d.set(1, 0.0 );
            wgt.set(1,  2.0e0 );
            break;

        case 2:

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            wgt.set(1,   1.0 );
            wgt.set(2,   1.0 );
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

    // Tensor product Gaussian Quadrature
    //  See note at top of code in how mpoints are arranged here...
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
        for (int m2=1; m2<=(mpoints1d); m2++)
        {

            k = k+1;

            //save gauss quad grid point location on interval [-1,1]^2
            spts->set(k,2, x1d.get(m1) );
            spts->set(k,1, x1d.get(m2) );
        }

    //evaluate the legendre polynomials at sample points
    for(int m=1; m <= mpoints; m++)
    {

        double xi, xi2,xi3,xi4, eta, eta2,eta3,eta4;

        //find xi and eta (point to be evaluated)
        xi = spts->get(m,1);
        eta = spts->get(m,2);
        xi2 = xi*xi;
        xi3 = xi*xi2;
        xi4 = xi*xi3;
        eta2 = eta*eta;
        eta3 = eta*eta2;
        eta4 = eta*eta3;

        // Legendre basis functions evaluated at (xi,eta) in the
        // interval [-1,1]x[-1,1].
        switch( mpoints1d )
        {
            case 5:  // fifth order
                phi.set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi.set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi.set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 4:  // fourth order
                phi.set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 3:  // third order
                phi.set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set( m,5,  sq5*(1.5*xi2 - 0.5) );
                phi.set( m,4,  3.0*xi*eta );            

            case 2:  // second order            
                phi.set( m,3, sq3*eta );
                phi.set( m,2, sq3*xi  );

            case 1:  // first order
                phi.set( m,1, 1.0 );
        }

    }//end of evaluating legendre polys at sample grid points indexed by spts
    delete spts;
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

    const int numel = qnew.numel();

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
            cout << " Error in DogSolveAdvec.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
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

                case 0:  // used for testing - dont' take any time steps!
                    BeforeStep (dt, aux, qnew, *this);
                    AfterStep  (dt, aux, qnew, *this);
                    perror( "    case 0: Not taking a time step! " );
                    break;

                case 1: // 1st order in time, no corrections

                    sl_state.t     = tn;

                    BeforeStep(dt, aux, qstar, *this);

                    SetAdvecSpeed(phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    SetAdvecSpeed(phi, qstar, aux, smax, 2, u1, u2, sl_state);

                    StepAdvec(dt, qold, qstar, aux, u1, u2, 1, sl_state); 
                    StepAdvec(dt, qstar, qnew, aux, u1, u2, 2, sl_state);

                    break;

                case 2:  // 2nd order in time (strang split method)

                    StrangStep( dt, tn, dt_stages, phi, u1, u2, sl_state, smax, 
                        qold, aux, qstar, qnew, *this );

                    break;

                case 3: // Third-order splitting (with positive coefficients)
                    
                    // Initial guess (apply method with time step h/2)
                    StrangStep( 0.5*dt, tn, dt_stages, phi, u1, u2, sl_state, smax, 
                        qold, aux, qstar, qstar2, *this );
                    StrangStep( 0.5*dt, tn+0.5*dt, dt_stages, phi, u1, u2, sl_state, smax, 
                        qstar2, aux, qstar, qnew, *this );

                    // Coarse step (apply method with time step h)
                    StrangStep(   dt, tn, dt_stages, phi, u1, u2, sl_state, smax, 
                        qold, aux, qstar, qstar2, *this );

#pragma omp parallel for                    
                        for( int k=0; k < qnew.numel(); k++ )
                        {
                            qnew.vset( k, ( (4.)*qnew.vget(k) - qstar2.vget(k) )/3.);
                        }

                    break;

                case 4: // 4th order method (Yoshida Splitting)

                    // initial setup ... Save all appropriate times into SL_state
                    SetSplitTime(dt, 4, tn, sl_state, dt_stages );

                    ///////////////////////////////////////////////////////////
                    ///// There are 7 stages for Yoshida Splitting        /////
                    ///////////////////////////////////////////////////////////


                    ///////// stage 1: A //////////////////////////////////////
                    sl_state.stage_num = 1;
                    sl_state.t         = sl_state.split_time->get(1,1);
                    SetAdvecSpeed( phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(1), qstar, qnew,  aux, u1, u2,  1, sl_state ); 

                    ///////// stage 2: B //////////////////////////////////////
                    sl_state.stage_num = 2;
                    sl_state.t         = sl_state.split_time->get(2,2);
                    SetAdvecSpeed( phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(2), qnew, qstar,  aux, u1, u2, 2, sl_state); 

                    ///////// stage 3: A //////////////////////////////////////
                    sl_state.stage_num = 3;
                    sl_state.t         = sl_state.split_time->get(3,1);
                    SetAdvecSpeed( phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(3), qstar , qnew, aux, u1, u2, 1, sl_state);

                    ///////// stage 4: B //////////////////////////////////////
                    sl_state.stage_num = 4;
                    sl_state.t         = sl_state.split_time->get(4,2);
                    SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(4), qnew, qstar,  aux, u1, u2, 2, sl_state); 

                    ///////// stage 5: A //////////////////////////////////////
                    sl_state.stage_num = 5;
                    sl_state.t         = sl_state.split_time->get(5,1);
                    SetAdvecSpeed( phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(5), qstar, qnew, aux, u1, u2, 1, sl_state);

                    ///////// stage 6: B //////////////////////////////////////
                    sl_state.stage_num = 6;
                    sl_state.t         = sl_state.split_time->get(6,2);
                    SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(6), qnew, qstar, aux, u1, u2, 2, sl_state); 

                    ///////// stage 7: A //////////////////////////////////////
                    sl_state.stage_num = 7;
                    sl_state.t         = sl_state.split_time->get(7,1);
                    SetAdvecSpeed( phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    StepAdvec(sl_state.dt_stages->get(7), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //////////////// --- Experimental New Time Stepping --- ///////////////////

                    //                  ///////// stage 1: A //////////////////////////////////////
                    //                  sl_state.stage_num = 1;
                    //                  sl_state.t         = sl_state.split_time->get(1,1);
                    //                  SetAdvecSpeed( phi, qstar, aux, smax, 1, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(1), qstar, qnew,  aux, u1, u2,  1, sl_state ); 

                    //                  ///////// stage 2: B //////////////////////////////////////
                    //                  sl_state.stage_num = 2;
                    //                  sl_state.t         = sl_state.split_time->get(2,2);
                    //                  SetAdvecSpeed( phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(2), qnew, qstar,  aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 3: A //////////////////////////////////////
                    //                  sl_state.stage_num = 3;
                    //                  sl_state.t         = sl_state.split_time->get(3,1);
                    //                  StepAdvec(sl_state.dt_stages->get(3), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //                  ///////// stage 4: B //////////////////////////////////////
                    //                  sl_state.stage_num = 4;
                    //                  sl_state.t         = sl_state.split_time->get(4,2);
                    //                  SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(4), qnew, qstar,  aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 5: A //////////////////////////////////////
                    //                  sl_state.stage_num = 5;
                    //                  sl_state.t         = sl_state.split_time->get(5,1);
                    //                  StepAdvec(sl_state.dt_stages->get(5), qstar, qnew, aux, u1, u2, 1, sl_state);

                    //                  ///////// stage 6: B //////////////////////////////////////
                    //                  sl_state.stage_num = 6;
                    //                  sl_state.t         = sl_state.split_time->get(6,2);
                    //                  SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(6), qnew, qstar, aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 7: A //////////////////////////////////////
                    //                  sl_state.stage_num = 7;
                    //                  sl_state.t         = sl_state.split_time->get(7,1);
                    //                  StepAdvec(sl_state.dt_stages->get(7), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //                  ///////// stage 8: B //////////////////////////////////////
                    //                  sl_state.stage_num = 8;
                    //                  sl_state.t         = sl_state.split_time->get(8,2);
                    //                  SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(8), qnew, qstar, aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 9: A //////////////////////////////////////
                    //                  sl_state.stage_num = 9;
                    //                  sl_state.t         = sl_state.split_time->get(9,1);
                    //                  StepAdvec(sl_state.dt_stages->get(9), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //                  ///////// stage 10: B //////////////////////////////////////
                    //                  sl_state.stage_num = 10;
                    //                  sl_state.t         = sl_state.split_time->get(10,2);
                    //                  SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(10), qnew, qstar, aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 11: A //////////////////////////////////////
                    //                  sl_state.stage_num = 11;
                    //                  sl_state.t         = sl_state.split_time->get(11,1);
                    //                  StepAdvec(sl_state.dt_stages->get(11), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //                  ///////// stage 12: B //////////////////////////////////////
                    //                  sl_state.stage_num = 12;
                    //                  sl_state.t         = sl_state.split_time->get(12,2);
                    //                  SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
                    //                  StepAdvec(sl_state.dt_stages->get(12), qnew, qstar, aux, u1, u2, 2, sl_state); 

                    //                  ///////// stage 13: A //////////////////////////////////////
                    //                  sl_state.stage_num = 13;
                    //                  sl_state.t         = sl_state.split_time->get(13,1);
                    //                  StepAdvec(sl_state.dt_stages->get(13), qstar , qnew, aux, u1, u2, 1, sl_state);

                    //////////////// --- Experimental New Time Stepping --- ///////////////////

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
                cout << setprecision(3);
                cout << "DogSolve2D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t << endl;
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
                    cout<<"DogSolve2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
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
        cout << "         DogSolveAdvec Running Time = " 
            << setw(3) << time2 - time1 << "  seconds" << endl << endl;        
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

void SetSplitTime(double dt, int morder, double tn, SL_state& sl_state,
        dTensor1* dt_stages )
{

    sl_state.tn = sl_state.t = tn;
    if( morder == 2 )
    {
        dt_stages->set(1, 0.5*dt );
        dt_stages->set(2,     dt );
        dt_stages->set(3, 0.5*dt );

        sl_state.split_time->set(1,1,tn);
        sl_state.split_time->set(1,2,tn);
        sl_state.split_time->set(2,1,tn+0.5*dt);
        sl_state.split_time->set(2,2,tn);
        sl_state.split_time->set(3,1,tn+0.5*dt);
        sl_state.split_time->set(3,2,tn+1.0*dt);
        return;
    }
    else if( morder == 4)
    {

        double gamma1 =  1.351207191959658;
        double gamma2 = -1.702414383919315;

        dt_stages->set(1, 0.5*(gamma1)*dt);
        dt_stages->set(2,      gamma1*dt);
        dt_stages->set(3, 0.5*(gamma1+gamma2)*dt);
        dt_stages->set(4,      gamma2*dt);
        dt_stages->set(5, 0.5*(gamma1+gamma2)*dt);
        dt_stages->set(6,      gamma1*dt);
        dt_stages->set(7, 0.5*(gamma1)*dt);

        double t1, t2;
        // stage 1
        sl_state.tn = t1 = t2 = tn;
        sl_state.split_time->set(1,1,t1);
        sl_state.split_time->set(1,2,t2);
        t1 += dt_stages->get(1);

        // stage 2
        sl_state.split_time->set(2,1,t1);
        sl_state.split_time->set(2,2,t2);
        t2 += dt_stages->get(2);

        // stage 3
        sl_state.split_time->set(3,1,t1);
        sl_state.split_time->set(3,2,t2);
        t1 += dt_stages->get(3);

        // stage 4
        sl_state.split_time->set(4,1,t1);
        sl_state.split_time->set(4,2,t2);
        t2 += dt_stages->get(4);

        // stage 5
        sl_state.split_time->set(5,1,t1);
        sl_state.split_time->set(5,2,t2);
        t1 += dt_stages->get(5);

        // stage 6
        sl_state.split_time->set(6,1,t1);
        sl_state.split_time->set(6,2,t2);
        t2 += dt_stages->get(6);

        // stage 7
        sl_state.split_time->set(7,1,t1);
        sl_state.split_time->set(7,2,t2);
        t1 += dt_stages->get(7);

        // --- New (13-stage) Split Method --- //
        //
        // coefficients selected from:
        //
        // http://www.sciencedirect.com/science/article/pii/S0377042701004927
        //
        // "Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström
        // methods"

        //      dTensor1 b_even( 7 );
        //      dTensor1 b_odd ( 6 );

        //      b_even.set(1,  0.0829844064174052 );
        //      b_even.set(2,   0.396309801498368 );
        //      b_even.set(3, -0.0390563049223486 );
        //      b_even.set(4,  1. - 2.*(b_even.get(1)+b_even.get(2)+b_even.get(3) ) );
        //      b_even.set(5,  b_even.get(3) );
        //      b_even.set(6,  b_even.get(2) );
        //      b_even.set(7,  b_even.get(1) );

        //      b_odd.set(1,  0.245298957184271 );
        //      b_odd.set(2,  0.604872665711080 );
        //      b_odd.set(3,  0.5-(b_odd.get(1)+b_odd.get(2) ) );
        //      b_odd.set(4,  b_odd.get(3) );
        //      b_odd.set(5,  b_odd.get(2) );
        //      b_odd.set(6,  b_odd.get(1) );

        //      dt_stages->set(1, b_even.get(1) * dt);
        //      dt_stages->set(2, b_odd.get(1)  * dt);
        //      dt_stages->set(3, b_even.get(2) * dt);
        //      dt_stages->set(4, b_odd.get(2)  * dt);
        //      dt_stages->set(5, b_even.get(3) * dt);
        //      dt_stages->set(6, b_odd.get(3)  * dt);
        //      dt_stages->set(7, b_even.get(4) * dt);
        //      dt_stages->set(8, b_odd.get(4)  * dt);
        //      dt_stages->set(9, b_even.get(5) * dt);
        //      dt_stages->set(10, b_odd.get(5)  * dt);
        //      dt_stages->set(11, b_even.get(6) * dt);
        //      dt_stages->set(12, b_odd.get(6)  * dt);
        //      dt_stages->set(13, b_even.get(7) * dt);

        //      // stage 1
        //      sl_state.tn = t1 = t2 = tn;
        //      sl_state.split_time->set(1,1,t1);
        //      sl_state.split_time->set(1,2,t2);
        //      t1 += dt_stages->get(1);

        //      // stage 2
        //      sl_state.split_time->set(2,1,t1);
        //      sl_state.split_time->set(2,2,t2);
        //      t2 += dt_stages->get(2);

        //      // stage 3
        //      sl_state.split_time->set(3,1,t1);
        //      sl_state.split_time->set(3,2,t2);
        //      t1 += dt_stages->get(3);

        //      // stage 4
        //      sl_state.split_time->set(4,1,t1);
        //      sl_state.split_time->set(4,2,t2);
        //      t2 += dt_stages->get(4);

        //      // stage 5
        //      sl_state.split_time->set(5,1,t1);
        //      sl_state.split_time->set(5,2,t2);
        //      t1 += dt_stages->get(5);

        //      // stage 6
        //      sl_state.split_time->set(6,1,t1);
        //      sl_state.split_time->set(6,2,t2);
        //      t2 += dt_stages->get(6);

        //      // stage 7
        //      sl_state.split_time->set(7,1,t1);
        //      sl_state.split_time->set(7,2,t2);
        //      t1 += dt_stages->get(7);

        //      sl_state.split_time->set(8,1,t1);
        //      sl_state.split_time->set(8,2,t2);
        //      t2 += dt_stages->get(8);

        //      sl_state.split_time->set(9,1,t1);
        //      sl_state.split_time->set(9,2,t2);
        //      t1 += dt_stages->get(9);

        //      sl_state.split_time->set(10,1,t1);
        //      sl_state.split_time->set(10,2,t2);
        //      t2 += dt_stages->get(10);

        //      sl_state.split_time->set(11,1,t1);
        //      sl_state.split_time->set(11,2,t2);
        //      t1 += dt_stages->get(11);

        //      sl_state.split_time->set(12,1,t1);
        //      sl_state.split_time->set(12,2,t2);
        //      t2 += dt_stages->get(12);

        //      sl_state.split_time->set(13,1,t1);
        //      sl_state.split_time->set(13,2,t2);
        //      t1 += dt_stages->get(13);

    }
}


void StrangStep( double dt, double tn, dTensor1* dt_stages, const dTensor2& phi, 
   dTensorBC2& u1, dTensorBC2& u2, SL_state& sl_state,
   dTensorBC3& smax, 
   dTensorBC4& qold, dTensorBC4& aux, dTensorBC4& qstar, dTensorBC4& qnew,
   DogSolverCart2& solver )
{

    SetSplitTime(dt, 2, tn, sl_state, dt_stages );
    sl_state.stage_num = 1;
    sl_state.t = tn;
    SetAdvecSpeed( phi, qold, aux, smax, 1, u1, u2, sl_state);  
    StepAdvec(0.5*dt, qold, qnew, aux, u1, u2, 1, sl_state ); 

    // Poisson solve called in BeforeStep for VP system  
    sl_state.stage_num = 2;
    BeforeStep(dt, aux, qnew, solver);

    sl_state.t  = tn;
    SetAdvecSpeed(phi, qnew, aux, smax, 2, u1, u2, sl_state);
    StepAdvec(dt, qnew, qstar, aux, u1, u2, 2, sl_state ); 

    sl_state.stage_num = 3;
    sl_state.t = tn + 0.5*dt;
    StepAdvec(0.5*dt, qstar, qnew, aux, u1, u2, 1, sl_state);

}
