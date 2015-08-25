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

//////////////////////////////////////////////////////////////////////////////
// Solver written specifically for the constant coefficient advection equation.
//
//  q_t + u q_x + v q_y = 0
//
//////////////////////////////////////////////////////////////////////////////

// semi-lagrangian solver
void DogSolveAdvecCC(DogSolverCart2& solver, double tstart, double tend)
{

    const edge_data& EdgeData = solver.get_EdgeData();
    dTensorBC3& smax = solver.fetch_smax();
    dTensorBC4& qold = solver.fetch_qold();
    DogStateCart2* dogStateCart2 = &solver.fetch_dogStateCart2();

    dTensorBC4& qnew = dogStateCart2->fetch_q();
    dTensorBC4& aux  = dogStateCart2->fetch_aux();
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
    double dt = dogStateCart2->get_initial_dt(); /// ???? ////
    const double CFL_max = cflv[1];
    const double CFL_target = cflv[2];
    double cfl = 0.0;
    double dtmin = dt;
    double dtmax = dt;

    dTensorBC4   qstar(mx,my,meqn,kmax,mbc);    //temporary place holder
    dTensorBC4 aux_old(mx,my,maux,kmax,mbc);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // sample grid points (gauss points)

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

            CopyQ(qold, qstar);
            /////////////////////////////////////
            // Perform a full time step
            /////////////////////////////////////
            
            // case 1: // 1st order in time, no corrections
                    
            *t_ptr     = tn;
            *(t_ptr+1) = tn+dt;

            sl_state.t     = tn;
            sl_state.ts    = tn;

            BeforeStep(dt, aux, qstar, solver);

            SetAdvecSpeed(phi, qstar, aux, smax, 1, u1, u2, sl_state);
            SetAdvecSpeed(phi, qstar, aux, smax, 2, u1, u2, sl_state);

            StepAdvec(dt, qstar, qstar, aux, u1, u2, 1, sl_state); 
            StepAdvec(dt, qstar, qnew, aux, u1, u2, 2, sl_state);

            // compute cfl number
            double GetCFL(double dt, const dTensorBC4& aux, 
                    const dTensorBC3& smax);
            cfl = GetCFL(dt,aux,smax);

            // output time step information
//          if (dogParams.get_verbosity()>0) 
//          {
//              cout << setprecision(3);
//              cout << "DogSolve2D ... Step" << setw(5) << n_step;
//              cout << "   CFL =" << setw(6) << fixed << cfl;
//              cout << "   dt =" << setw(11) << scientific << dt;
//              cout << "   t =" << setw(11) << scientific << t << endl;
//              //cout << " CFL_target ="<<setw(6)<<fixed<<CFL_target<<endl;
//          }

            if (cfl>0.0)
            {   dt = Min(dogParams.get_max_dt(),dt*CFL_target/cfl); }
            else
            { dt = dogParams.get_max_dt(); }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
            {   // accept 
                m_accept = 1;
                dogStateCart2->set_time(t);
                sl_state.tnm1 = tn;
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

        void AfterFullSLTimeStep(
            dTensorBC4& aux, dTensorBC4& qnew, double t );
        AfterFullSLTimeStep(aux, qnew, t );

        // apply the limiter - This way global integrals stay positive
        void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
        if(dogParams.using_moment_limiter())
        { ApplyPosLimiter(aux, qnew); }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    // set initial time step for next call to DogSolveAdvec
    dogStateCart2->set_initial_dt(dt);

}
///////////////////////////////////////////////////////////////////////////////
