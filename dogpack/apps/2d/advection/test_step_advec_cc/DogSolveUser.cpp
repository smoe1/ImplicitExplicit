#include <stdlib.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"

// -- This routine is used in place of 'FluxFunc' -- //
void SetSpeeds( dTensor2* vel_vec )
{
    // forced three equations:

    // 1st equation:
    vel_vec->set(1,1, 1.0 );
    vel_vec->set(1,2, 0.0 );

    // 2nd equation:
    vel_vec->set(2,1, 0.0 );
    vel_vec->set(2,2, 2.0 );

    // 3rd equation:
    vel_vec->set(3,1,-2.0 );
    vel_vec->set(3,2,-1.0 );

}

// Example used to specifically test StepAdvecCC
//
void DogSolverCart2::DogSolveUser(double tstart, double tend)
{
    dTensorBC4& qnew = fetch_state().fetch_q();
    dTensorBC4& aux  = fetch_state().fetch_aux();
    dTensorBC4& qold = fetch_state_old().fetch_q();
    dTensorBC4& auxold = fetch_state_old().fetch_aux();
    dTensorBC3& smax = fetch_smax();
    const int nv = dogParams.get_nv();
    const double* cflv = dogParams.get_cflv();

    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // ------------------------------------------------------------
    // Function definitions
    void StepAdvecCC(double dt, dTensorBC4& qold, 
        const dTensor2& speeds, dTensorBC4& qnew );
    // ------------------------------------------------------------

    // -- Hackish way to set the speeds here! -- //
    dTensor2 speeds( meqn, 2 );
    void SetSpeeds( dTensor2* vel_vec );
    SetSpeeds( &speeds );
//  for( int me=1; me <= meqn; me++ )
//  {
//      printf("speeds(%d,1) = %f\n", me, speeds.get(me,1) );
//      printf("speeds(%d,2) = %f\n", me, speeds.get(me,2) );
//  }

    // define local variables
    int n_step = 0;
    double t = tstart;
    double dt = get_dt();
    const double CFL_max = cflv[1];
    const double CFL_target = cflv[2];
    
    dTensorBC4   qstar(mx, my, meqn, kmax, mbc);
    dTensorBC4 auxstar(mx, my, maux, kmax, mbc);

    // Set initialize qstar and auxstar values
    qstar.copyfrom( qold );
    auxstar.copyfrom( aux );

    // User-defined time stepping
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveUser.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // copy qnew into qold
        qold.copyfrom( qnew );
        auxold.copyfrom( aux );

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            set_time_hack(told);
            set_dt(dt);

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // Take a single step of the method
            StepAdvecCC( dt, qold, speeds, qnew );

            // compute cfl number
            //double cfl = GetCFL(dt);
            double cfl = 0.;
            for( int me = 1; me <= meqn; me++ )
            {
                cfl = Max( 
                    Max( cfl, fabs( speeds.get(me,1)*dt/dx ) ),
                              fabs( speeds.get(me,2)*dt/dy ) );
            }

            // output time step information
            if (dogParams.get_verbosity()>0) 
            {
                printf("DogSolve2D ... Step %5d"
                        "   CFL =%6.3f"
                        "   dt =%11.3e"
                        "   t =%11.3e\n",
                        n_step,cfl,dt,t);
            }

            // choose new time step
            if (cfl>0.0)
            {   dt = Min(dogParams.get_max_dt(), dt*CFL_target/cfl); }
            else
            { dt = dogParams.get_max_dt(); }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
            // accept
            { 
                m_accept = 1; 

                // do any extra work
                ::AfterFullTimeStep(fetch_solver());
            }
            else 
            //reject
            {   
                t = told;
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
                }

                // copy qold into qnew
                qnew.copyfrom( qold  );
                aux.copyfrom( auxold );
            }      
        }

        // apply the limiter - This way global integrals stay positive
        void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
        if(dogParams.using_moment_limiter())
        { ApplyPosLimiter(aux, qnew); }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    // set initial time step for next call to DogSolveUser
    set_dt(dt);

}
