#include "DogSolveSL.h"

// ------------------------------------------------------------
// parameters
//
// node: list of x position of node.  node(i,1) is the left endpoint of grid
// cell i.
//        1st index = which grid cell
//        2nd index == 1 (for 1d case)
//
// prim_vol: length of each grid cell
//
// aux (auxiliary arrays): helper arrays for updating q.  
//      this will be the velocity for the advection eqn
//        (../../advection/smooth_example/SetWaveSpd.cpp)
//
// qnew: q values
//
// smax: maximum wave speed.  need to make sure CFL condition is satisfied.
//
// tstart, tend: start and end time for updating
//
// nv: max number of time steps allowed
//
// method[]: read any dogpack.data file to see each option
//
// dtv[1] = initial time step
// dtv[2] = max allowed time step
//
// cflv[1] = max allowed CFL number
// cflv[2] = desired CFL number 
//
// outputdir = output directory.  default is ./output/
// ------------------------------------------------------------
// simple wrapper function for the desired time stepping method
void DogSolveUser(const dTensor2& node, const dTensor1& prim_vol,
        dTensorBC3& aux, dTensorBC3& qold,
        dTensorBC3& qnew, dTensorBC1& smax,
        double tstart, double tend,int nv, const int method[],
        double dtv[], const double cflv[],string outputdir)
{  
    DogSolveSL(node, prim_vol, aux, qold, qnew, smax, tstart, tend, nv, 
            method, dtv, cflv, outputdir);
}

void DogSolveSL(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
        dTensorBC3& qold, dTensorBC3& qnew, dTensorBC1& smax, double tstart, 
        double tend, int nv, const int method[], double dtv[], 
        const double cflv[], string outputdir)
{
    //double CFL_max;
    double t = tstart;
    double dt = dtv[1];            // initial time step used
    double CFL_max = cflv[1];
    double CFL_target = cflv[2];
    double cfl = 0.0;
    const double dx = node.get(2,1) - node.get(1,1);

    const int melems = qnew.getsize(1); //number of elements
    const int meqn   = qnew.getsize(2); // number of equations
    const int kmax   = qnew.getsize(3); 
    const int mbc    = qnew.getmbc();   //number of ghost cells
    const int maux   = aux.getsize(2);

    double *t_ptr = new double[2];
    dTensorBC3 aux_old(melems, maux, method[1], mbc);
    dTensorBC3 qstar  (melems, meqn, kmax, mbc);

    // Apply the Limiter - this way global integrals stay positive //
    void ApplyPosLimiter(const dTensorBC3& aux, dTensorBC3& q);
    if(dogParams.using_moment_limiter())
    { ApplyPosLimiter(aux, qnew); }

    ////////////////////////////////////////////////////////////////////////////
    // Main Time-Stepping loop  ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    int n_step = 0;
    while( t < tend )
    {    

        //initialize data
        int m_accept = 0;
        n_step++;

        //check if max number of time steps is exceeded
        if(n_step > nv)
        {
            cout << " Error in DogSolveAdvec.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }

        // copy qnew into qold (in order to save data in case we reject step)
        CopyQ(qnew,qold);
        CopyQ(aux,aux_old);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            double tn = *t_ptr = t;
            if (told+dt > tend) { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            for (int j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // take one full time step
            *(t_ptr+1) = tn;
            BeforeStep(dt, node, aux, qold, (void*)t_ptr );

            StepAdvec (dt, node, aux, smax, qold, qnew);
            //StepAdvecFluxForm(dt, node, aux, smax, qold, qnew);
            //StepAdvecNonCons (dt, node, aux, smax, qold, qnew);


            AfterStep (dt, node, aux, qnew);

            // Apply the Limiter - this way global integrals stay positive //
            if(dogParams.using_moment_limiter())
            { ApplyPosLimiter(aux, qnew); }

            // compute cfl number
            cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);

            // output time step information
            if (method[4]>0) 
            {
                cout << setprecision(3);
                cout << "DogSolve1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            double dt_old = dt;
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max) // accept  
            { m_accept = 1; }      
            else 
                //reject
            {   
                t = told;
                if (method[4]>0)
                {
                    cout<<"DogSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew, and aux_old into aux
                CopyQ(qold,    qnew);
                CopyQ(aux_old, aux);
            }

        }

        // compute conservation and print to file
        ConSoln(method,node,aux,qnew,t,outputdir);
    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

    delete t_ptr;

}
