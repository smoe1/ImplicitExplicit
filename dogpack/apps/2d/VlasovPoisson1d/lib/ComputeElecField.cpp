#include <math.h>
#include "dogdefs.h"
#include "DogParams.h"       // for morder
#include "DogParamsCart2.h"
#include "VlasovParams.h"    // integration constant is stored here

///////////////////////////////////////////////////////////////////////////////
//
//  Function to compute the electric field for the Vlasov Poisson Eqn.
//
//  For a fixed time t0, setting
//
//          frhs(x) = -CONST + \int_v f(x,v,t0)\ dv
//
//  We have from Maxwell's equations:
//
//          phi_xx = -frhs(x)
//
//  And the electric field is given by
//
//          E(x,t_0) = -phi_x.
//
//  Parameters:
//
//     node(mx+1,my+1, mdim) - the space-velocity nodal points (mdim=2)
//
//     node1d(mx+1, mdim1d) - the space nodal points (mdim1d=1)
//
//     qvals(mx,my,meqn,kmax, mbc)  - 2d Legendre weights for prob function f
//
//  Returns:
//
//     Evals(melems1d, meqn, kmax1d, mbc) = 
//            1d legendre coefficient weights for electric field
//
///////////////////////////////////////////////////////////////////////////////

void ComputeElecField(double t, const dTensor2& node1d, const dTensorBC4& qvals, 
        dTensorBC3& Evals)
{
    void ComputeElecField(double t, 
            const dTensor2* node1d, 
            const dTensorBC4* qvals, 
            dTensorBC3* Evals);
    ComputeElecField(t,&node1d,&qvals,&Evals);
}

void ComputeElecField(double t, 
        const dTensor2* node1d, 
        const dTensorBC4* qvals, 
        dTensorBC3* Evals)
{

    // --------------------- Function Declarations -------------------------//
    void ConvertBC(const dTensor2* node, 
            const int mbc1, 
            const int mbc2, 
            const double tmp1,
            const double tmp2, 
            const dTensorBC3& Fvals, 
            double& gamma, 
            double& beta);

    void PoissonSolve(const int mstart,
            const dTensor2* node, 
            const double gamma,
            const double beta, 
            const dTensorBC3& Fvals,
            dTensorBC3* qvals);

    void IntegrateQ1d(const int mopt,
            const dTensorBC4* q2d, 
            dTensorBC3* q1d);
    // --------------------- Function Declarations -------------------------//

    // --------------------- Local Variables -------------------------------//
    const int mx   = qvals->getsize(1);
    const int my   = qvals->getsize(2);
    const int meqn = qvals->getsize(3);
    const int kmax = qvals->getsize(4);
    const int mbc  = qvals->getmbc();
    const int kmax1d  = dogParams.get_space_order();
    const int mpoints = kmax1d*kmax1d;

    dTensorBC3 fvals(mx, 2 , kmax1d, mbc, 1);  // right hand side function

    // phi(mx,1,kmax1d) = 'electric' field
    // phi(mx,2,kmax1d) = potential  (not unique!)
    dTensorBC3 phi(mx, 2, kmax1d, mbc, 1);  // Electric field and potential

    //////////////////////////////////
    // compute \int q(x,v)\ dv
    //////////////////////////////////
    IntegrateQ1d(1, qvals, &fvals);

    for(int i=1; i<= mx; i++)
    { 
        double tmp = fvals.get(i,1,1);
        fvals.set(i, 1, 1, tmp - vlasovParams.get_integration_constant() ); 
    }

    //////////////////////////////////////////////////////////////////////////
    // Solve for electric field and potential
    //////////////////////////////////////////////////////////////////////////

    // periodic boundary conditions: alpha = beta = 0;
    //       gamma = delta will be forced provided \int f = 0.
    double gamma, beta;
    ConvertBC(node1d, 1, 1, 0.0, 0.0, fvals, gamma, beta);
    PoissonSolve(1, node1d, gamma, beta, fvals, &phi);

    //////////////////////////////////////////////////////////////////////////
    // Save results into Evals
    //////////////////////////////////////////////////////////////////////////
    for(int i = 1; i <= mx; i++ )
        for(int k=1; k <= kmax1d; k++ )
        { Evals->set(i, 1, k, phi.get(i,  1, k) ); }

}// end of Compute Electric Field


void ComputeElecPotential(double t, 
        const dTensor2& node1d, 
        const dTensorBC4& qvals, 
        dTensorBC3& phi)
{
    void ComputeElecPotential(double t, 
            const dTensor2* node1d, 
            const dTensorBC4* qvals, 
            dTensorBC3* phi);
    ComputeElecPotential(t,&node1d,&qvals,&phi);
}

void ComputeElecPotential(double t, 
        const dTensor2* node1d, 
        const dTensorBC4* qvals, 
        dTensorBC3* phi)
{
    /////////////////////////////////////////////////////////////////////////
    // This function is only ever called in Output.cpp -- the purpose is to
    // print the moments for both the electric potential, as well as the 
    // electric field
    /////////////////////////////////////////////////////////////////////////

    // --------------------- Function Declarations -------------------------//
    void ConvertBC(const dTensor2* node, 
            const int mbc1, 
            const int mbc2, 
            const double tmp1, 
            const double tmp2, 
            const dTensorBC3& Fvals, 
            double& gamma,
            double& beta);

    void PoissonSolve(const int mstart, 
            const dTensor2* node, 
            const double gamma, 
            const double beta, 
            const dTensorBC3& Fvals, 
            dTensorBC3* qvals);

    void IntegrateQ1d(const int mopt,
            const dTensorBC4* q2d, 
            dTensorBC3* q1d);
    // --------------------- Function Declarations -------------------------//

    // --------------------- Local Variables -------------------------------//
    const int mx   = qvals->getsize(1);
    const int my   = qvals->getsize(2);
    const int meqn = qvals->getsize(3);
    const int kmax = qvals->getsize(4);
    const int mbc  = qvals->getmbc();

    const int kmax1d  = dogParams.get_space_order();
    const int mpoints = kmax1d*kmax1d;

    dTensorBC3 fvals(mx, 2 , kmax1d, mbc, 1);  // right hand side function

    /////////////////////////////////////////////
    // compute \int q(x,v)\ dv - rho_0
    /////////////////////////////////////////////
    IntegrateQ1d(1, qvals, &fvals);
    for(int i=1; i<= mx; i++)
    { 
        double tmp = fvals.get(i,1,1);
        fvals.set(i, 1, 1, tmp - vlasovParams.get_integration_constant() ); 
    }

    //////////////////////////////////////////////////////////////////////////
    // Solve for electric field and potential
    //////////////////////////////////////////////////////////////////////////

    // periodic boundary conditions: alpha = beta = 0;
    //       gamma = delta will be forced provided \int f = 0.
    double gamma, beta;
    ConvertBC(node1d, 1, 1, 0.0, 0.0, fvals, gamma, beta);
    PoissonSolve(1, node1d, gamma, beta, fvals, phi);

}

// This tester function simply checks if the right hand side is a valid right
// hand side ...
void CheckSolvability(const dTensorBC3& fvals)
{
    const double dx = dogParamsCart2.get_dx();
    double tmp = 0.0;
    for(int i = 1; i <= fvals.getsize(1); i++ )
    {
        tmp += dx * fvals.get(i,1,1);
    }
    if( fabs(tmp) > 1e-10 )
    {
        printf("   Bad right hand side function");
        printf("   integral = %g\n", tmp );
        exit(1);
    }
}


