#include <math.h>
#include "dogdefs.h"
#include "DogParams.h"       // for morder
#include "DogParamsCart1.h"
#include "QuadMomentParams.h"    // integration constant is stored here

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
//     node(mx+1, mdim) - the space nodal points (mdim=1)
//
//     qvals(mx,meqn,kmax, mbc)  - 1d Legendre weights for prob function f
//
//  Returns:
//
//     Evals(melems, meqn, kmax, mbc) = 
//            1d legendre coefficient weights for electric field
//
///////////////////////////////////////////////////////////////////////////////

void ComputeElecField(double t, const dTensor2& node, const dTensorBC3& qvals, 
                      dTensorBC3& aux, dTensorBC3& Evals)
{

    void ConvertBC(const dTensor2& node, const int mbc1, const int mbc2, 
           const double tmp1, const double tmp2, const dTensorBC3& Fvals, 
           double& gamma, double& beta);
    
    void PoissonSolve(const int mstart, const dTensor2& node, 
        const double gamma, const double beta, 
        const dTensorBC3& Fvals, dTensorBC3& qvals);


    const int melems = qvals.getsize(1);
    const int   meqn = qvals.getsize(2);
    const int   kmax = qvals.getsize(3);
    const int   mbc  = qvals.getmbc();

    dTensorBC3 fvals(melems, 1 , kmax, mbc);  // right hand side function

    // phi(mx,1,kmax1d) = 'electric' field
    // phi(mx,2,kmax1d) = potential  (not unique!)
    dTensorBC3 phi(melems, 2, kmax, mbc);  // Electric field and potential


    for(int i=1; i<= melems; i++) 
       for(int k=1; k<=kmax; k++)
    { 
        fvals.set(i, 1, k, qvals.get(i,1,k) - aux.get(i,1,k) ); 
    }

    //////////////////////////////////////////////////////////////////////////
    // Solve for electric field and potential
    //////////////////////////////////////////////////////////////////////////
    
    // periodic boundary conditions: alpha = beta = 0;
    //       gamma = delta will be forced provided \int f = 0.
    double gamma, beta;
    ConvertBC(node, 1, 1, 0.0, 0.0, fvals, gamma, beta);
    PoissonSolve(1, node, gamma, beta, fvals, phi);

    //////////////////////////////////////////////////////////////////////////
    // Save results into Evals
    //////////////////////////////////////////////////////////////////////////
    for(int i = 1; i <= melems; i++ )
       for(int k=1; k <= kmax; k++ )
         { Evals.set(i, 1, k, phi.get(i, 1, k) ); }

}// end of Compute Electric Field

void ComputeElecPotential(double t, const dTensor2& node, const dTensorBC3& qvals,
                          dTensorBC3& aux, dTensorBC3& phi)
{
    /////////////////////////////////////////////////////////////////////////
    // This function is only ever called in Output.cpp -- the purpose is to
    // print the moments for both the electric potential, as well as the 
    // electric field
    /////////////////////////////////////////////////////////////////////////

    // --------------------- Function Declarations -------------------------//
    void ConvertBC(const dTensor2& node, const int mbc1, const int mbc2, 
           const double tmp1, const double tmp2, const dTensorBC3& Fvals, 
           double& gamma, double& beta);
    
    void PoissonSolve(const int mstart, const dTensor2& node, 
        const double gamma, const double beta, 
        const dTensorBC3& Fvals, dTensorBC3& qvals);
    // --------------------- Function Declarations -------------------------//

    // --------------------- Local Variables -------------------------------//
    const int melems = qvals.getsize(1);
    const int   meqn = qvals.getsize(2);
    const int   kmax = qvals.getsize(3);
    const int   mbc  = qvals.getmbc();

    dTensorBC3 fvals(melems, 2 , kmax, mbc, 1);  // right hand side function

    for(int i=1; i<= melems; i++)
      for(int k=1; k<=kmax; k++)
    { 
        double tmp = qvals.get(i,1,k);
        fvals.set(i, 1, k, tmp - aux.get(i,1,k) ); 
    }

    //////////////////////////////////////////////////////////////////////////
    // Solve for electric field and potential
    //////////////////////////////////////////////////////////////////////////
    
    // periodic boundary conditions: alpha = beta = 0;
    //       gamma = delta will be forced provided \int f = 0.
    double gamma, beta;
    ConvertBC(node, 1, 1, 0.0, 0.0, fvals, gamma, beta);
    PoissonSolve(1, node, gamma, beta, fvals, phi);

}

// This tester function simply checks if the right hand side is a valid right
// hand side ...
void CheckSolvability(const dTensorBC3& fvals)
{
    const double dx = dogParamsCart1.get_dx();
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


