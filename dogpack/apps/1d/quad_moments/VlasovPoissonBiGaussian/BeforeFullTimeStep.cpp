#include "tensors.h"
#include "DogParams.h"
//#include "DogState1d.h"
#include "QuadMomentParams.h" 
#include <cmath>

// Function that is called before a full time step
void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold, dTensorBC3& q)
{
  // Compute the Electric Field E(x) from q. This is saved in aux(1:melems, 2, 1:kmax)  
  
    void ComputeElecField(double t, const dTensor2& node,
                          const dTensorBC3& qvals, dTensorBC3& aux, dTensorBC3& Evals);
    void SetBndValues(const dTensor2& node, dTensorBC3& q, dTensorBC3& aux);

    const int melems  = q.getsize(1);
    const int meqn    = q.getsize(2);
    const int kmax    = q.getsize(3);
    const int mbc     = q.getmbc();
    const int maux    = aux.getsize(2);

    dTensorBC3 Evals(melems, 1, kmax, mbc);

printf("ERROR in BeforeFullTimeStep.  What time should we use?\n");
exit(1);
    ComputeElecField(fetch_dogState().get_time(), node, q, aux, Evals);

    // save electric field into aux array (velocity):
    for(int i=1; i<= melems; i++)
      for(int k=1; k<= kmax; k++)
        {
           aux.set(i,2,k, -Evals.get(i,1,k) );
        }
    SetBndValues(node, q, aux);

/////////////////////////////Solve ODE with dt/2////////////////////////////

    void L2Project(int mopt, int istart, int iend,
                   const dTensor2& node,
                   const dTensorBC3& qin, 
                   const dTensorBC3& auxin,  
                   dTensorBC3& Fout,
                   void (*Func)(const dTensor1&, const dTensor2&, 
                                const dTensor2&, dTensor2&));

    void FokkerPlanck(const dTensor1& xpts, 
                     const dTensor2& Q, 
                     const dTensor2& auxvals, 
                     dTensor2& qvals);

    L2Project(0,1-mbc,melems+mbc,node,q,aux,q,&FokkerPlanck);
  

}
