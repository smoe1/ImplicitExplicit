#include "tensors.h"
#include "DogParams.h"
#include "DogState1d.h"
#include "QuadMomentParams.h" 
#include <cmath>

// Function that is called after a full time step (i.e., after all stages are complete)

void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
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
    ComputeElecField(fetch_dogState().get_time(), node, q, aux, Evals);

    // save electric field into aux array (velocity):
    for(int i=1; i<=melems; i++)
      for(int k=1; k<=kmax; k++)
        {
	  aux.set(i,2,k, -Evals.get(i,1,k) );
        }
    SetBndValues(node, q, aux);

/////////////////////////////Solve ODE with dt/2 ////////////////////////////

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

    dTensorBC3 qstar(melems,meqn,kmax,mbc);
    L2Project(0,1-mbc,melems+mbc,node,q,aux,qstar,&FokkerPlanck);
    for (int i=(1-mbc); i<=(melems+mbc); i++)
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    q.set(i,m,k, qstar.get(i,m,k) );
	  }

}

void FokkerPlanck(const dTensor1& xpts, 
                 const dTensor2& Q, 
                 const dTensor2& auxvals, 
                 dTensor2& qvals)
{
    const int numpts = xpts.getsize();

    const double epsilon = quadMomentParams.epsilon;
    const double dt = fetch_dogState().get_dt();

    for (int i=1; i<=numpts; i++)
      {
        const double M0 = Q.get(i,1);
        const double M1 = Q.get(i,2);
        const double M2 = Q.get(i,3);
        const double M3 = Q.get(i,4);
        const double M4 = Q.get(i,5);
        const double E  = auxvals.get(i,2);
        const double E2 = E*E;
        const double E3 = E2*E;
        const double E4 = E3*E;
        const double Z  = exp(-0.5*dt/epsilon);
        //printf("dt = %22.16e \n",dt);
        //printf("epsilon = %f\n",epsilon);
        //printf("Z = %f\n",Z);
        const double Z2 = Z*Z;
        const double Z3 = Z2*Z;
        const double Z4 = Z3*Z;

        qvals.set(i,1,M0);
        qvals.set(i,2,Z*(E*M0+M1)-E*M0);
        qvals.set(i,3,Z2*((E2-1.0)*M0+2.0*E*M1+M2)-2.0*E*Z*(E*M0+M1)+(1.0+E2)*M0);
        qvals.set(i,4,Z3*(E*(E2-3.0)*M0+3.0*(E2-1.0)*M1+3.0*E*M2+M3)
                      -3.0*E*Z2*((E2-1.0)*M0+2.0*E*M1+M2)
                      +3.0*(E2+1.0)*Z*(E*M0+M1)-E*(E2+3.0)*M0);
        qvals.set(i,5,Z4*((E4-6.0*E2+3.0)*M0+4.0*E*(E2-3.0)*M1+6.0*M2*(E2-1.0)+4.0*E*M3+M4)
                      -4.0*E*Z3*(E*(E2-3.0)*M0+3.0*M1*(E2-1.0)+3.0*E*M2+M3)
                      +6.0*(E2+1.0)*Z2*((E2-1.0)*M0+2.0*E*M1+M2)
                      -4.0*E*(E2+3.0)*Z*(E*M0+M1)+(E4+6.0*E2+3.0)*M0);
     
      }   
}

