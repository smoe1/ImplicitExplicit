#include "tensors.h"
#include "DogState1d.h" 
#include "dogdefs.h"
#include <cmath>

// Function that is called before initial conditions are set
// Initial density and background density are set here
void BeforeQinit(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   mbc  = q.getmbc();
  const int   maux = aux.getsize(2);

  dTensorBC3  Evals(melems, 1, kmax, mbc); 

  void L2Project(int mopt, int istart, int iend,
                 const dTensor2& node,
                 const dTensorBC3& qin, 
                 const dTensorBC3& auxin,  
                 dTensorBC3& Fout,
                 void (*Func)(const dTensor1&, const dTensor2&, 
                              const dTensor2&, dTensor2&));

  void SetBackground(const dTensor1& xpts, 
                     const dTensor2& Q, 
                     const dTensor2& qvals, 
                     dTensor2& auxvals);

  void InitDensity(const dTensor1& xpts, 
                   const dTensor2& Q, 
                   const dTensor2& auxvals, 
                   dTensor2& qvals);
  

  void ComputeElecField(double t, const dTensor2& node, const dTensorBC3& qvals, 
                        dTensorBC3& aux, dTensorBC3& Evals);


   // save background density into aux(1:melems,1,1:kmax):
  L2Project(0,1-mbc,melems+mbc,node,q,aux,aux,&SetBackground);


   // Initialize density
  L2Project(0,1-mbc,melems+mbc,node,q,aux,q,&InitDensity);

  
   // save electric field into aux(1:melems,2,1:kmax):
  ComputeElecField(fetch_dogState().get_time(), node, q, aux, Evals);
    
  for(int i=1; i<= melems; i++)
    for(int k=1; k<= kmax; k++)
      {
         aux.set(i,2,k, -Evals.get(i,1,k) );
      }
  
}

void InitDensity(const dTensor1& xpts, 
                 const dTensor2& NOT_USED_1, 
                 const dTensor2& NOT_USED_2, 
                 dTensor2& qvals)
{
   const int numpts = xpts.getsize();
   for(int i=1;i<=numpts;i++)
     {
        double x = xpts.get(i);
        double rho; 

         // equilibrium distribution with rho = sqrt(2*pi)/2*(2+cos(2*pi*x))
        rho = 0.5*sqrt(2.0*pi)*(2.0+cos(2.0*pi*x));

        qvals.set(i,1,rho);
     }
}

void SetBackground(const dTensor1& xpts, 
                   const dTensor2& NOT_USED_1, 
                   const dTensor2& NOT_USED_2, 
                   dTensor2& auxvals)
{
   const int numpts = xpts.getsize();
   for(int i=1;i<=numpts;i++)
     {
        double x = xpts.get(i);
        double rho0;

         // equilibrium distribution with rho = sqrt(2*pi)/2*(2+cos(2*pi*x))
        rho0=(sqrt(2.0*pi)/1.2660658777520083)*exp(cos(2.0*pi*x));

        auxvals.set(i,1,rho0);
     }
}

