#include "tensors.h"
#include "dogdefs.h"
#include <cmath>

// Function that is called after initial conditions are set
void AfterQinit(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   mbc  = q.getmbc();
  const int   maux = aux.getsize(2);


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


//  initial equilibrium distribution with constant rho

  for(int i=1;i<=melems;i++)
     for(int k=1;k<=kmax;k++)
        {
           // Save rho0 in aux(1:melems,1,1:kmax)
           aux.set(i,1,k,q.get(i,1,k));

        }

//  initial equilibrium distribution with rho = sqrt(2*pi)/2*(2+cos(2*pi*x))
//  L2Project(0,1-mbc,melems+mbc,node,q,aux,aux,&SetBackground);
  
  
}


void SetBackground(const dTensor1& xpts, 
                   const dTensor2& Q, 
                   const dTensor2& qvals, 
                   dTensor2& auxvals)
{
   const int numpts = xpts.getsize();
   for(int i=1;i<=numpts;i++)
     {
        double x = xpts.get(i);
         // initial equilibrium distribution with rho = sqrt(2*pi)/2*(2+cos(2*pi*x))
        double rho0=sqrt(2.0*pi)/1.2661*exp(cos(2*pi*x));
        auxvals.set(i,1,rho0);
     }
}
