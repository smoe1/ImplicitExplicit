#include "tensors.h"
#include "QuadMomentParams.h"
#include <cmath>
#include <string>
#include <iostream>
using namespace std;

// Function that is called after each time step
void AfterStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int   mx   = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   maux = aux.getsize(2); 
  const int   mbc  = q.getmbc();

  const double epsilon = quadMomentParams.epsilon;
  

  void L2Project(int mopt, int istart, int iend,
        const dTensor2& node,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
        dTensorBC3& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
  void Maxwellian(const dTensor1& xpts,
                const dTensor2& qvals,
                const dTensor2& auxvals,
                dTensor2& fvals);
/*
  dTensorBC3 Feq(mx,meqn,kmax,mbc);
  L2Project(0,1-mbc,mx+mbc,node,q,aux,Feq,&Maxwellian);

  for (int i=1-mbc; i<=mx+mbc; i++)
    for ( int k=1; k<=kmax; k++)
      {
        double qold = q.get(i,4,k);
        double qeq  = Feq.get(i,4,k);
        double temp = qeq + (qold - qeq)*(-dt/epsilon);
        q.set(i,4,k,temp);
      }
*/
  
}

void Maxwellian(const dTensor1& xpts,
                const dTensor2& qvals,
                const dTensor2& auxvals,
                dTensor2& fvals)

{
  const int numpts = xpts.getsize();
  
  const int use_collision = quadMomentParams.use_collision;
  const double epsilon = quadMomentParams.epsilon;
  const string closure = quadMomentParams.closure;

  for (int i=1; i<=numpts; i++)
  {
    /////////////////// M4 closure in 1D case////////////////////
    if(closure=="M4")
      {
        // Variables
        const double rho    = qvals.get(i,1);
        const double u      = qvals.get(i,2)/rho;
        const double energy = qvals.get(i,3);
        const double m3     = qvals.get(i,4);
        const double press  = energy-rho*pow(u,2);

        fvals.set(i,1,  0.0 );
        fvals.set(i,2,  0.0 );
        fvals.set(i,3,  0.0 );
        if (use_collision == 1)
        {
          fvals.set(i,4,  (3.0*u*energy-2.0*rho*pow(u,3)-m3)/epsilon ); 
        }
        else 
        {
          fvals.set(i,4,  0.0);
        }
     }
   else if(closure=="M6")
     {  
        fvals.set(i,1,  0.0 );
        fvals.set(i,2,  0.0 );
        fvals.set(i,3,  0.0 );
        fvals.set(i,4,  0.0 );
        fvals.set(i,5,  0.0 );
        fvals.set(i,6,  0.0 );

     }

  }
}
