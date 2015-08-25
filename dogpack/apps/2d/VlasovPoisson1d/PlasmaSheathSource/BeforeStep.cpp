#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"

//////////////////////////////////////////////////////////////////////////////
//
//  Compute the Electric Field E(x) from q.  This is saved in 
//  aux(1:mx, 1:my, 2, 1:kmax)
//
//////////////////////////////////////////////////////////////////////////////
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver)
{
  // function declarations:
  void ComputeElecField(double t, const dTensor2& node1d,
                        const dTensorBC4& qvals, dTensorBC3& Evals);
  void ConvertQ1dToQ2d(const int &mopt, int istart, int iend, 
                       int jstart, int jend, 
                       const dTensorBC3& q1d, dTensorBC4& q2d);
  void SetBndValues(dTensorBC4& q, dTensorBC4& aux);

  // local variables
  const int mx      = q.getsize(1);
  const int my      = q.getsize(2);
  const int meqn    = q.getsize(3);
  const int kmax    = q.getsize(4);
  const int mbc     = q.getmbc();
  const int maux    = aux.getsize(3);
  const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
  const int kmax1d  = int(sqrt(mpoints));
  
  dTensorBC3 Evals(mx, meqn, kmax1d, mbc,1);
  dTensorBC4 Evals2d(mx, my, meqn, kmax, mbc);
  dTensor2 node1d(mx+1, 1);
  
  const int me = 1;  // equation number TODO me > 1 ???
  
  // save 1d grid points ( this is used for poisson solve )
  for(int i=1; i<=(mx+1); i++)
    { node1d.set(i, 1, dogParamsCart2.get_xl(i)); }
  
  //ComputeElecField(solver.get_time(), node1d, q, Evals);
  ComputeElecField(solver.get_time_hack(),node1d,q,Evals);
  ConvertQ1dToQ2d(1, 1, mx, 1, my, Evals, Evals2d );
  
  // save electric field into aux array (velocity):
  //  double Emax = 0.0;
  //#pragma omp parallel for
  for(int i=1; i<= mx; i++)
    for(int j=1; j<= my; j++)
      for(int k=1; k<= kmax; k++)
        {
          //      Emax = Max( Emax, fabs( Evals2d.get(i,j,me,k) ) );
          aux.set(i,j,2,k, Evals2d.get(i, j, me, k) );
        }
  SetBndValues(q,aux);
  
}

// Simple wrapper function.  We have no source term in this problem, so no need
// to use void* data
void BeforeStep(double dt,
                dTensorBC4& aux,
                dTensorBC4& q,
                DogSolverCart2& solver,
                void* data)
{
  BeforeStep(dt,aux,q,solver);
}
