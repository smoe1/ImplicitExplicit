#include <cmath>
#include "dogdefs.h"
#include "DogSolverCart2.h"

void ElectricField(const dTensor2& xpts, dTensor2& e, void* data);

void ElectricFieldWrap(const dTensor2& xpts, 
		       const dTensor2& NOT_USED_1,
		       const dTensor2& NOT_USED_2,
		       dTensor2& e, 
		       void* data)
{
  ElectricField(xpts,e,data);
}

void L2Project_extra(const int istart,
                     const int iend,
                     const int jstart,
                     const int jend,
                     const int QuadOrder,
                     const int BasisOrder_qin,
                     const int BasisOrder_auxin,
                     const int BasisOrder_fout,    
                     const dTensorBC4* qin,
                     const dTensorBC4* auxin,
                     dTensorBC4* fout,
                     void (*Func)(const dTensor2&,const dTensor2&,
                                  const dTensor2&,dTensor2&, void* data),
                     void* data);

void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver, void* data);
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver);

// Function that is called after each time stage
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver )
{
  int mx   = q.getsize(1);
  int my   = q.getsize(2);
  int meqn = q.getsize(3);
  int kmax = q.getsize(4);
  int mbc  = q.getmbc();
  int maux = aux.getsize(3);
}

// Function that is called after each time stage
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver, void* data)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
  
  //return;
  
  // Compare with exact electric and print error
  dTensorBC4 Evals2d(mx,my,1,kmax,mbc);
  const int istart = 1;
  const int iend   = mx;
  const int jstart = 1;
  const int jend   = 1;
  const int QuadOrder = 5;

  L2Project_extra(istart,iend,jstart,jend,QuadOrder,
		  QuadOrder,QuadOrder,QuadOrder,
		  &q,&aux,&Evals2d, &ElectricFieldWrap, data);

  // do a poisson solve on q
  dTensorBC4 aux_extra(mx,my,maux,kmax,mbc);
  BeforeStep(dt, aux_extra, q, solver);
  
  double tmp = 0;
  for( int i=1; i <= mx; i++)
    for( int k=1; k <= kmax; k++)
      {
        tmp += fabs( Evals2d.get(i,1,1,k) - aux_extra.get(i,1,2,k) );
      }
  printf("   error in electric field = %2.8e\n", tmp);
 
}
