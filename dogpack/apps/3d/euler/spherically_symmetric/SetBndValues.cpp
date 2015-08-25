#include "dogdefs.h"
#include "DogParams.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH-ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues(dTensorBC5& q, 
		  dTensorBC5& aux)
{
  void SetBndValues(dTensorBC5* q);
  SetBndValues(&q);
  if (dogParams.get_maux()>0)
    {  SetBndValues(&aux);  }
}  

void SetBndValues(dTensorBC5* q)
{
  const int mx   = q->getsize(1);
  const int my   = q->getsize(2);
  const int mz   = q->getsize(3);
  const int meqn = q->getsize(4);
  const int kmax = q->getsize(5);
  const int mbc  = q->getmbc();
  
  // ***********************************************
  // LEFT BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(1-mbc); i<=(0); i++)
    for (int j=1; j<=(my); j++)
      for (int k=1; k<=(mz); k++)
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(1,j,k,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }
  // ***********************************************
  

  // ***********************************************
  // RIGHT BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=(my); j++)
      for (int k=1; k<=(mz); k++)    
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(mx,j,k,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }
  // ***********************************************


  // ***********************************************
  // FRONT BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(0); j++)    
      for (int k=1; k<=(mz); k++)
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(i,1,k,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }
  // ***********************************************
  

  // ***********************************************
  // BACK BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(my); j<=(my+mbc); j++)    
      for (int k=1; k<=(mz); k++)
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(i,my,k,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }	
  // ***********************************************


  // ***********************************************
  // BOTTOM BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)    
      for (int k=(1-mbc); k<=(0); k++)
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(i,j,1,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }
  // ***********************************************
  

  // ***********************************************
  // TOP BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)    
      for (int k=(mz+1); k<=(mz+mbc); k++)
        for (int m=1; m<=meqn; m++)
          for (int ell=1; ell<=kmax; ell++)
            {
              double tmp = q->get(i,j,mz,m,ell);
              q->set(i,j,k,m,ell, tmp );
            }
  // ***********************************************
}
