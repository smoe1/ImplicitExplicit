#include "dogdefs.h"
#include "DogParams.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
//
//      ZEROTH ORDER EXTRAPOLATION IN Y DIRECTION
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
  void SetBndValues(dTensorBC4& q);
  SetBndValues(q);
  if (dogParams.get_maux()>0)
    {  SetBndValues(aux);  }
}

void SetBndValues(dTensorBC4& q)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();

  // -----------------------
  // BOUNDARY CONDITION LOOP
  // -----------------------

  // ***********************************************
  // BOTTOM BOUNDARY
  // ***********************************************
  for (int j=0; j>=(1-mbc); j--)
    for (int i=1; i<=mx; i++)
      for (int m=1; m<=meqn; m++)
	for (int ell=1; ell<=kmax; ell++)
	  {
	    double tmp = q.get(i,-j+1,m,ell);
	    q.set(i,j,m,ell, tmp );
	  }                           
  // ***********************************************      
      

  // ***********************************************
  // TOP BOUNDARY
  // ***********************************************
  for (int j=(my+1); j<=(my+mbc); j++)
    for (int i=1; i<=mx; i++)
      for (int m=1; m<=meqn; m++)
	for (int ell=1; ell<=kmax; ell++)
	  {
	    double tmp = q.get(i, 2*my-j+1,m,ell);
	    q.set(i,j,m,ell, tmp );
	  }
  // ***********************************************
  

  // ***********************************************
  // LEFT BOUNDARY (Periodic)
  // ***********************************************
  for (int i=0; i>=(1-mbc); i--)
    for (int j=(1-mbc); j<=(my+mbc); j++)
      for (int m=1; m<=meqn; m++)
	for (int ell=1; ell<=kmax; ell++)
	  {
	    double tmp = q.get(i+mx,j,m,ell);
	    q.set(i,j,m,ell, tmp );
	  }
  // ***********************************************
  
  
  // ***********************************************
  // RIGHT BOUNDARY (Periodic)
  // ***********************************************
  for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)
      for (int m=1; m<=meqn; m++)
	for (int ell=1; ell<=kmax; ell++)
	  {
	    double tmp = q.get(i-mx,j,m,ell);
	    q.set(i,j,m,ell, tmp );
	  }
  // ***********************************************
  
}
