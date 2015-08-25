#include "dogdefs.h"
#include "DogParams.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      TRIPLE PERIODIC BOUNDARY CONDITIONS
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
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=(1-mbc); j<=(my+mbc); j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(1-mbc); k<=(mz+mbc); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		   const double tmp = q->get(imod,jmod,kmod,m,ell);
		   q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************

  // ***********************************************
  // RIGHT BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=(mx+1); i<=(mx+mbc); i++)
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=(1-mbc); j<=(my+mbc); j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(1-mbc); k<=(mz+mbc); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		    const double tmp = q->get(imod,jmod,kmod,m,ell);
		    q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************

  // ***********************************************
  // FRONT BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=(1-mbc); j<=(0); j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(1-mbc); k<=(mz+mbc); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		    const double tmp = q->get(imod,jmod,kmod,m,ell);
		    q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************

  // ***********************************************
  // BACK BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=(my+1); j<=(my+mbc); j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(1-mbc); k<=(mz+mbc); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		    const double tmp = q->get(imod,jmod,kmod,m,ell);
		    q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************
  
  // ***********************************************
  // BOTTOM BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=1; j<=my; j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(1-mbc); k<=(0); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		    const double tmp = q->get(imod,jmod,kmod,m,ell);
		    q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************

  // ***********************************************
  // TOP BOUNDARY
  // ***********************************************
#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    {
      const int imod = i + mx*( (i<1) - (i>mx) );

      for (int j=1; j<=my; j++)
	{
	  const int jmod = j + my*( (j<1) - (j>my) );
	  
	  for (int k=(mz+1); k<=(mz+mbc); k++)
	    {
	      const int kmod = k + mz*( (k<1) - (k>mz) );

	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=kmax; ell++)
		  {	      
		    const double tmp = q->get(imod,jmod,kmod,m,ell);
		    q->set(i,j,k,m,ell, tmp );
		  }
	    }
	}
    }
  // ***********************************************
}
