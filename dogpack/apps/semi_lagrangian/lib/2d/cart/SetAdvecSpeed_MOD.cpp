#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

///////////////////////////////////////////////////////////////////////////////
// Routine for setting advection speeds for advection equation of the form
// q_t + u(y) q_x + v(x) q_y = 0.
//
// Function also sets the maximum wave speed for each cell.
//
//    The advection speeds sampled at the gaussian quadrature points 
//    1:mpoints1d are provided in:
//
//    u1(my, mpoints1d) = u(y) - the 1d gauss point values for speeds
//    u2(mx, mpoints1d) = v(x) - the 1d gauss point values for speeds
//
///////////////////////////////////////////////////////////////////////////////
void SetAdvecSpeed_MOD(const dTensor2& phi1d, 
		       const dTensorBC4& q,
		       const dTensorBC3* Efield,
		       const dTensorBC3* v1d,
		       dTensorBC3& smax, 
		       const int &vel_dir,
		       dTensorBC2& u1, 
		       dTensorBC2& u2)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int maux = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int kmax1d = dogParams.get_space_order();
  const int mpoints1d = u1.getsize(2);
  const double dx = dogParamsCart2.get_dx();
  const double dy = dogParamsCart2.get_dy();

  switch( vel_dir )
    {
    case 1:
      
      // save the velocity u1
      for(int j=1; j<=my; j++)
	for(int m1d=1; m1d<=mpoints1d; m1d++)
	  {
	    double tmp1 = 0.0;
	    for(int k=1; k<=kmax1d; k++)
	      {
		tmp1 += v1d->get(j,1,k)*phi1d.get(m1d,k);
	      }
	    u1.set(j,m1d, tmp1);                // velocity, u(y)
	  }
      break;
      
    case 2:
      
      // save the velocity u2
      for(int i=1; i<=mx; i++)	
	for(int m1d=1; m1d<=mpoints1d; m1d++)
	  {
	    double tmp2 = 0.0;
	    for(int k=1; k<=kmax1d; k++)
	      {
		tmp2 += Efield->get(i,1,k)*phi1d.get(m1d,k);
	      }
	    
	    u2.set(i,m1d, tmp2);                 // velocity, v(x)
	  }	
      break;
    }
  
  // set maximum wave speed  
  for(int j=1; j<=my; j++)
    {
      double tmp1 = 0.0;
      
      for(int m=1; m<=mpoints1d; m++)
	{  tmp1 = Max(tmp1, fabs( dy*u1.get(j,m) ) );  }

      for(int i=1; i<=mx; i++)
	{  smax.set(i,j,1, tmp1 );  }
    }

  for(int i=1; i<=mx; i++)
    {
      double tmp2 = 0.0;
      
      for(int m=1; m<=mpoints1d; m++)
	{  tmp2 = Max(tmp2, fabs( dx*u2.get(i,m) ) );  }

      for(int j=1; j<=my; j++)
	{  smax.set(i,j,2, tmp2 );  }
    }
  
}// end of function SetAdvecSpeed
///////////////////////////////////////////////////////////////////////////////
