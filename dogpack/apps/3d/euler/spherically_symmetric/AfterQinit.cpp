#include "dogdefs.h"
#include "DogSolverCart3.h"
#include "EulerParams.h"
#include "dog_math.h"

void AfterQinit(DogSolverCart3& solver)
{
  dTensorBC5&  q = solver.fetch_state().fetch_q();
  const int   mx = q.getsize(1);
  const int   my = q.getsize(2);
  const int   mz = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);

  // Don't do anything in first-order case
  if (kmax==1)
    { return; }

  // points where to test positivity
  const int N  = 5;
  const int N3 = N*N*N;
  dTensor1 pts1d(N);
  dTensor2 phi(N3,kmax);

  // take 5-point Gauss-Lobatto rule
  pts1d.set(1, -1.0      );
  pts1d.set(2, -sq3*osq7 );
  pts1d.set(3,  0.0      );
  pts1d.set(4,  sq3*osq7 );
  pts1d.set(5,  1.0      );

  int nsum = 0;
  for (int i=1; i<=N; i++)
    for (int j=1; j<=N; j++)
      for (int k=1; k<=N; k++)
	{
	  nsum = nsum + 1;

	  const double xi   = pts1d.get(i);
	  const double eta  = pts1d.get(j);
	  const double zeta = pts1d.get(k);
	  const double xi2  = xi*xi;
	  const double xi3  = xi*xi2;
	  const double xi4  = xi*xi3;
	  const double eta2 = eta*eta;
	  const double eta3 = eta*eta2;
	  const double eta4 = eta*eta3;
	  const double zeta2 = zeta*zeta;
	  const double zeta3 = zeta*zeta2;
	  const double zeta4 = zeta*zeta3;
      
	  switch( kmax )
	    {
	    case 20:  // fourth order
	      phi.set( nsum,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
	      phi.set( nsum,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
	      phi.set( nsum,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
	      phi.set( nsum,17, 3.0*sq3*xi*eta*zeta                 );
	      phi.set( nsum,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
	      phi.set( nsum,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
	      phi.set( nsum,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
	      phi.set( nsum,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
	      phi.set( nsum,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
	      phi.set( nsum,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
	  
	    case 10:  // third order
	      phi.set( nsum,10, onehalf*sq5*(3.0*zeta2-1.0)         );
	      phi.set( nsum,9,  onehalf*sq5*(3.0*eta2-1.0)          );
	      phi.set( nsum,8,  onehalf*sq5*(3.0*xi2-1.0)           );
	      phi.set( nsum,7,  3.0*eta*zeta                        );
	      phi.set( nsum,6,  3.0*xi*zeta                         );
	      phi.set( nsum,5,  3.0*xi*eta                          );
	      
	    case 4:  // second order  
	      phi.set( nsum,4,  sq3*zeta                            );
	      phi.set( nsum,3,  sq3*eta                             );
	      phi.set( nsum,2,  sq3*xi                              );
	      
	    case 1:  // first order
	      phi.set( nsum,1, 1.0                                  );
	      
	      break;
	      
	    default:
	      unsupported_value_error(kmax);
	    }
	}

  const double gamma = eulerParams.gamma;
  const double MIN_RHO = 1.0;
  const double MAX_RHO = 2.0;
  const double MIN_ENG = 1.0/(gamma-1.0);
  const double MAX_ENG = 2.0/(gamma-1.0);

#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
      for (int k=1; k<=mz; k++)
	{
	  double qavR = q.get(i,j,k,1,1);
	  double qavE = q.get(i,j,k,5,1);

	  double minR = 0.0;
	  double maxR = 0.0; 
	  double minE = 0.0;
	  double maxE = 0.0;
	  
	  for (int m=1; m<=kmax; m++)
	    {
	      minR = minR + phi.get(1,m)*q.get(i,j,k,1,m);
	      minE = minE + phi.get(1,m)*q.get(i,j,k,5,m);
	    }
	  maxR = minR;
	  maxE = minE;
	  
	  for (int n=2; n<=N3; n++)
	    {
	      double valR = 0.0;
	      double valE = 0.0;
	      
	      for (int m=1; m<=kmax; m++)
		{
		  valR = valR + phi.get(n,m)*q.get(i,j,k,1,m);
		  valE = valE + phi.get(n,m)*q.get(i,j,k,5,m);
		}
	      
	      minR = Min(valR,minR);
	      maxR = Max(valR,maxR);

	      minE = Min(valE,minE);
	      maxE = Max(valE,maxE);	      
	    }

	  double thetaR = Min(Min(fabs((MAX_RHO-qavR)/(maxR-qavR)),
				  fabs((MIN_RHO-qavR)/(minR-qavR))),1.0);
	  double thetaE = Min(Min(fabs((MAX_ENG-qavE)/(maxE-qavE)),
				  fabs((MIN_ENG-qavE)/(minE-qavE))),1.0);
	  double theta = Min(thetaR,thetaE);

	  if (theta<1.0)
	    {
	      for (int me=1; me<=meqn; me++)
		for (int m=2; m<=kmax; m++)
		  {
		    q.set(i,j,k,me,m, theta*q.get(i,j,k,me,m) );
		  }		  
	    }	  
	}
}
