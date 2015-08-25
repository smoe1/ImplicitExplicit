#include <cmath>
#include "tensors.h"


// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
		    const dTensor1& Q_ave,
		    const dTensor2& Qvals,
		    dTensor2& Wvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  
  double GetAlpha(const double rho,
		  const double p,
		  const double q,
		  const double r);
   
  // Average states
  const double rho = Q_ave.get(1);
  const double u   = Q_ave.get(2)/rho;
  const double u2  =  u*u;
  const double u3  = u2*u;
  const double u4  = u3*u;
  const double p   = Q_ave.get(3)-rho*u2;
  const double q   = Q_ave.get(4)-3.0*p*u-rho*u3;
  const double r   = Q_ave.get(5)-4.0*q*u-6.0*p*u2-rho*u4;
  
  const double alpha = GetAlpha(rho,p,q,r);
  
  double qt;      
  if (alpha < 1.0e-12)
    {  qt = 0.0;  }
  else
    {  qt = q/alpha;  }
  
  
  double z1 = u + 0.4*qt/p - sqrt((5.0+sqrt(10.0))*p*(1.0-alpha)/rho);
  double z2 = u + 0.4*qt/p - sqrt((5.0-sqrt(10.0))*p*(1.0-alpha)/rho);
  double z3 = u + 0.4*qt/p;
  double z4 = u + 0.4*qt/p + sqrt((5.0-sqrt(10.0))*p*(1.0-alpha)/rho);
  double z5 = u + 0.4*qt/p + sqrt((5.0+sqrt(10.0))*p*(1.0-alpha)/rho);
  
  double d1 = (z1-z2)*(z1-z3)*(z1-z4)*(z1-z5);
  double d2 = (z2-z1)*(z2-z3)*(z2-z4)*(z2-z5);
  double d3 = (z3-z1)*(z3-z2)*(z3-z4)*(z3-z5);
  double d4 = (z4-z1)*(z4-z2)*(z4-z3)*(z4-z5);
  double d5 = (z5-z1)*(z5-z2)*(z5-z3)*(z5-z4);
  
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Wvals.set(1,k, z2*z3*z4*z5/d1*Qvals.get(1,k) 
                - (z2*z3*z4+z2*z3*z5+z2*z4*z5+z3*z4*z5)/d1*Qvals.get(2,k)
                + (z2*z3+z2*z4+z2*z5+z3*z4+z3*z5+z4*z5)/d1*Qvals.get(3,k)
                - (z2+z3+z4+z5)/d1*Qvals.get(4,k) + Qvals.get(5,k)/d1 );
      
      Wvals.set(2,k, z1*z3*z4*z5/d2*Qvals.get(1,k)
                - (z1*z3*z4+z1*z3*z5+z1*z4*z5+z3*z4*z5)/d2*Qvals.get(2,k)
                + (z1*z3+z1*z4+z1*z5+z3*z4+z3*z5+z4*z5)/d2*Qvals.get(3,k)
                - (z1+z3+z4+z5)/d2*Qvals.get(4,k) + Qvals.get(5,k)/d2 );
      
      Wvals.set(3,k, z1*z2*z4*z5/d3*Qvals.get(1,k)
                - (z1*z2*z4+z1*z2*z5+z1*z4*z5+z2*z4*z5)/d3*Qvals.get(2,k)
                + (z1*z2+z1*z4+z1*z5+z2*z4+z2*z5+z4*z5)/d3*Qvals.get(3,k)
                - (z1+z2+z4+z5)/d3*Qvals.get(4,k) + Qvals.get(5,k)/d3 );
      
      Wvals.set(4,k, z1*z2*z3*z5/d4*Qvals.get(1,k)
                - (z1*z2*z3+z1*z2*z5+z1*z3*z5+z2*z3*z5)/d4*Qvals.get(2,k)
                + (z1*z2+z1*z3+z1*z5+z2*z3+z2*z5+z3*z5)/d4*Qvals.get(3,k)
                - (z1+z2+z3+z5)/d4*Qvals.get(4,k) + Qvals.get(5,k)/d4 );
      
      Wvals.set(5,k, z1*z2*z3*z4/d5*Qvals.get(1,k)
                - (z1*z2*z3+z1*z2*z4+z1*z3*z4+z2*z3*z4)/d5*Qvals.get(2,k)
                + (z1*z2+z1*z3+z1*z4+z2*z3+z2*z4+z3*z4)/d5*Qvals.get(3,k)
                - (z1+z2+z3+z4)/d5*Qvals.get(4,k) + Qvals.get(5,k)/d5 );
    }
  
   
}
