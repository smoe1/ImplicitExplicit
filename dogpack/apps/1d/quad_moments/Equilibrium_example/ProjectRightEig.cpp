#include <cmath>
#include "tensors.h"


// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
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

  // Project onto right eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Qvals.set(1,k, Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(3,k) 
                + Wvals.get(4,k) + Wvals.get(5,k)  );
      
      Qvals.set(2,k, z1*Wvals.get(1,k) + z2*Wvals.get(2,k) + z3*Wvals.get(3,k)
                + z4*Wvals.get(4,k) + z5*Wvals.get(5,k)  );
      
      Qvals.set(3,k, pow(z1,2)*Wvals.get(1,k) + pow(z2,2)*Wvals.get(2,k)
                + pow(z3,2)*Wvals.get(3,k) + pow(z4,2)*Wvals.get(4,k) 
                + pow(z5,2)*Wvals.get(5,k)  );
      
      Qvals.set(4,k, pow(z1,3)*Wvals.get(1,k) + pow(z2,3)*Wvals.get(2,k)
                + pow(z3,3)*Wvals.get(3,k) + pow(z4,3)*Wvals.get(4,k) 
                + pow(z5,3)*Wvals.get(5,k)  );
      
      Qvals.set(5,k, pow(z1,4)*Wvals.get(1,k) + pow(z2,4)*Wvals.get(2,k)
                + pow(z3,4)*Wvals.get(3,k) + pow(z4,4)*Wvals.get(4,k) 
                + pow(z5,4)*Wvals.get(5,k)  );
    }
  

}
