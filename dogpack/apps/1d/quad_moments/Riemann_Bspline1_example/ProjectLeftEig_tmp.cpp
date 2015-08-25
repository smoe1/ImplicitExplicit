#include <cmath>
#include "tensors.h"
#include "constants.h"

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
   
  // Average states
  const double rho = Q_ave.get(1);
  const double u   = Q_ave.get(2)/rho;
  const double u2  =  u*u;
  const double u3  = u2*u;
  const double p   = Q_ave.get(3)-rho*u2;
  const double q   = Q_ave.get(4)-3.0*p*u-rho*u3;  
  const double alpha = 3.0/13.0;
  
  double z1 = u + (13.0/6.0)*q/p - osq13*sqrt(p/rho*(33.0 + 2.0*sqrt(165.0)));
  double z2 = u + (13.0/6.0)*q/p - osq13*sqrt(p/rho*(33.0 - 2.0*sqrt(165.0)));
  double z3 = u + (13.0/6.0)*q/p + osq13*sqrt(p/rho*(33.0 - 2.0*sqrt(165.0)));
  double z4 = u + (13.0/6.0)*q/p + osq13*sqrt(p/rho*(33.0 + 2.0*sqrt(165.0)));
  
  double detR = (-z4+z3)*(-z4+z2)*(-z3+z2)*(z1-z4)*(-z3+z1)*(-z2+z1);
  dTensor2 Lmat(4,4);

  Lmat.set(1,1, -z3*z4*z2*(-z4 + z3)*(-z4 + z2)*(-z3 + z2)/detR                );
  Lmat.set(1,2, (z2*z3 + z2*z4 + z4*z3)*(-z4 + z3)*(-z4 + z2)*(-z3 + z2)/detR  );
  Lmat.set(1,3, -(z2 + z3 + z4)*(-z4 + z3)*(-z4 + z2)*(-z3 + z2)/detR          );
  Lmat.set(1,4, (-z4 + z3)*(-z4 + z2)*(-z3 + z2)/detR                          );

  Lmat.set(2,1, z3*z4*z1*(-z3 + z1)*(z1 - z4)*(-z4 + z3)/detR                  );
  Lmat.set(2,2, -(z3*z1 + z4*z1 + z4*z3)*(-z3 + z1)*(z1 - z4)*(-z4 + z3)/detR  );
  Lmat.set(2,3, (z1 + z3 + z4)*(-z3 + z1)*(z1 - z4)*(-z4 + z3)/detR            );
  Lmat.set(2,4, -(-z3 + z1)*(z1 - z4)*(-z4 + z3)/detR                          );

  Lmat.set(3,1, -z4*z2*z1*(-z2 + z1)*(z1 - z4)*(-z4 + z2)/detR                 );
  Lmat.set(3,2, (z2*z1 + z4*z1 + z2*z4)*(-z2 + z1)*(z1 - z4)*(-z4 + z2)/detR   );
  Lmat.set(3,3, -(z1 + z2 + z4)*(-z2 + z1)*(z1 - z4)*(-z4 + z2)/detR           );
  Lmat.set(3,4, (-z2 + z1)*(z1 - z4)*(-z4 + z2)/detR                           );

  Lmat.set(4,1, z3*z2*z1*(-z2 + z1)*(-z3 + z1)*(-z3 + z2)/detR                 );
  Lmat.set(4,2, -(z3*z1 + z2*z1 + z2*z3)*(-z2 + z1)*(-z3 + z1)*(-z3 + z2)/detR );
  Lmat.set(4,3, (z1 + z3 + z2)*(-z2 + z1)*(-z3 + z1)*(-z3 + z2)/detR           );
  Lmat.set(4,4, -(-z2 + z1)*(-z3 + z1)*(-z3 + z2)/detR                         );

  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    for (int m=1; m<=4; m++)
      {
  	Wvals.set(m,k, Lmat.get(m,1)*Qvals.get(1,k) + Lmat.get(m,2)*Qvals.get(2,k)
		  + Lmat.get(m,3)*Qvals.get(3,k) + Lmat.get(m,4)*Qvals.get(4,k) );
      }
  
   
}
