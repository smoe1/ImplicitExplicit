#include "dogdefs.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
{    
  int m,k;
  int meqn = Qvals.getsize(1);
  int kmax = Qvals.getsize(2)+1;

  // Project onto left eigenvectors
  for (k=1; k <= kmax; k++)
    {
      Qvals.set(1,k, (Wvals.get(2,k)-Wvals.get(1,k) ) );
      
      Qvals.set(2,k, (Wvals.get(2,k)+Wvals.get(1,k) ) );
    }

}
