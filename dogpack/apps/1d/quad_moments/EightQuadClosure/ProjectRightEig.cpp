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
}
