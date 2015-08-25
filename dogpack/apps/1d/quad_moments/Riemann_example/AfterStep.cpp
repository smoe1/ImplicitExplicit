#include "tensors.h"
#include <cmath>

// Function that is called after each time step
void AfterStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int   mx   = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   maux = aux.getsize(2); 
  const int   mbc  = q.getmbc();
 

}

