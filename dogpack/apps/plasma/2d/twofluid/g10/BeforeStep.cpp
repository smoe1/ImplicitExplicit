#include "tensors.h"
#include "DogSolverCart2.h"
#include "Limiters.h"

// Function that is called before each time stage
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver)
{
  eprintf("you forgot to override me!");

  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);

  // positivity of cell averages is assumed to hold at this point

  // enforce that positivity holds at the Gauss-Lobatto quadrature points
  PositivityLimiter limiter(q,PositivityPointType::FLUX);
  limiter.apply();
  // otherwise at least enforce positivity at Riemann points ...
}
