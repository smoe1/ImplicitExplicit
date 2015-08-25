#include "tensors.h"
#include "DogParams.h"
#include "Legendre2d.h"
#include <cart/Limiters.h>
#include "Limiters.h"

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
                  void (*ProjectRightEig)(int,const dTensor1&,
                                          const dTensor1&,const dTensor2&,
                                          dTensor2&),
                  void (*ProjectLeftEig)(int,const dTensor1&,
                                         const dTensor1&,const dTensor2&,
                                         dTensor2&))
{
    // assumes that cell averages satisfy positivity
    //
    ApplyLimiterKrivodonova(aux, q, ProjectRightEig, ProjectLeftEig);

    // positivity limiters are independent of the divergence-free
    // condition, so we can do this before or after projecting
    // onto the locally divergence-free subspace.
    //ApplyPositivityLimiter(q);

    // Divergence-free fix-up (moment-limiting does not preserve div(B)=0, this
    //                         fix-up will restore div(B)=0 preservation)
    if (dogParams.get_use_divfree()==1)
      project_onto_locally_divergence_free_subspace(q);
}
