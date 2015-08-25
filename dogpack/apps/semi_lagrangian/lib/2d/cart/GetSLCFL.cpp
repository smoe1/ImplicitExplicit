#include "tensors.h"
#include "dog_math.h"
#include "DogParamsCart2.h"

// Wrapper function.  Applications / libraries that with to replace this should
// do so in their own Makefiles
double GetSLCFL(double dt, const dTensorBC4& aux, const dTensorBC3& smax)
{

    double GetCFL(double dt, const dTensorBC4& aux, const dTensorBC3& smax);
    return GetCFL(dt, aux, smax);

}
