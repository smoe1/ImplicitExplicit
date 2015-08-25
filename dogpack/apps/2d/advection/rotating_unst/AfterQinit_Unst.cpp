#include "dogdefs.h"
#include "mesh.h"
#include "DogParams.h"

// Function that is called after the intial conditions are set
void AfterQinit_Unst(const mesh& Mesh, dTensor3& aux, dTensor3& q)
{

    // ---------------------------------------------------------
    // Post-processing of initial conditions
    // ---------------------------------------------------------

    void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q);
    if( dogParams.using_moment_limiter() )
    { ApplyPosLimiter_Unst(Mesh, aux, q); }
    // ---------------------------------------------------------

}
