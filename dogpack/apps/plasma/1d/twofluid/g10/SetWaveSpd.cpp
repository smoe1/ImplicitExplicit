#include "tensors.h"
#include "Params.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Two-fluid plasma equations (10-moment)
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{
    // should really compute the sound wave speeds to
    // ensure that they do not exceed the light speed

    s1 =-appParams.cs_light; 
    s2 = appParams.cs_light;    
}
