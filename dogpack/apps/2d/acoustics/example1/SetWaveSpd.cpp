#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "AcousticParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
		const dTensor1& Ql, const dTensor1& Qr,
		const dTensor1& Auxl, const dTensor1& Auxr,
		double& s1,double& s2)
{

    // sound speed
    const double c = acousticParams.get_c();

    // Minimum speed
    s1 = -c;
    
    // Maximum speed
    s2 = c;
    
}
