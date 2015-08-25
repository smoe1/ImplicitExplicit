#include "dog_math.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Simple advection equation
//
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
		const dTensor1& Ql, const dTensor1& Qr,
                const dTensor1& Auxl, const dTensor1& Auxr,
		double& s1, double& s2)
{
    double un_left,un_right,un_av;

    // Minimum speed
    s1 = 1.0;
    
    // Maximum speed
    s2 = 1.0;
    
}
