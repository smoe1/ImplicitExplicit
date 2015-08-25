#include "tensors.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Two-fluid plasma equations (5-moment)
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{
    double cs_light;
    
    cs_light = 0.5*(Auxl.get(4)+Auxr.get(4));

    s1 =-cs_light; 
    s2 = cs_light;    
}
