#include "dogdefs.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Acoustics equations
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{
  // Minimum speed
  s1 = -Auxl.get(1);
  
  // Maximum speed
  s2 = Auxr.get(1);
}
