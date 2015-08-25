#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{       
  const double rhol = Ql.get(1);
  const double ul   = Ql.get(2)/rhol;
  const double ul2  = ul*ul;
  const double ul3  = ul2*ul;
  const double pl   = Ql.get(3) - rhol*ul2;
  const double ql   = Ql.get(4) - 3.0*pl*ul - rhol*ul3;
  
  double s1l = ul + 13.0*ql/(6.0*pl) - osq13*sqrt(pl*(33.0+2.0*sqrt(165.0))/rhol);
  double s2l = ul + 13.0*ql/(6.0*pl) + osq13*sqrt(pl*(33.0+2.0*sqrt(165.0))/rhol);

  const double rhor = Qr.get(1);
  const double ur   = Qr.get(2)/rhor;
  const double ur2  = ur*ur;
  const double ur3  = ur2*ur;
  const double pr   = Qr.get(3) - rhor*ur2;
  const double qr   = Qr.get(4) - 3.0*pr*ur - rhor*ur3;

  double s1r = ur + 13.0*qr/(6.0*pr) - osq13*sqrt(pr*(33.0+2.0*sqrt(165.0))/rhor);
  double s2r = ur + 13.0*qr/(6.0*pr) + osq13*sqrt(pr*(33.0+2.0*sqrt(165.0))/rhor);

  s1 = Min(s1l,s1r);
  s2 = Max(s2l,s2r);
}



