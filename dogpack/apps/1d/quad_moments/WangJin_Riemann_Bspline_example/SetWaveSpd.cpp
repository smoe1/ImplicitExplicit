#include <cmath>
#include "dog_math.h"
#include "tensors.h"


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
  double GetAlpha(const double rho,
		  const double p,
		  const double q,
		  const double r);
  
  const double rhol = Ql.get(1);
  const double ul   = Ql.get(2)/rhol;
  const double ul2  = ul*ul;
  const double ul3  = ul2*ul;
  const double ul4  = ul3*ul;
  const double pl   = Ql.get(3) - rhol*ul2;
  const double ql   = Ql.get(4) - 3.0*pl*ul - rhol*ul3;
  const double rl   = Ql.get(5) - 4.0*ql*ul - 6.0*pl*ul2 - rhol*ul4;
  
  double alphal = GetAlpha(rhol,pl,ql,rl);
  double qtl;
  
  if (alphal < 1.0e-5)
    {  qtl = 0.0;  }
  else
    {  qtl = ql/alphal;  }
  
  const double rhor = Qr.get(1);
  const double ur   = Qr.get(2)/rhor;
  const double ur2  = ur*ur;
  const double ur3  = ur2*ur;
  const double ur4  = ur3*ur;
  const double pr   = Qr.get(3) - rhor*ur2;
  const double qr   = Qr.get(4) - 3.0*pr*ur - rhor*ur3;
  const double rr   = Qr.get(5) - 4.0*qr*ur - 6.0*pr*ur2 - rhor*ur4;
  
  double alphar = GetAlpha(rhor,pr,qr,rr);
  double qtr;
  
  if (alphar < 1.0e-4)
    {  qtr = 0.0;  }
  else
    {  qtr = qr/alphar;  }
  
  double const MG = 2.856970013872806; //sqrt((5.0+sqrt(10.0))
  
  s1 = ul + 0.4*qtl/pl - MG*sqrt((1.0-alphal)*pl/rhol);
  s2 = ur + 0.4*qtr/pr + MG*sqrt((1.0-alphar)*pr/rhor);

  /*
  // Deterministic of Quadrature points
  const double detl  =  sqrt(pl/rhol + pow(0.5*ql/pl,2));
  const double detr  =  sqrt(pr/rhor + pow(0.5*qr/pr,2));
   
  // Maximum speed
  s1 = ul + ql/2.0e0/pl - detl;
  s2 = ur + qr/2.0e0/pr + detr;
  */

  //s1 = Min(ul-sqrt(3.0*pl/rhol), ur-sqrt(3.0*pr/rhor));
  //s2 = Min(ul+sqrt(3.0*pl/rhol), ur+sqrt(3.0*pr/rhor));
  
}



