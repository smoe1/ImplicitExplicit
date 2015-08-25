#include "tensors.h"
#include "TwoFluidParams.h"
#include "dog_math.h"

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
  double Min(double,double);
  double Max(double,double);
  const double light_speed = twofluidParams.light_speed;
  const double gamma = twofluidParams.gamma;

  dTensor1 Qave(10);
  for (int k=1; k<=10; k++)
    {
      Qave.set(k, 0.5*(Qr.get(k)+Ql.get(k)) );
    }

  double rho_i = Qave.get(1);
  double u1_i  = Qave.get(2)/rho_i;
  double u2_i  = Qave.get(3)/rho_i;
  double u3_i  = Qave.get(4)/rho_i;
  double En_i  = Qave.get(5);
  
  double rho_e = Qave.get(6);
  double u1_e  = Qave.get(7)/rho_e;
  double u2_e  = Qave.get(8)/rho_e;
  double u3_e  = Qave.get(9)/rho_e;
  double En_e  = Qave.get(10);      
  
  double press_i  = (gamma-1.0e0)*(En_i - 0.5e0*rho_i*(u1_i*u1_i + u2_i*u2_i + u3_i*u3_i));
  
  double press_e  = (gamma-1.0e0)*(En_e - 0.5e0*rho_e*(u1_e*u1_e + u2_e*u2_e + u3_e*u3_e));

  double cs_i = sqrt(gamma*press_i/rho_i);
  double cs_e = sqrt(gamma*press_e/rho_e);

  s1 = -2.85;//Min(-light_speed, Min(u1_i-cs_i, u1_e-cs_e) ); 
  s2 =  2.85;//Max( light_speed, Max(u1_i+cs_i, u1_e+cs_e) );
}
