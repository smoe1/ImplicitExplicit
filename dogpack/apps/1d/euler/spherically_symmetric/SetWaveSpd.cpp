#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Euler equations for gas dynamics
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,
		double& s2)
{ 
  // Gas constant
  double const gamma = eulerParams.gamma;
  
  // Left states
  double const rhol = Ql.get(1);
  double const u1l  = Ql.get(2)/rhol;
  double const energyl = Ql.get(3);
  double const pressl = (gamma-1.0e0)*(energyl-0.5e0*rhol*(u1l*u1l));
  
  // Right states
  double const rhor = Qr.get(1);
  double const u1r  = Qr.get(2)/rhor;
  double const energyr = Qr.get(3);
  double const pressr = (gamma-1.0e0)*(energyr-0.5e0*rhor*(u1r*u1r));
  
  // Average states
  double const rho    = 0.5e0*(rhol+rhor);
  double const u1     = 0.5e0*(u1l+u1r);
  double const press  = 0.5e0*(pressl+pressr);
   
  // Sound speeds
  double const cl = sqrt(fabs(gamma*pressl/rhol));
  double const cr = sqrt(fabs(gamma*pressr/rhor));
  double const c  = sqrt(fabs(gamma*press/rho));
  
  // Minimum speed
  s1 = Min(u1l-cl,u1-c);
  
  // Maximum speed
  s2 = Max(u1r+cr,u1+c);
  
}
