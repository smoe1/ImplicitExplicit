#include "tensors.h"
#include "dog_math.h"
#include <iostream>

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

  /*
  // Gas constant
  double const gamma = 5.0/3.0;
    
  // --------------------
  // ION WAVES
  // --------------------

  // Left states
  double rhol = Ql.get(1);
  double u1l  = Ql.get(2)/rhol;
  double u2l  = Ql.get(3)/rhol;
  double u3l  = Ql.get(4)/rhol;
  double energyl = Ql.get(5);
  double pressl = (gamma-1.0e0)*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l+u3l*u3l));

  // Right states
  double rhor = Qr.get(1);
  double u1r  = Qr.get(2)/rhor;
  double u2r  = Qr.get(3)/rhor;
  double u3r  = Qr.get(4)/rhor;
  double energyr = Qr.get(5);
  double pressr = (gamma-1.0e0)*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r+u3r*u3r));
  
  // Average states
  double rho    = 0.5e0*(rhol+rhor);
  double u1     = 0.5e0*(u1l+u1r);
  double u2     = 0.5e0*(u2l+u2r);
  double u3     = 0.5e0*(u3l+u3r);
  double press  = 0.5e0*(pressl+pressr);
   
  // Sound speeds
  double cl = sqrt(fabs(gamma*pressl/rhol));
  double cr = sqrt(fabs(gamma*pressr/rhor));
  double c  = sqrt(fabs(gamma*press/rho));
  
  // Minimum speed
  double s1ion = Min(u1l-cl,u1-c);
  
  // Maximum speed
  double s2ion = Max(u1r+cr,u1+c);    


  // --------------------
  // ELECTRON WAVES
  // --------------------

  // Left states
  rhol = Ql.get(1);
  u1l  = Ql.get(2)/rhol;
  u2l  = Ql.get(3)/rhol;
  u3l  = Ql.get(4)/rhol;
  energyl = Ql.get(5);
  pressl = (gamma-1.0e0)*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l+u3l*u3l));

  // Right states
  rhor = Qr.get(1);
  u1r  = Qr.get(2)/rhor;
  u2r  = Qr.get(3)/rhor;
  u3r  = Qr.get(4)/rhor;
  energyr = Qr.get(5);
  pressr = (gamma-1.0e0)*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r+u3r*u3r));
  
  // Average states
  rho    = 0.5e0*(rhol+rhor);
  u1     = 0.5e0*(u1l+u1r);
  u2     = 0.5e0*(u2l+u2r);
  u3     = 0.5e0*(u3l+u3r);
  press  = 0.5e0*(pressl+pressr);
   
  // Sound speeds
  cl = sqrt(fabs(gamma*pressl/rhol));
  cr = sqrt(fabs(gamma*pressr/rhor));
  c  = sqrt(fabs(gamma*press/rho));
  
  // Minimum speed
  double s1electron = Min(u1l-cl,u1-c);
  
  // Maximum speed
  double s2electron = Max(u1r+cr,u1+c);    
  */
  
  // --------------------
  // LIGHT WAVES
  // --------------------
  double cs_light = 0.5*(Auxl.get(4)+Auxr.get(4));

  // Minimum speed
  double s1light =-cs_light; 
  
  // Maximum speed
  double s2light = cs_light;   


  // --------------------
  // Find MIN and MAX
  // --------------------
  s1 = s1light;//Min(Min(s1ion,s1electron),s1light);
  s2 = s2light;//Max(Max(s2ion,s2electron),s2light);

}
