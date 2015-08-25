#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "QuadMomentParams.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;


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
  const string closure = quadMomentParams.closure;

  int cl=0;
  if (closure == "M4")
    {  cl = 1;  }
  else if (closure == "M6")
    {  cl = 2;  }
  else if (closure == "M8")
    {  cl = 3;  }

  switch(cl)
    {
       case 1:
      {
   /////////////////////////////////// M4 closre /////////////////////////////
      //Left states
      const double rhol = Ql.get(1);
      const double ul   = Ql.get(2)/rhol;
      const double energyl = Ql.get(3);
      const double pressl = energyl-rhol*pow(ul,2);
      const double ql   = Ql.get(4)-rhol*pow(ul,3)-3.0e0*pressl*ul;
  
      // Right states
      const double rhor = Qr.get(1);
      const double ur   = Qr.get(2)/rhor;
      const double energyr = Qr.get(3);
      const double pressr = energyr-rhor*pow(ur,2);
      const double qr   = Qr.get(4)-rhor*pow(ur,3)-3.0e0*pressr*ur;

      // Deterministic of Quadrature points
      const double detl  =  sqrt(pressl/rhol + pow(0.5*ql/pressl,2));
      const double detr  =  sqrt(pressr/rhor + pow(0.5*qr/pressr,2));
   
      // Maximum speed
      s1 = ul + ql/2.0e0/pressl - detl;
      s2 = ur + qr/2.0e0/pressr + detr;

       break;
      }

       case 2:
      {
  /////////////////////////////////// M6 closure ////////////////////////////////

      const double miu1l = Auxl.get(1);
      const double miu2l = Auxl.get(2);
      const double miu3l = Auxl.get(3);  
      const double w1l   = Auxl.get(4);
      const double w2l   = Auxl.get(5);
      const double w3l   = Auxl.get(6);
      
      const double miu1r = Auxr.get(1);
      const double miu2r = Auxr.get(2);
      const double miu3r = Auxr.get(3);  
      const double w1r   = Auxr.get(4);
      const double w2r   = Auxr.get(5);
      const double w3r   = Auxr.get(6);

      // Minimum Speed
      s1 = Min(Min(miu1l,miu2l),miu3l);

      // Maximum Speed
      s2 = Max(Max(miu1r,miu2r),miu3r);

       break;
      }

       case 3:
      {
  /////////////////////////////////// M8 closure ////////////////////////////////
      
      const double miu1l = Auxl.get(1);
      const double miu2l = Auxl.get(2);
      const double miu3l = Auxl.get(3); 
      const double miu4l = Auxl.get(4); 
      const double w1l   = Auxl.get(5);
      const double w2l   = Auxl.get(6);
      const double w3l   = Auxl.get(7);
      const double w4l   = Auxl.get(8);
  
      const double miu1r = Auxr.get(1);
      const double miu2r = Auxr.get(2);
      const double miu3r = Auxr.get(3); 
      const double miu4r = Auxr.get(4); 
      const double w1r   = Auxr.get(5);
      const double w2r   = Auxr.get(6);
      const double w3r   = Auxr.get(7);
      const double w4r   = Auxr.get(8);
  
      // Minimum Speed
      s1 = Min(Min(Min(miu1l,miu2l),miu3l),miu4l);

      // Maximum Speed
      s2 = Max(Max(Max(miu1r,miu2r),miu3r),miu4r);

       break;
      }
      
    }  
  
}



