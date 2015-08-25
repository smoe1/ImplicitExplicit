#include "dogdefs.h"
#include "math.h"
#include <cmath>
#include "dog_math.h"
#include "QuadMomentParams.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//    
//

void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();
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
      
      ////////////////////////// M4 closure/////////////////////////////
      for (int i=1; i<=numpts; i++)
	{
	  // Variables
	  const double rho    = Q.get(i,1);
	  const double u      = Q.get(i,2)/rho;
	  const double energy = Q.get(i,3);
	  const double m3     = Q.get(i,4);
	  const double press  = energy-rho*pow(u,2);
	  const double q      = m3-rho*pow(u,3)-3.0e0*press*u;
	  
	  const double miu1 = u + q/2.0/press - sqrt(press/rho+q*q/(4.0*press*press));
	  const double miu2 = u + q/2.0/press + sqrt(press/rho+q*q/(4.0*press*press));
	  const double w1   = rho/2.0 + rho*rho*q/2.0/sqrt(rho*rho*q*q   +4.0*rho*press*press*press);
	  const double w2   = rho/2.0 - rho*rho*q/2.0/sqrt(rho*rho*q*q+4.0*rho*press*press*press);
	  
	  // Flux function
	  flux.set(i,1, Q.get(i,2) );
	  flux.set(i,2, Q.get(i,3) );
	  flux.set(i,3, Q.get(i,4) ); 	
	  flux.set(i,4, w1*pow(miu1,4) + w2*pow(miu2,4) );	
	}
      break;
      
    case 2:
      
      /////////////////////////// M6 closure ///////////////////////////
      for (int i=1; i<=numpts; i++)
	{
          const double miu1 = Aux.get(i,1);
          const double miu2 = Aux.get(i,2);
          const double miu3 = Aux.get(i,3);  
          const double w1   = Aux.get(i,4);
          const double w2   = Aux.get(i,5);
          const double w3   = Aux.get(i,6);
                      
	  flux.set(i,1, Q.get(i,2) );
	  flux.set(i,2, Q.get(i,3) );
	  flux.set(i,3, Q.get(i,4) ); 
	  flux.set(i,4, Q.get(i,5) );
	  flux.set(i,5, Q.get(i,6) );
	  flux.set(i,6, w1*pow(miu1,6)+w2*pow(miu2,6)+w3*pow(miu3,6)); 
        }
      break;     
 
    case 3:

      /////////////////////////// M8 closure ///////////////////////////
      for (int i=1; i<=numpts; i++)
        {
          const double miu1 = Aux.get(i,1);
          const double miu2 = Aux.get(i,2);
          const double miu3 = Aux.get(i,3);
          const double miu4 = Aux.get(i,4);
          const double w1   = Aux.get(i,5);
          const double w2   = Aux.get(i,6);
          const double w3   = Aux.get(i,7);
          const double w4   = Aux.get(i,8);
          /*printf("miu1 is %f\n",miu1);
          printf("miu2 is %f\n",miu2);
          printf("miu3 is %f\n",miu3);
          printf("miu4 is %f\n",miu4);
          printf("w1 is %f\n",w1);
          printf("w2 is %f\n",w2);
          printf("w3 is %f\n",w3);
          printf("w4 is %f\n",w4);*/

          flux.set(i,1, Q.get(i,2) );
	  flux.set(i,2, Q.get(i,3) );
	  flux.set(i,3, Q.get(i,4) ); 
	  flux.set(i,4, Q.get(i,5) );
	  flux.set(i,5, Q.get(i,6) );
          flux.set(i,6, Q.get(i,7) );
          flux.set(i,7, Q.get(i,8) );
          flux.set(i,8, w1*pow(miu1,8)+w2*pow(miu2,8)+w3*pow(miu3,8)+w4*pow(miu4,8));
        }
       break;
          
    }

}



