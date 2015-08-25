#include "dogdefs.h"  //2d library defs
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//  Function for converting boundary conditions to gamma and beta
//  
//  Parameters:
//
//  tmp1 = "value" of u at x=a.  This can be either u_x(a) or u(a).  "value" is
//  specified by parameter mbc1.
//
//  tmp2 = "value" of u at x=b.  This can be either u_x(b) or u(b).  "value is
//  specified by parameter mbc2.
//
//  mbc1 == 1 means tmp1 = alpha = u(a)
//  mbc1 == 2 means tmp1 = gamma = u_x(a)
//
//  mbc2 == 1 means tmp2 = beta  = u(b)
//  mbc2 == 2 means tmp2 = delta = u_x(b)
//
//  When boundary conditions are full neumman, we set beta = u(b) = 0
//
///////////////////////////////////////////////////////////////////////////////
void ConvertBC(const dTensor2* node, 
	       const int mbc1, 
	       const int mbc2, 
	       const double tmp1, 
	       const double tmp2, 
	       const dTensorBC3& Fvals, 
	       double& gamma, 
	       double& beta) 
{
  int i;
  int melems = Fvals.getsize(1);
  int   meqn = Fvals.getsize(2);
  int morder = Fvals.getsize(3);
  double a = node->get(1,1);
  double b = node->get(melems+1,1);
  double dx = node->get(2,1)-node->get(1,1);
  double xc,dcount1,dcount2;

  //nothing to do in this case
  // u_x(a) = gamma; u(b) = beta
  if(mbc1==2 && mbc2==1)
    {
      gamma = tmp1;
      beta = tmp2;
      return;
    }
  
  //reverse boundary conditions.
  //u(a) = alpha; u_x(b) = delta
  if(mbc1==1 && mbc2==2)
    {
      dcount1 = tmp2;
      for(i=1; i<=melems; i++) 
        {
	  dcount1 = dcount1 + dx*Fvals.get(i,1,1);
        } 
      gamma = dcount1;
      
      dcount2 = 0.0;
      if (morder==1)
        {
	  for(i=1; i<=melems; i++) 
            {
	      xc = node->get(i,1) + 0.5*dx;
	      dcount2 = dcount2 + dx*(xc-b)*Fvals.get(i,1,1);
            }
        }
      else
        {
	  for(i=1; i<=melems; i++) 
            {
	      xc = node->get(i,1) + 0.5*dx;
	      dcount2 = dcount2 + dx*(xc-b)*Fvals.get(i,1,1)
		+ dx*dx/(2.0*sq3)*Fvals.get(i,1,2);
            }
        }
      beta = tmp1 + (b-a)*gamma +dcount2;
      
      return;
    }
  
  // full Dirichlet boundary conditions
  // u(a) = tmp1 = alpha; u(b) = tmp2 = beta
  if(mbc1==1 && mbc2==1)
    {		
      beta = tmp2;
      
      dcount1 = (tmp2-tmp1)/dx;	
      for(i=1; i<=melems; i++) 
        {
	  xc = node->get(i,1) + 0.5*dx;
	  dcount1 = dcount1 + (b-xc)*Fvals.get(i,1,1);
        }
      dcount1 = dcount1 * dx/(b-a);
      
      dcount2 = 0.0;
      if (morder>1)
        {
	  for(i=1; i<=melems; i++) 
            {
	      dcount2 = dcount2 + Fvals.get(i,1,2);
            }
	  dcount2 = dcount2 *( -dx*dx / (2.0*sq3*(b-a)) );
        }
      
      gamma = dcount1 + dcount2;
      return;
    }
  
  //full Neumman boundary conditions
  //solution is only unique up to an additive constant.
  //this is the solution with u(a) = alpha = 0
  // tmp1 = u_x(a) = gamma;  tmp2 = u_x(b) = delta.
  if(mbc1==2 && mbc2==2)
    {
      gamma = tmp1;
      dcount1 = (b-a)*tmp2;
      for(i=1; i<=melems; i++)
        {
	  xc = node->get(i,1) + 0.5*dx;
	  dcount1 = dcount1 + dx*(xc-a)*Fvals.get(i,1,1);
        }
      
      dcount2 = 0.0;
      if (morder>1)
        {
	  for(i=1; i<=melems; i++)
            {
	      dcount2 = dcount2 + dx*dx/(2.0*sq3)*Fvals.get(i,1,2);
            }
        }
      beta = dcount1 + dcount2;
      return;
    }

  cout << "Error in ConvertBC.cpp.  Didn't find appropriate boundary conditions" << endl;
  cout << "    mbc1 = " << mbc1 << "    mbc2 = " << mbc2 << endl;
  exit(1);
  
}//end of ConvertBCs
