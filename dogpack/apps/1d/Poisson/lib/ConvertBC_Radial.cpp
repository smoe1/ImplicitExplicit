#include <iostream>
#include "dogdefs.h"

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
//  mbc2 == 1 means tmp2 = beta = u(b)
//  mbc2 == 2 means tmp2 = delta = u_x(b)
//
//  When boundary conditions are full neumman, we set beta = u(b) = 0
//
///////////////////////////////////////////////////////////////////////////////

void ConvertBC_Radial(const dTensor2& node, const int mbc1, const int mbc2, 
        const double tmp1, const double tmp2, const dTensorBC3& Fvals, 
        double& gamma, double& beta) 
{

    const int melems = Fvals.getsize(1);
    const int   meqn = Fvals.getsize(2);
    const int morder = Fvals.getsize(3);
    const double a   = node.get(1,1);
    const double b   = node.get(melems+1,1);
    const double dx  = node.get(2,1)-node.get(1,1);

    //nothing to do in this case
    // u_x(a) = gamma; u(b) = beta
    if(mbc1==2 && mbc2==1)
    {
        gamma = tmp1;
        beta  = tmp2;
        return;
    }

    printf("We haven't implemented anything other than mbc1 == 2\n");
    exit(1);

}//end of ConvertBCs
