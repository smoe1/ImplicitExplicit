#include <stdlib.h>
#include <math.h>
#include "dogdefs.h"
#include "tensors.h"
#include "dog_math.h"
#include "DogParamsCart2.h"

double GetSLCFL(double dt, const dTensorBC4& aux, const dTensorBC3& smax)
{

    const double prim_vol = dogParamsCart2.get_prim_vol();
    const int mx = smax.getsize(1);
    const int my = smax.getsize(2);

    double cfl2 = -100.0;
    //for (int j=1; j<=my; j++)
    for (int i=1; i<=mx; i++)
    {
        double tmp = smax.get(i,1,2);
        cfl2 = Max(dt*tmp/prim_vol, cfl2);
    }

    return cfl2;

//  double cfl = -100.0;
//  for (int j=1; j<=my; j++)
//  {
        //double tmp = Max(smax.get(i,j,1),smax.get(i,j,2));
//      double tmp = smax.get(1,j,1);
//      cfl = Max(dt*tmp/prim_vol, cfl);
//  }

    // double check we have the same cfl number ...
//  if( fabs( cfl - cfl2) > 1e-10 ) 
//  {

//      printf("Warning:  cfl numbers do not agree: \n");
//      printf("      cfl = %2.8e;  cfl2 = %2.8e;", cfl, cfl2);
//      printf("      fabs(cfl - cfl2) = %2.12e\n", fabs(cfl-cfl2) );
//      exit(1);
//  }

}
