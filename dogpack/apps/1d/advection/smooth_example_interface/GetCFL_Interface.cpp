#include <iostream>
#include "dog_math.h"
#include "tensors.h"
#include "stdlib.h"
#include "stdio.h"
using namespace std;

double GetCFL_Interface(double dt,double dtmax,
        const dTensor1& prim_vol,
        const int method[],
        const dTensorBC3& aux,
        const dTensorBC1& smax,const dTensor2& dxi,const dTensorBC1 global2interf)
{
    const int melems =aux.getsize(1);
    double cfl=-100.0;

    if (dt>dtmax)
    {
        cout << endl;
        cout << " Error: dt is out of bounds ... " << endl;
        cout << "     dt = " << dt << endl;
        cout << endl;
        exit(1);
    }

// This can't be done in parallel!!!  (-DS)
//#pragma omp parallel for
    for (int j=1; j<=melems; j++)
    {
      if(abs(global2interf.get(j))<=1.0e-1)
      {
        cfl = Max(dt*smax.get(j)/prim_vol.get(j),cfl);  
      }
      else
      {
        int iint=int(global2interf.get(j));
        cfl = Max(dt*smax.get(j)/prim_vol.get(j),cfl);  
        //cfl = Max(Max(dt*smax.get(j)/dxi.get(iint,1),dt*smax.get(j)/dxi.get(iint,2)),cfl);
      }

    }

    if (cfl>1.0e8)
    {
        cout << endl;
        cout << " Error: CFL number is out of bounds ... " << endl;
        cout << "     CFL = " << cfl << endl;
        cout << endl;
        exit(1);
    }

    return cfl;
}
