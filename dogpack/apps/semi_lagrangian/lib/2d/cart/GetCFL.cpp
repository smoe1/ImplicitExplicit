#include "tensors.h"
#include "dog_math.h"
#include "DogParamsCart2.h"

// we should parallelize this. -eaj
double GetCFL(double dt, const dTensorBC4& aux, const dTensorBC3& smax)
{
    int i,j,k,mcapa;
    double capa,tmp;
    double cfl=-100.0;
    int mx = aux.getsize(1);
    int my = aux.getsize(2);
   
    const double prim_vol = dogParamsCart2.get_prim_vol();
    for (j=1; j<=my; j++)
      for (i=1; i<=mx; i++)
        {
	  tmp = Max(smax.get(i,j,1),smax.get(i,j,2));
	  cfl = Max(tmp*(dt/prim_vol), cfl);
	}

    return cfl;
}
