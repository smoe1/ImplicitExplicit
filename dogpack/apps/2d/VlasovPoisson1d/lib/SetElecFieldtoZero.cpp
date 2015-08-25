#include <math.h>
#include "dogdefs.h"
#include "DogParams.h"       // for morder
#include "DogParamsCart2.h"
#include "VlasovParams.h"    // integration constant is stored here

void ComputeElecField(double t, const dTensor2& node1d, const dTensorBC4& qvals, dTensorBC3& Evals)
{
  const int  mbc = Evals.getmbc();
  const int   mx = Evals.getsize(1);
  const int meqn = Evals.getsize(2);
  const int kmax = Evals.getsize(3);

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  Evals.set(i,m,k, 0.0 );
	}
}


void ComputeElecPotential(double t, const dTensor2& node1d, const dTensorBC4& qvals, dTensorBC3& phi)
{
  const int  mbc = phi.getmbc();
  const int   mx = phi.getsize(1);
  const int meqn = phi.getsize(2);
  const int kmax = phi.getsize(3);

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  phi.set(i,m,k, 0.0 );
	}
}
