#define MAX_DEBUG_LEVEL 3
#include "gasCart2.h"
#include "dogdefs.h"
#include "L2Project.h"
#include "DogParams.h"
//#include "DogSolverCart2.h"

void UpdateEntropies(
    double dt,dTensorBC4& aux,
    dTensorBC4& q,
    int firstEntropyIdx,
    int numEntropyIndices,
    void (*reset_entropy)(
       const dTensor2& xpts,
       const dTensor2& qvals,
       const dTensor2& auxvals,
       dTensor2& q_new)
    )
{
    if(q.getsize(3)<firstEntropyIdx) return;

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);
    assert_eq(maux,numEntropyIndices);

    // record the rate of entropy production
    // (problem: in the case of thermal exchange between species
    // the total entropy production hereby calculated seems
    // difficult to interpret in a way that makes it higher than
    // first-order-accurate, since each species is advected with
    // its own velocity)

    // compute the entropy of each species,
    // reset the entropy, and store the entropy
    // production in aux.
    
    // store the entropy
    int idxOffset = firstEntropyIdx-1;
    dTensor4 entropies(mx,my,numEntropyIndices,kmax);
    #pragma omp parallel for
    for(int i=1;i<=mx;i++)
    for(int j=1;j<=my;j++)
    for(int k=1;k<=kmax;k++)
    for(int m=1;m<=numEntropyIndices;m++)
    {
      entropies.set(i,j,m,k, q.get(i,j,m+idxOffset,k));
    }
    // reset the entropy
    const int space_order = dogParams.get_space_order();
    L2Project(2-mbc,mx+mbc-1,2-mbc,my+mbc-1,
	      q,aux,q,reset_entropy);
	      //space_order,space_order,space_order,space_order,
	      //&q,&aux,&q,reset_entropy);
    // save the difference
  //eprintf("numEntropyIndices=%d",numEntropyIndices);
    #pragma omp parallel for
    for(int i=1;i<=mx;i++)
    for(int j=1;j<=my;j++)
    for(int m=1;m<=numEntropyIndices;m++)
    for(int k=1;k<=kmax;k++)
    {
      double entropy_production_rate
        = (q.get(i,j,m+idxOffset,k)-entropies.get(i,j,m,k))/dt;
      aux.set(i,j,m,k,  entropy_production_rate);
      {
        dprintf3("entropies.get(%d,%d, %d,%d) = %g",
          i,j,m,k,entropies.get(i,j,m,k));
        dprintf3("        q.get(%d,%d,%d,%d) = %g",
          i,j,m+idxOffset,k, q.get(i,j,m+idxOffset,k));
        dprintf3("entropy_production_rate = %g",entropy_production_rate);
      }
    }
}

// For passive convection simply copy
//
void ProjectRightConvectedScalar(int idx,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
    if(Q.getsize(1)<idx) return;

    assert_eq(Q.getsize(2),W.getsize(2));
    for(int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(idx,k, W.get(idx,k) );
    }
}
void ProjectLeftConvectedScalar(int idx,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
    if(Q.getsize(1)<idx) return;
    
    assert_eq(Q.getsize(2),W.getsize(2));
    for(int k=1; k<=Q.getsize(2); k++)
    {
        W.set(idx,k, Q.get(idx,k) );
    }
}
