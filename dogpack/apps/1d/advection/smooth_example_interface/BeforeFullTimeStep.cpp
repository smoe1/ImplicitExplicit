#include "DogSolver.h"
#include "tensors.h"
#include "DogParams.h"
#include <cmath>
using namespace std;

// Function that is called before a full time step
void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold, dTensorBC3& q)
{



  const int melems = q.getsize(1);
  const int meqn   = q.getsize(2);
  const int kmax   = q.getsize(3);
  const int mbc    = q.getmbc();
  const int maux   = aux.getsize(4);

  double time=dogParams.get_time();
  //double time=DogSolver::get_time_hack();
 
  const int mx     = melems;


      for(int i=1-mbc;i<=mx+mbc;i++)
      {
              for(int ma=1;ma<=maux;ma++)
              for(int k=1;k<=kmax;k++)
              {
                  //printf("i=%d, k=%d, Aux=%e \n",i,k,aux.get(i,ma,k));
              }
      }



  if(std::abs(time-0.75)<dt/2.0)
  {
      for(int i=1-mbc;i<=mx+mbc;i++)
      {
              for(int ma=1;ma<=maux;ma++)
              for(int k=1;k<=kmax;k++)
              {
                  //aux.set(i,ma,k,-aux.get(i,ma,k));
                  //auxold.set(i,ma,k,-auxold.get(i,ma,k));
              }
      }
  }



}
