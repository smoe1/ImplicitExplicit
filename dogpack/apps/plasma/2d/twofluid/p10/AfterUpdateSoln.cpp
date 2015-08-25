#if 0
#include "dogdefs.h"
#include "PlasmaParams.h"
#include "Components.h"
#include "Relax.h"
#include "DogParams.h"

void relax(
   const dTensor2& xpts,
   const dTensor2& qvals,
   const dTensor2& auxvals,
   dTensor2& q_new,
   void* data)
{
    int numpts=xpts.getsize(1);
    double d_time = *(double*)data;

    int top_copy_idx;
    bool did_ion_relax = Relax::isotropize_spc(
        q_new,
        qvals,
        numpts,
        d_time,
        plasmaParams.get_iso_period_type(),
        plasmaParams.get_ion_base_iso_period(),
        plasmaParams.get_ion_mass(),
        _rho_i,
        _M1_i,
        _M2_i,
        _M3_i,
        _N11_i,
        _N12_i,
        _N13_i,
        _N22_i,
        _N23_i,
        _N33_i,
        _B1,
        _B2,
        0);
    top_copy_idx = 10;
    if(did_ion_relax) top_copy_idx = 4;
    for (int i=1; i<=numpts; i++)
    {
        for(int k=1;k<=top_copy_idx;k++) q_new.set(i,k, qvals.get(i,k));
    }
    for (int i=1; i<=numpts; i++)
    {
        for(int k=11;k<=q_new.getsize(2);k++) q_new.set(i,k, qvals.get(i,k));
    }
}

// This is done *before* limiters are applied.
void AfterUpdateSoln(
  const dTensorBC4& aux,
  dTensorBC4& q,
  double dt,
  double beta)
{
    if(plasmaParams.get_relaxation_type()!=RelaxationType::INSTANTANEOUS)
      return;
    int mxelems = q.getsize(1);
    int myelems = q.getsize(2);
    int    mbc = q.getmbc();

    void L2Project_extra(const int istart, 
			 const int iend, 
			 const int jstart, 
			 const int jend,
			 const int QuadOrder, 
			 const int BasisOrder_qin,
			 const int BasisOrder_auxin,
			 const int BasisOrder_fout,		  
			 const dTensorBC4* qin,
			 const dTensorBC4* auxin, 
			 dTensorBC4* fout,
			 void (*Func)(const dTensor2&,const dTensor2&,
				      const dTensor2&,dTensor2&,void* data),
			 void* data);
    double d_time = dt*beta;
    const int space_order = dogParams.get_space_order();
    L2Project_extra(2-mbc,mxelems+mbc-1,2-mbc,myelems+mbc-1,
		    space_order,space_order,space_order,space_order,
		    &q,&aux,&q,&relax,&d_time);
}
#endif
