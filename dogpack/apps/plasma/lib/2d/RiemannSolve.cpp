#include "debug.h"
#include "assert.h"
#include "tensors.h"
#include "dog_math.h"
#include "PlasmaParams.h"
#include <cmath>
#include "RiemannSolve.h"

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc1
// * FluxFunc2
// * FluxFunc
// * SetWaveSpds

// by incorporating this method into a class the user can elect
// to avoid memory allocation on each call.
//
double RiemannSolver::solve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
  dTensor1& ffl   = fetch_ffl  ();
  dTensor1& ffr   = fetch_ffr  ();
  dTensor2& xedge_tmp = fetch_xedge();
  dTensor2& Ql_tmp    = fetch_Ql   ();
  dTensor2& Qr_tmp    = fetch_Qr   ();
  dTensor2& Auxl_tmp  = fetch_Auxl ();
  dTensor2& Auxr_tmp  = fetch_Auxr ();

  const int meqn = Ql.getsize();
  const int maux = Auxl.getsize();
  double speeds1[MAX_NUM_SUBSYSTEMS],speeds2[MAX_NUM_SUBSYSTEMS];

  // Reshape Ql,Qr,Auxl,Auxr,xedge to conform with FluxFunc.cpp
  for (int m=1; m<=meqn; m++)
  {  
      Ql_tmp.set(1,m, Ql.get(m) );  
      Qr_tmp.set(1,m, Qr.get(m) );
  }
  for (int m=1; m<=maux; m++)
  {  
      Auxl_tmp.set(1,m, Auxl.get(m) );  
      Auxr_tmp.set(1,m, Auxr.get(m) );
  }
  for (int m=1; m<=2; m++)
  {  xedge_tmp.set(1,m, xedge.get(m) );  }

  // Evaluate flux function at Ql and Qr
  const double n1 = nvec.get(1);
  const double n2 = nvec.get(2);
  assert_almost_eq(n1*n1+n2*n2,1.);
  // check for the cheap cases first
  if(n2==0.)
  {
    dTensor2& ftmp = fetch_flux1();

    void FluxFunc1(
      const dTensor2& xpts,
      const dTensor2& Q,
      const dTensor2& Aux,
      dTensor2& flux);
    FluxFunc1(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m) );  }

    FluxFunc1(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m) );  }
  }
  else if(n1==0.)
  {
    dTensor2& ftmp = fetch_flux1();

    void FluxFunc2(
      const dTensor2& xpts,
      const dTensor2& Q,
      const dTensor2& Aux,
      dTensor2& flux);
    FluxFunc2(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m) );  }

    FluxFunc2(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m) );  }
  }
  else
  {
    dTensor3& ftmp  = fetch_flux2 ();

    void FluxFunc(const dTensor2& xpts, const dTensor2& Q,
        const dTensor2& Aux, dTensor3& flux);
    FluxFunc(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m,1)*n1 + ftmp.get(1,m,2)*n2 );  }

    FluxFunc(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m,1)*n1 + ftmp.get(1,m,2)*n2 );  }
  }

  void SetWaveSpds(
    const dTensor1& nhat,
    const dTensor1& xedge,
    const dTensor1& Ql,
    const dTensor1& Qr,
    const dTensor1& Auxl,
    const dTensor1& Auxr,
    double* speeds1,double* speeds2);
  // Calculate minimum and maximum HLLE speeds
  SetWaveSpds(nvec,xedge,Ql,Qr,Auxl,Auxr,speeds1,speeds2);

  double smax_edge = 0.;
  // set fluxes for each subsystem
  for(int nsystem=1;
      nsystem<=plasmaParams.get_num_riemann_subsystems();nsystem++)
  {
    const double s1 = speeds1[nsystem];
    const double s2 = speeds2[nsystem];
    // dprintf("using nsystem=%d speeds s1=%f, s2=%f",nsystem,s1,s2);
    const int firstindex
      = plasmaParams.get_riemann_subsystem_firstindex(nsystem);
    const int lastindex
      = plasmaParams.get_riemann_subsystem_lastindex(nsystem);
    // Calculate Fluxes
    int mcase=0;
    if (fabs(s1)<=1.0e-12 && fabs(s2)<=1.0e-12)
    { mcase = 1; }
    else
    {
        if (s1*s2 > 0.0)
        {
            if (s1>0)
            { mcase = 1; }
            else
            { mcase = 2; }
        }
        else
        { mcase = 3; }
    }

    switch ( mcase )
    {
        default:
            invalid_value_error(mcase);
        case 1:  // both s1 and s2 are positive (or both zero)

            for (int m=firstindex; m<=lastindex; m++)
            { 
                Fl.set(m, ffl.get(m) );
                Fr.set(m, ffl.get(m) );
            }	
            break;

        case 2:  // both s1 and s2 are negative

            for (int m=firstindex; m<=lastindex; m++)
            { 
                Fl.set(m, ffr.get(m) );
                Fr.set(m, ffr.get(m) );
            }	
            break;

        case 3:  // s1 is negative and s2 is positive

            for (int m=firstindex; m<=lastindex; m++)
            {
                const double Qstar = (ffr.get(m)-ffl.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }	
            break;
    }

    smax_edge = Max(smax_edge, Max(fabs(s1),fabs(s2)));
  }
  return smax_edge;
}

// support the old RiemannSolve interface for backwards compatibility
//
double RiemannSolve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
  const int meqn = Ql.getsize();
  const int maux = Auxl.getsize();
  RiemannSolver riemannSolver(meqn,maux);
  return riemannSolver.solve(nvec, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);
}
