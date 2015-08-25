#include "dogdefs.h"

// add generic weak term diffusive term to the right-hand side:
//   q_t = Div(g(q,grad(h(q))))
//
// This method is a generic framework for weak
// diffusive terms and should be moved into lib/2d/cart
// 
// In a typical application h(q) might be the temperature
// and g might be the heat flux (diffusive thermal energy flux).
// (We allow g to depend also on q to allow the
// the heat conductivity to depend on the state.)
//
// In the five-moment case viscosity is another such term.
// In this framework h(q)=u (i.e. the velocity)
// and g(q,grad h) is the stress tensor
// (diffusive momentum flux).
//
// The user is expected to supply callbacks for h and g.
//
// This framework assumes that the diffusive flux g
// is small.  I think it also should depend weakly
// on u in comparison to grad h, since if there is no
// dependence on grad h then this could be a hyperbolic term.
//
// To handle this term we make the definition
//
//   w = grad h,
//
// multiply the equations
//
//   q_t = Div(g(u,w))
//   w = grad h,
//
// by ((reciprocal) basis) test functions phi,
// and integrate by parts:
//
//   D_t(int(q*phi)) = oint(n.g(u,w)*phi) - int(g(u,w).grad phi),
//   int(w*psi) = oint(n*h(u)*psi) - int(h(u)*grad(phi))
//
void LstarWeakDiffusion(
  dTensorBC4& q, dTensorBC4& aux, dTensorBC4& Lstar,
  /*g=*/fluxFunc,
  /*h=*/stateFunc,
  // number of components of h
  int numeq)
{
  const int mx   = q->getsize(1);
  const int my   = q->getsize(2);
  const int kmax = q->getsize(4);
  const int mbc  = q->getmbc();
  // need to know number of components of h
  // in order to allocate the appropriate amount
  // of space for v.
  dTensorBC4(mx,my,numeq,kmax,mbc) vx;
  dTensorBC4(mx,my,numeq,kmax,mbc) vy;
  //
  // int(v*psi) = oint(n*h(u)*psi) - int(h(u)*grad(phi))
  //
  // for edge integrals we take the average from each side.
  //
  // We need to use the same flux for the cells on both
  // sides of a given boundary, so for consistency and
  // computational efficiency we traverse over the
  // boundaries and add or subtract the boundary integral
  // from the cell on each side.
  for (int i=(2-mbc); i<=(mx+mbc); i++)
  for (int j=(2-mbc); j<=(my+mbc-1); j++)
  for (int k=1; k<=space_order; k++)
  // [why does ConstructL look at space_order?]
  // because it is the number of quadrature points
}

void LstarWeakDiffusionTerm(q,aux,Lstar,fluxFunc,stateFunc)
{
  // allocate space for 
}

void ComputeGradients(q,aux,stateFunc,wxx,wxy,wyx,wyy)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
  // add contribution of interior integrals
  //
  // L2ProjectGrad has the desired functionality
  // but (1) it forms an unwanted dot product
  // and (2) it computes the entire state at every
  // sample point, which is inefficient for
  // a large, loosely coupled system if we only
  // wish to modify a few specific components.
  //
  // Unfortunately the representation of cells
  // is not modularized, making it difficult
  // to pass a cell's state to a generic routine
  // in a way that allows the routine to access
  // only the state of that particular cell.
  // We could create a Cart2Cell class which
  // adopts the cell information from an underlying
  // dTensor4 array.
  //
  // It will suffice if we can specify a *range*
  // of state indices.

  // In this algorithm keep in mind that we need
  // to evolve one layer of ghost cells (so that
  // we can limit the boundary cells).

  // for the interior of each cell:
  //
  // evaluate h at each interior quadrature point.
  // integrate against phi,x and subtract from wxx and wxy
  // integrate against phi,y and subtract from wyy and wyx
  {
    // iterate over cells (including one layer of boundary cells)
    # pragma omp parallel for
    for(int i=0;i<=mx+1;i++)
    for(int j=0;j<=my+1;j++)
    {
    }
  }

  // for each boundary perpendicular to x axis:
  //
  // evaluate h on high side, compute edge integral,
  // and add full value to wxx
  // and half the value to wyx
  //
  // if wxy is not null:
  // evaluate h on low side, compute edge integral,
  // and add half the value to wyx

  // for each boundary perpendicular to y axis:
  //
  // evaluate h on high side, compute edge integral,
  // and add full value to wyy
  // and half the value to wxy
  //
  // if wyx is not null:
  // evaluate h on low side, compute edge integral,
  // and add half the value to wxy

  // add contribution of boundary integrals
  wxx.set(i,j,m,k, )
}

void LstarVstress(const dTensorBC4& q, const dTensorBC4& aux, dTensorBC4& Lstar)
{
  // get w := Grad u, where u = fluid velocity
  //
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
  dTensorBC4 wxx(mx,my,3,kmax,mbc);
  dTensorBC4 wyx(mx,my,3,kmax,mbc); // 1/3 waste for simplicity
  dTensorBC4 wyx(mx,my,3,kmax,mbc); // 1/3 waste for simplicity
  dTensorBC4 wyy(mx,my,3,kmax,mbc);
  //
  ComputeGradients(q,aux,velocityFunc,wxx,wxy,wyx,wyy);
  ComputeDiv(q,aux,wxx,wxy,wyx,wyy,vStressFunc,Lstar);
}

void LstarHeat(const dTensorBC4& q, const dTensorBC4& aux, dTensorBC4& Lstar)
{
  // get w:= Grad T, where T = fluid temperature
  //
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
  dTensorBC4 wxx(mx,my,1,kmax,mbc);
  dTensorBC4 null(mx,my,0,kmax,mbc);
  dTensorBC4 wyy(mx,my,1,kmax,mbc);
  //
  ComputeGradients(q,aux,temperatureFunc,wxx,null,null,wyy);
  ComputeDiv(q,au,wxx,null,null,wyy,heatFluxFunc,Lstar);
}

void LstarExtra(dTensorBC4& q, dTensorBC4& aux, dTensorBC4& Lstar)
{
  LstarVstress(q,aux,Lstar);
  LstarHeat(q,aux,Lstar);
  //LstarWeakDiffusionTerm(q,aux,Lstar,stressFunc,velFunc);
  //LstarWeakDiffusionTerm(q,aux,Lstar,heatFluxFunc,tempFunc);
}
