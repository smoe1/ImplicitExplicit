#include <cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"

// Non-conservative advection equation:
//
//       q_t + u(x,y) q_x + v(x,y) q_y = 0
//
// DG weak formulation
//
//              /                   //
//      Q_t = - | \phi F . n ds  +  ||  grad(\phi) . F dx dy  
//              /                   //
//
//              //
//           +  || \phi q div(u) dx dy  
//              //
//
//    NOTE:  The first two integrals have already been taken 
//           care of by the standard "ConstructL" code. The
//           last term (i.e., the non-conservative term) still
//           needs to be approximated. That is the job of this
//           function.
//
void LstarExtra(const dTensorBC4& q, 
		const dTensorBC4& aux, 
		dTensorBC4& Lstar)
{
  int i,j,k;
  int mx   = q.getsize(1);
  int my   = q.getsize(2);
  int meqn = q.getsize(3);
  int kmax = q.getsize(4);
  int mbc  = q.getmbc();
  int maux = aux.getsize(3);        
  int morder = int((sqrt(1+8*kmax)-1)/2);
  const double dx = dogParamsCart2.get_dx();
  const double dy = dogParamsCart2.get_dy();
  void L2Project(const int istart, 
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
			      const dTensor2&,dTensor2&));
  void DivFunc(const dTensor2& xpts, const dTensor2& q, 
               const dTensor2& divU, dTensor2& outvals);
  dTensorBC4   divU(mx,my,1,kmax,mbc);
  dTensorBC4 Lextra(mx,my,meqn,kmax,mbc);
  dTensor1 U(15),V(15);

  // Initialize U and V
  for (k=1; k<=15; k++)
    {
      U.set(k, 0.0 );
      V.set(k, 0.0 );
    }

  // Compute div(U) in each element
  switch(morder)
    {

    case 5:
      
      for (i=(1-mbc); i<=(mx+mbc); i++)
        for (j=(1-mbc); j<=(my+mbc); j++)
          {
            for (k=1; k<=kmax; k++)
              {
                U.set(k, aux.get(i,j,1,k) );
                V.set(k, aux.get(i,j,2,k) );
              }

            divU.set(i,j,1,1,  (2.0*(U.get(2)*sq3+U.get(9)*sq7))/dx 
                     + (2.0*(V.get(3)*sq3+V.get(10)*sq7))/dy );
            divU.set(i,j,1,2,  2.0*sq3*(U.get(5)*sq5+3*U.get(14))/dx 
                     + (2.0*(V.get(12)*sq7+sq3*V.get(4)))/dy );
            divU.set(i,j,1,3,  (2.0*(U.get(11)*sq7+sq3*U.get(4)))/dx 
                     + 2.0*sq3*(3*V.get(15)+V.get(6)*sq5)/dy );
            divU.set(i,j,1,4,  2.0*U.get(7)*sq3*sq5/dx + 2.0*sq3*V.get(8)*sq5/dy );
            divU.set(i,j,1,5,  2.0*sq5*U.get(9)*sq7/dx + 2.0*V.get(7)*sq3/dy );
            divU.set(i,j,1,6,  2.0*sq3*U.get(8)/dx + 2.0*sq5*V.get(10)*sq7/dy );
            divU.set(i,j,1,7,  2.0*sq5*sq7*U.get(11)/dx + 2.0*sq5*V.get(13)*sq3/dy );
            divU.set(i,j,1,8,  2.0*sq5*U.get(13)*sq3/dx + 2.0*sq5*sq7*V.get(12)/dy );
            divU.set(i,j,1,9,  6.0*U.get(14)*sq7/dx + 2.0*sq3*V.get(11)/dy );
            divU.set(i,j,1,10, 2.0*sq3*U.get(12)/dx + 6.0*V.get(15)*sq7/dy );
            divU.set(i,j,1,11, 0.0 );
            divU.set(i,j,1,12, 0.0 );
            divU.set(i,j,1,13, 0.0 );
            divU.set(i,j,1,14, 0.0 );
            divU.set(i,j,1,15, 0.0 );
          }
      break;

    case 4:
      
      for (i=(1-mbc); i<=(mx+mbc); i++)
        for (j=(1-mbc); j<=(my+mbc); j++)
          {
            for (k=1; k<=kmax; k++)
              {
                U.set(k, aux.get(i,j,1,k) );
                V.set(k, aux.get(i,j,2,k) );
              }

            divU.set(i,j,1,1,  (2.0*(U.get(2)*sq3+U.get(9)*sq7))/dx 
                     + (2.0*(V.get(3)*sq3+V.get(10)*sq7))/dy );
            divU.set(i,j,1,2,  2.0*sq3*(U.get(5)*sq5+3*U.get(14))/dx 
                     + (2.0*(V.get(12)*sq7+sq3*V.get(4)))/dy );
            divU.set(i,j,1,3,  (2.0*(U.get(11)*sq7+sq3*U.get(4)))/dx 
                     + 2.0*sq3*(3*V.get(15)+V.get(6)*sq5)/dy );
            divU.set(i,j,1,4,  2.0*U.get(7)*sq3*sq5/dx + 2.0*sq3*V.get(8)*sq5/dy );
            divU.set(i,j,1,5,  2.0*sq5*U.get(9)*sq7/dx + 2.0*V.get(7)*sq3/dy );
            divU.set(i,j,1,6,  2.0*sq3*U.get(8)/dx + 2.0*sq5*V.get(10)*sq7/dy );
            divU.set(i,j,1,7,  0.0 );
            divU.set(i,j,1,8,  0.0 );
            divU.set(i,j,1,9,  0.0 );
            divU.set(i,j,1,10, 0.0 );
          }
      break;

    case 3:
           
      for (i=(1-mbc); i<=(mx+mbc); i++)
        for (j=(1-mbc); j<=(my+mbc); j++)
          {
            for (k=1; k<=kmax; k++)
              {
                U.set(k, aux.get(i,j,1,k) );
                V.set(k, aux.get(i,j,2,k) );
              }

            divU.set(i,j,1,1,  2.0*sq3*(U.get(2)/dx + V.get(3)/dy) );
            divU.set(i,j,1,2,  2.0*sq3*(sq5*U.get(5)/dx + V.get(4)/dy) );
            divU.set(i,j,1,3,  2.0*sq3*(U.get(4)/dx + sq5*V.get(6)/dy) );     
            divU.set(i,j,1,4,  0.0 );
            divU.set(i,j,1,5,  0.0 );
            divU.set(i,j,1,6,  0.0 );
          }
      break;

    case 2:
      
      for (i=(1-mbc); i<=(mx+mbc); i++)
        for (j=(1-mbc); j<=(my+mbc); j++)
          {
            for (k=1; k<=kmax; k++)
              {
                U.set(k, aux.get(i,j,1,k) );
                V.set(k, aux.get(i,j,2,k) );
              }

            divU.set(i,j,1,1, 2.0*sq3*(U.get(2)/dx + V.get(3)/dy) );     
            divU.set(i,j,1,2,  0.0 );
            divU.set(i,j,1,3,  0.0 );
          }
      break;

    case 1:

      for (i=(1-mbc); i<=(mx+mbc); i++)
        for (j=(1-mbc); j<=(my+mbc); j++)
          {
            divU.set(i,j,1,1,  0.0e0 );
          }
      break;
    }

  // Project  q * div(U) onto the Legendre basis
  L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
	    morder,morder,morder,morder,
	    &q,&divU,&Lextra,&DivFunc);

  // update solution
  for (i=(2-mbc); i<=(mx+mbc-1); i++)
    for (j=(2-mbc); j<=(my+mbc-1); j++)
      for (k=1; k<=kmax; k++)
        {
	  Lstar.set(i,j,1,k, Lstar.get(i,j,1,k) + Lextra.get(i,j,1,k) );
        }
  
}

// Function to evaluate q * div(U)
void DivFunc(const dTensor2& xpts, const dTensor2& q, 
             const dTensor2& divU, dTensor2& outvals)
{
  int i,m;
  int numpts=xpts.getsize(1);
  double x,y,qtmp,divUtmp;
  
  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i,1);
      y = xpts.get(i,2);
      
      qtmp = q.get(i,1);
      divUtmp = divU.get(i,1);
      
      outvals.set(i,1, qtmp*divUtmp );
    }
}

