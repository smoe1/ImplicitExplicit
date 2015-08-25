#include "dogdefs.h"
#include "dog_math.h"
#include "L2Project.h"
#include "DogParamsCart3.h"
#include "InitialParams.h"
#include "L2ProjectInline3d.h"

void inline SetXYZ(const int i,
		   const int j,
		   const int k,
		   dTensor1& x,
		   dTensor1& y,
		   dTensor1& z)
{
  double   dx = dogParamsCart3.get_dx();
  double   dy = dogParamsCart3.get_dy();
  double   dz = dogParamsCart3.get_dz();

  double xlow = dogParamsCart3.get_xlow();
  double ylow = dogParamsCart3.get_ylow();
  double zlow = dogParamsCart3.get_zlow();

  double   xl = xlow + (double(i)-1.0)*dx;
  double   yl = ylow + (double(j)-1.0)*dy;
  double   zl = zlow + (double(k)-1.0)*dz;	  
  
  x.set(1, xl      );
  y.set(1, yl      );
  z.set(1, zl      );
	  
  x.set(2, xl + dx );
  y.set(2, yl      );
  z.set(2, zl      );
  
  x.set(3, xl      );
  y.set(3, yl + dy );
  z.set(3, zl      );
  
  x.set(4, xl      );
  y.set(4, yl      );
  z.set(4, zl + dz );

  x.set(5, xl + dx );
  y.set(5, yl + dy );
  z.set(5, zl      );
  
  x.set(6, xl + dx );
  y.set(6, yl      );
  z.set(6, zl + dz );
  
  x.set(7, xl      );
  y.set(7, yl + dy );
  z.set(7, zl + dz );
  
  x.set(8, xl + dx );
  y.set(8, yl + dy );
  z.set(8, zl + dz );
}

void inline LevelSetInitialCond(const int mlength,
				const dTensor1& x,
				const dTensor1& y,
				const dTensor1& z,
				iTensor1& zphi)
{
  for (int m=1; m<=mlength; m++)
    {
      const double rad = sqrt(pow(x.get(m),2)+pow(y.get(m),2)+pow(z.get(m),2));
      if (rad<0.3e0)
	{ zphi.set(m, 0 ); }
      else 
	{ zphi.set(m, 1 ); }
    }
}

void L2ProjectInitialCond(const int istart,
			  const int iend,
			  const int jstart,
			  const int jend,
			  const int kstart,
			  const int kend,
			  const int QuadOrder,
			  const int BasisOrder_qin,
			  const int BasisOrder_auxin,
			  const int BasisOrder_fout,    
			  const dTensorBC5* qin,
			  const dTensorBC5* auxin,
			  dTensorBC5* fout,
			  void (*Func)(const dTensor2&,
				       const dTensor2&,
				       const dTensor2&,
				       dTensor2&))
{
  // dx, dy, dz
  const double dx = dogParamsCart3.get_dx();
  const double dy = dogParamsCart3.get_dy();
  const double dz = dogParamsCart3.get_dz();

  // mbc
  const int mbc = qin->getmbc();
  assert_eq(mbc,auxin->getmbc());
  assert_eq(mbc,fout->getmbc());

  // qin variable
  const int       mx = qin->getsize(1);
  const int       my = qin->getsize(2);
  const int       mz = qin->getsize(3);
  const int     meqn = qin->getsize(4);
  const int kmax_qin = qin->getsize(5);
  int ktmp = (BasisOrder_qin*(BasisOrder_qin+2)*(BasisOrder_qin+1))/6;
  assert_eq(kmax_qin,ktmp);

  // auxin variable
  assert_eq(mx,auxin->getsize(1));
  assert_eq(my,auxin->getsize(2));
  assert_eq(mz,auxin->getsize(3));
  const int       maux = auxin->getsize(4);
  const int kmax_auxin = auxin->getsize(5);
  ktmp = (BasisOrder_auxin*(BasisOrder_auxin+2)*(BasisOrder_auxin+1))/6;
  assert_eq(kmax_auxin,ktmp);

  // fout variables
  //  TODO - why assume this has the same size?  What if we want to only
  //  project onto a single line, e.g. for padding boundary cell data? (-DS)
  assert_eq(mx,fout->getsize(1));
  assert_eq(my,fout->getsize(2));
  assert_eq(mz,fout->getsize(3));
  const int mcomps_out = fout->getsize(4);
  const int  kmax_fout = fout->getsize(5);
  assert_eq(kmax_fout,(BasisOrder_fout*(BasisOrder_fout+2)*(BasisOrder_fout+1))/6);

  // starting and ending indeces
  assert_ge(istart,1-mbc);
  assert_le(iend,mx+mbc);
  assert_ge(jstart,1-mbc);
  assert_le(jend,my+mbc);
  assert_ge(kstart,1-mbc);
  assert_le(kend,mz+mbc);

  // number of quadrature points
  assert_ge(QuadOrder,1);
  assert_le(QuadOrder,20);
  int QuadOrder_MOD = QuadOrder;
  switch(QuadOrder)
    {
    case 7:
      QuadOrder_MOD = 8;
      break;
    case 9:
      QuadOrder_MOD = 10;
      break;
    case 11:
      QuadOrder_MOD = 12;
      break;
    case 13:
      QuadOrder_MOD = 14;
      break;
    case 15:
      QuadOrder_MOD = 16;
      break;
    case 17:
      QuadOrder_MOD = 18;
      break;
    case 19:
      QuadOrder_MOD = 20;
      break;
    }
  const int mpoints = QuadOrder_MOD*QuadOrder_MOD*QuadOrder_MOD;

  // set quadrature point and weight information
  void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
  dTensor1 wgt(mpoints);
  dTensor2 spts(mpoints, 3);
  SetQuadWgtsPts(QuadOrder_MOD, wgt, spts);

  // Loop over each quadrature point to construct Legendre polys
  const int kmax = iMax(iMax(kmax_qin,kmax_auxin),kmax_fout);
  dTensor2 phi(mpoints, kmax);
  void SetLegendrePolys(const int, const int, const dTensor2&, dTensor2&);
  SetLegendrePolys(mpoints, kmax, spts, phi);

  // For efficiency compute weight*phi and then take transpose
  dTensor2 wgt_phi_transpose(kmax,mpoints);
  for(int mp=1;mp<=mpoints;mp++)
    for(int k=1;k<=kmax;k++)
      {
        wgt_phi_transpose.set(k,mp, wgt.get(mp)*phi.get(mp,k) );
      }

  // ------------------------------------------------------------- //
  // Loop over every grid cell indexed by user supplied parameters //
  // described by istart...iend, jstart...jend                     //
  // ------------------------------------------------------------- //
#pragma omp parallel for
  for (int i=istart; i<=iend; i++)
    for (int j=jstart; j<=jend; j++)
      for (int k=kstart; k<=kend; k++)
	{
	  dTensor1 x(8);
	  dTensor1 y(8);
	  dTensor1 z(8);
	  iTensor1 zphi(8);

	  SetXYZ(i,j,k,x,y,z);

	  LevelSetInitialCond(8,x,y,z,zphi);

	  int sumphi = zphi.get(1);
	  for (int m=2; m<=8; m++)
	    {  sumphi = sumphi + zphi.get(m); }
	  
	  switch(sumphi)
	    {
	    case 0:
	      fout->set(i,j,k,1,1, initialParams.rhol );    // density
	      fout->set(i,j,k,2,1, initialParams.m1l  );    // 1-momentum
	      fout->set(i,j,k,3,1, initialParams.m2l  );    // 2-momentum
	      fout->set(i,j,k,4,1, initialParams.m3l  );    // 3-momentum
	      fout->set(i,j,k,5,1, initialParams.El   );    // energy
	      
	      for (int me=1; me<=mcomps_out; me++)
		for (int mk=2; mk<=kmax_fout; mk++)
		  {  fout->set(i,j,k,me,mk, 0.0 );  }
	      break;

	    case 8:
	      fout->set(i,j,k,1,1, initialParams.rhor );    // density
	      fout->set(i,j,k,2,1, initialParams.m1r  );    // 1-momentum
	      fout->set(i,j,k,3,1, initialParams.m2r  );    // 2-momentum
	      fout->set(i,j,k,4,1, initialParams.m3r  );    // 3-momentum
	      fout->set(i,j,k,5,1, initialParams.Er   );    // energy
	      
	      for (int me=1; me<=mcomps_out; me++)
		for (int mk=2; mk<=kmax_fout; mk++)
		  {  fout->set(i,j,k,me,mk, 0.0 );  }
	      break;

	    default:
	      dTensor2   qvals(mpoints, meqn);
	      dTensor2 auxvals(mpoints, maux);
	      dTensor2    xpts(mpoints, 3);
	      
	      // Flux function and its Jacobian:
	      dTensor2   fvals(mpoints, mcomps_out);
	      
	      //find center of current cell
	      const double xc = dogParamsCart3.get_xc(i);
	      const double yc = dogParamsCart3.get_yc(j);
              const double zc = dogParamsCart3.get_zc(k);

              // Compute q, aux and fvals at each Gaussian quadrature point
              // for this current cell indexed by (i,j,k)
              // Save results into dTensor2 qvals, auxvals and fvals.
              L2ProjectInline3d::set_vals_at_each_Gaussian_quadrature_point(i, j, k,
                                                                            mpoints, meqn, maux,
                                                                            kmax_qin, kmax_auxin,
                                                                            xc, yc, zc, dx, dy, dz,
                                                                            spts, phi, qin, auxin,
                                                                            xpts, qvals, auxvals);
              // Evaluate Func at Gaussian quadrature point
              Func(xpts, qvals, auxvals, fvals);

              // Evaluate integral on current cell (project onto Legendre basis)
              // using Gaussian Quadrature for the integration
              L2ProjectInline3d::integrate_on_current_cell(false, i, j, k,
                                                           mcomps_out, kmax_fout, mpoints,
                                                           wgt_phi_transpose,
                                                           fvals, fout);
	    }
	}
}
