#include "dogdefs.h"
#include "fcmp.h"
#include "diffusionCart2.h"
#include "Components.h"
#include "PlasmaParams.h"
#include "DogParamsCart2.h"
#include "DogParams.h"

const int fxy_map[7] = {
  0, 0, 1, 0, 2, 3, 0
//0, 1, 2, 3, 4, 5, 6
};

const int fyx_map[7] = {
  0, 1, 2, 3, 0, 0, 0
//0, 1, 2, 3, 4, 5, 6
};

bool do_checks = true;

// components of heat flux flux
//
enum HeatFFid{
    _fxx11=1,
    _fxx12=2,
    _fxx13=3,
    _fxx22=4,
    _fxx23=5,
    _fxx33=6,

    _fxy12=7,
    _fxy22=8,
    _fxy23=9,

    _fyx11=10,
    _fyx12=11,
    _fyx13=12,

    _fyy11=13,
    _fyy12=14,
    _fyy13=15,
    _fyy22=16,
    _fyy23=17,
    _fyy33=18,
};
const int fxy_offset = _fxy12-1;
const int fyx_offset = _fyx11-1;
const int fyy_offset=_fyy11-1;
const int mheatff=18;

#if 0
#define assert_equal(a,b) assert_equal_func(a,b,__FILE__,__LINE__)

bool assert_equal_func(double a, double b, char* file, int line)
{
  double size = abs(a)+abs(b);
  double diff = abs(a-b);
  if(size < 1)
  {
    if(diff<EPSILON) return true;
  }
  else
  {
    double ratio = abs(b-a)/size;
    if(ratio<EPSILON) return true;
  }
  eprintf("assertion failed in file %s line %d:\n"
           " %f==%f\n", file, line,a,b);
  return false;
}
#endif

// Calculate f = heatff (the "heat flux fluxes").
//
// Q111 = fxxi11,x
// Q112 = fxxi12,x + fxyi12,y
// Q113 = fxxi13,x
// Q122 = fxxi22,x + fxyi22,y
// Q123 = fxxi23,x + fxyi23,y
// Q133 = fxxi33,x
//
// Q211 = fyxi11,x + fyyi11,y
// Q212 = fyxi12,x + fyyi12,y
// Q213 = fyxi13,x + fyyi13,y
// Q222 =          + fyyi22,y
// Q223 =          + fyyi23,y
// Q233 =          + fyyi33,y
//
// fxxi11 = 3 T11*\nu
// fxxi12 = 2 T12*\nu, fxyi12 =   T11*\nu
// fxxi13 = 2 T13*\nu
// fxxi22 =   T22*\nu, fxyi22 = 2 T12*\nu
// fxxi23 =   T23*\nu, fxyi23 =   T13*\nu
// fxxi33 =   T33*\nu
//
// fyxi11 = 2 T12*\nu, fyyi11 =   T11*\nu
// fyxi12 =   T22*\nu, fyyi12 = 2 T12*\nu
// fyxi13 =   T23*\nu, fyyi13 =   T13*\nu
//                     fyyi22 = 3 T22*\nu 
//                     fyyi23 = 2 T23*\nu 
//                     fyyi33 =   T33*\nu 
//
void heatffFunc(const dTensor2& xpts, const dTensor2& q,
             const dTensor2& aux, dTensor2& heatff)
{
    int numpts=q.getsize(1);
    //int meqn=heatff.getsize(2);
        
    const double ion_mass = plasmaParams.get_ion_mass();
    const double nu = plasmaParams.get_ion_heat_conductivity();
    for (int i=1; i<=numpts; i++)
    {
        //double x = xpts.get(i,1);
        //double y = xpts.get(i,2);

        double T11, T12, T13, T22, T23, T33;
        // get the temperature tensor
        {
            double rho = q.get(i,_rho_i);
            double M1 = q.get(i,_M1_i);
            double M2 = q.get(i,_M2_i);
            double M3 = q.get(i,_M3_i);
            double u1 = M1/rho;
            double u2 = M2/rho;
            double u3 = M3/rho;
            assert_printf( rho==rho && M1==M1 && M2==M2 && M3==M3,
              "failure for i=%d\n"
              "rho=%f, M1=%f, M2=%f, M3=%f\n",
              i,rho,M1,M2,M3);

            // get the pressure tensor
            double P11 = q.get(i,_N11_i) - M1*u1;
            double P12 = q.get(i,_N12_i) - M1*u2;
            double P13 = q.get(i,_N13_i) - M1*u3;
            double P22 = q.get(i,_N22_i) - M2*u2;
            double P23 = q.get(i,_N23_i) - M2*u3;
            double P33 = q.get(i,_N33_i) - M3*u3;
            if(do_checks)
            {
              assert_isnum(P11);
              assert_isnum(P12);
              assert_isnum(P13);
              assert_isnum(P22);
              assert_isnum(P23);
              assert_isnum(P33);

              assert_gt(rho,EPSILON);
            }

            // get the number density
            double n = rho/ion_mass;

            // divide the pressure by the number density
            T11 = P11/n;
            T12 = P12/n;
            T13 = P13/n;
            T22 = P22/n;
            T23 = P23/n;
            T33 = P33/n;
        }

        double TA11, TA12, TA13, TA22, TA23, TA33;
        double detT;
        // find adjugate and det of temperature tensor
        {
            // find the adjugate
            TA11 = T22*T33 - T23*T23;
            TA12 = T23*T13 - T12*T33;
            TA13 = T12*T23 - T22*T13;
            TA22 = T11*T33 - T13*T13;
            TA23 = T12*T13 - T11*T23;
            TA33 = T11*T22 - T12*T12;
            // find the determinant
            detT = TA11*T11 + TA12*T12 + TA13*T13;
            // verify that T is positive definite
            assert_gt(T11, EPSILON);
            assert_gt(TA33, EPSILON);
            assert_gt(detT, EPSILON);
            if(false) // verified on data
            // verify that TA*T = id*det(T)
            {
                // dot product of corresponding rows is determinant
                //bool test_equal(double a, double b);
                assert_almost_eq(detT,TA11*T11+TA12*T12+TA13*T13);
                assert_almost_eq(detT,TA12*T12+TA22*T22+TA23*T23);
                assert_almost_eq(detT,TA13*T13+TA23*T23+TA33*T33);
                // first column times other rows is zero
                assert_almost_eq(0.,  TA12*T11+TA22*T12+TA23*T13);
                assert_almost_eq(0.,  TA13*T11+TA23*T12+TA33*T13);
                // second column times other rows is zero
                assert_almost_eq(0.,  TA11*T12+TA12*T22+TA13*T23);
                assert_almost_eq(0.,  TA13*T12+TA23*T22+TA33*T23);
                // third column times other rows is zero
                assert_almost_eq(0.,  TA11*T13+TA12*T23+TA13*T33);
                assert_almost_eq(0.,  TA12*T13+TA22*T23+TA23*T33);
            }
        }

        // compute the heat flux flux

        // compute rescaled negated adjugate (TA*nu/detT = invT*nu):
        //const double m_nu_detT = nu/detT;
        const double m_nu_detT = -nu/detT;
        const double TAR11 = m_nu_detT*TA11;
        const double TAR12 = m_nu_detT*TA12;
        const double TAR13 = m_nu_detT*TA13;
        const double TAR22 = m_nu_detT*TA22;
        const double TAR23 = m_nu_detT*TA23;
        const double TAR33 = m_nu_detT*TA33;

        heatff.set(i,_fxx11, 3*TAR11);
        heatff.set(i,_fxx12, 2*TAR12);
        heatff.set(i,_fxx13, 2*TAR13);
        heatff.set(i,_fxx22,   TAR22);
        heatff.set(i,_fxx23,   TAR23);
        heatff.set(i,_fxx33,   TAR33);

        heatff.set(i,_fxy12,   TAR11);
        heatff.set(i,_fxy22, 2*TAR12);
        heatff.set(i,_fxy23,   TAR13);

        heatff.set(i,_fyx11, 2*TAR12);
        heatff.set(i,_fyx12,   TAR22);
        heatff.set(i,_fyx13,   TAR23);

        heatff.set(i,_fyy11,   TAR11);
        heatff.set(i,_fyy12, 2*TAR12);
        heatff.set(i,_fyy13,   TAR13);
        heatff.set(i,_fyy22, 3*TAR22);
        heatff.set(i,_fyy23, 2*TAR23);
        heatff.set(i,_fyy33,   TAR33);
    }
}

void LstarExtra_helper(double dx, double dy,
  const dTensorBC4& f, dTensorBC4& Lstar)
{
  const int mx   = f.getsize(1);
  const int my   = f.getsize(2);
  const int kmax = f.getsize(4);
  const int mbc  = f.getmbc();
  dTensorBC4 Vxx(mx,my,6,kmax,mbc);
  dTensorBC4 Vyy(mx,my,6,kmax,mbc);
  // allocate for 3 components to avoid wasting RAM
  dTensorBC4 Vxy(mx,my,3,kmax,mbc);
  dTensorBC4 Vyx(mx,my,3,kmax,mbc);
  
  // need to properly curtail boundary index references;
  // rework the loops and split them up.
  {
  // compute Vxx and Vyy
  #pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc-1); i++)
  for (int j=(1-mbc); j<=(my+mbc-1); j++)
  for (int k=1; k<=kmax; k++)
  {
    for (int m=1; m<=6; m++)
    {
      double BxpmFxx = 0; double AxmFxx = 0;
      double BypmFyy = 0; double AymFyy = 0;
      for (int n=1; n<=kmax; n++)
      {
        if(Bx_nonzero[k][n]) // use sparsity to accelerate
        {
          BxpmFxx    += Bxpm   [k][n]*f.get(i+1,j  ,m,n);
          AxmFxx     += Axm    [k][n]*f.get(i  ,j  ,m,n);
        }

        if(By_nonzero[k][n]) // use sparsity to accelerate
        {
          const int myy_idx = fyy_offset+m;
          BypmFyy    += Bypm   [k][n]*f.get(i  ,j+1,myy_idx,n);
          AymFyy     += Aym    [k][n]*f.get(i  ,j  ,myy_idx,n);
        }
      }
      Vxx.set(i,j,m,k, (BxpmFxx-AxmFxx)/dx);
      Vyy.set(i,j,m,k, (BypmFyy-AymFyy)/dy);
    }
  }
  // compute Vyx
  #pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc-1); i++)
  for (int j=(1-mbc); j<=(my+mbc); j++)
  for (int k=1; k<=kmax; k++)
  {
    for (int m=1; m<=6; m++)
    if(fyx_map[m])
    {
      double half_CxFyx = 0; double BxpmFyx = 0; double BxmpFyx = 0;
      for (int n=1; n<=kmax; n++)
      if(Bx_nonzero[k][n]) // use sparsity to accelerate
      {
            const int myx_idx = fyx_offset+fyx_map[m];
            if(half_Cx[k][n])
              half_CxFyx += half_Cx[k][n]*f.get(i  ,j  ,myx_idx,n);
            BxpmFyx      += Bxpm   [k][n]*f.get(i+1,j  ,myx_idx,n);
            BxmpFyx      += Bxmp   [k][n]*f.get(i-1,j  ,myx_idx,n);
      }
      Vyx.set(i,j,fyx_map[m],k, (half_CxFyx+0.5*(BxpmFyx-BxmpFyx))/dx);
    }
  }
  // compute Vxy
  #pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
  for (int j=(2-mbc); j<=(my+mbc-1); j++)
  for (int k=1; k<=kmax; k++)
  {
    for (int m=1; m<=6; m++)
    if(fxy_map[m])
    {
      double half_CyFxy = 0; double BypmFxy = 0; double BympFxy = 0;
      for (int n=1; n<=kmax; n++)
      if(By_nonzero[k][n]) // use sparsity to accelerate
      {
            const int mxy_idx = fxy_offset+fxy_map[m];
            if(half_Cy[k][n])
              half_CyFxy += half_Cy[k][n]*f.get(i  ,j  ,mxy_idx,n);
            BypmFxy      += Bypm   [k][n]*f.get(i  ,j+1,mxy_idx,n);
            BympFxy      += Bymp   [k][n]*f.get(i  ,j-1,mxy_idx,n);
      }
      Vxy.set(i,j,fxy_map[m],k, (half_CyFxy+0.5*(BypmFxy-BympFxy))/dy);
    }
  }
  }

  // fluxes to evaluate solution
  {
  #pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc-1); i++)
  for (int j=(2-mbc); j<=(my+mbc-1); j++)
  for (int m=1; m<=6; m++)
  for (int k=1; k<=kmax; k++)
  {
    double AxpVxx = 0.; double BxmpVxx = 0.;
    double half_CxVxy = 0.; double BxpmVxy = 0.; double BxmpVxy = 0.;
    double half_CyVyx = 0.; double BypmVyx = 0.; double BympVyx = 0.;
    double AypVyy = 0.; double BympVyy = 0.;
    for (int n=1; n<=kmax; n++)
    {
      if(Bx_nonzero[k][n]) // use sparsity to accelerate
      {
        AxpVxx     += Axp    [k][n]*Vxx.get(i  ,j  ,m,n);
        BxmpVxx    += Bxmp   [k][n]*Vxx.get(i-1,j  ,m,n);

        if(const int mxy_idx = fxy_map[m])
        {
          if(half_Cx[k][n])
            half_CxVxy += half_Cx[k][n]*Vxy.get(i  ,j  ,mxy_idx,n);
          BxpmVxy      += Bxpm   [k][n]*Vxy.get(i+1,j  ,mxy_idx,n);
          BxmpVxy      += Bxmp   [k][n]*Vxy.get(i-1,j  ,mxy_idx,n);
        }
      }

      if(By_nonzero[k][n]) // use sparsity to accelerate
      {
        if(const int myx_idx = fyx_map[m])
        {
          if(half_Cy[k][n]) // use sparsity to accelerate
            half_CyVyx += half_Cy[k][n]*Vyx.get(i  ,j  ,myx_idx,n);
          BypmVyx      += Bypm   [k][n]*Vyx.get(i  ,j+1,myx_idx,n);
          BympVyx      += Bymp   [k][n]*Vyx.get(i  ,j-1,myx_idx,n);
        }

        AypVyy     += Ayp    [k][n]*Vyy.get(i  ,j  ,m,n);
        BympVyy    += Bymp   [k][n]*Vyy.get(i  ,j-1,m,n);
      }
    }
    double xchange = AxpVxx-BxmpVxx+half_CxVxy + 0.5*(BxpmVxy-BxmpVxy);
    double ychange = AypVyy-BympVyy+half_CyVyx + 0.5*(BypmVyx-BympVyx);
    const double change = xchange/dx + ychange/dy;
    const int m_offset = _N11_i-1;
    const int m_idx = m+m_offset;
    Lstar.set(i,j,m_idx,k, Lstar.get(i,j,m_idx,k) + change);
  }
  }
}

// Compute explicit time-stepping diffusion operator
//
// called from ConstructL.cpp
//
void LstarExtra(dTensorBC4& q, dTensorBC4& aux, dTensorBC4& Lstar)
{
  // cout << "entering LstarExtra" << endl;
  bool no_diffusion_flux = plasmaParams.get_ion_heat_conductivity() < EPSILON;
  if(no_diffusion_flux) return;

  int mx   = q.getsize(1);
  int my   = q.getsize(2);
  int meqn = q.getsize(3);
  int kmax = q.getsize(4);
  int mbc  = q.getmbc();

  // check for NaN
  if(do_checks)
  {
    dprintf3("making checks\n");
    #pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc-1); i++)
    for (int j=(1-mbc); j<=(my+mbc-1); j++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {
      const double val = q.get(i,j,m,k);
      //assert_isnum(val);
      assert_printf(val==val,
        "i=%d, j=%d, m=%d, k=%d",
        i,j,m,k);
    }
  }
  
  // To avoid inverting the temperature tensor
  // multiple times, we store fxx, fxy, fyx, and fyy
  //
  dTensorBC4 heatff(mx,my,mheatff,kmax,mbc);

  // if we had multiple species we would call
  // L2Project_extra twice, passing in a species_code parameter
  // with values of ion_code or elc_code.
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
    const int space_order = dogParams.get_space_order();
    L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
	      space_order,space_order,space_order,space_order,
	      &q,&aux,&heatff,&heatffFunc);

  double dx = dogParamsCart2.get_dx();
  double dy = dogParamsCart2.get_dy();
  LstarExtra_helper(dx, dy, heatff, Lstar);
}

