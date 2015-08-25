#include <cmath>
#include "dogdefs.h"
#include "gas10.h"

void TenMomentFluxFunc(
    int n_offset,
    const dTensor2& Q,
    dTensor3& flux)
{
    using namespace TenMomentComponentID;
    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N11 = n_offset + _N11;
    const int n_N12 = n_offset + _N12;
    const int n_N13 = n_offset + _N13;
    const int n_N22 = n_offset + _N22;
    const int n_N23 = n_offset + _N23;
    const int n_N33 = n_offset + _N33;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        // Variables
        const double& rho = Q.get(j,n_rho);
        if(rho <= 0.)
        {
          eprintf("n_offset=%d, j=%d, rho=%24.16e", n_offset, j, rho);
        }
        const double  rho_inv = 1./rho;
        const double& M1  = Q.get(j,n_M1 );
        const double& M2  = Q.get(j,n_M2 );
        const double& M3  = Q.get(j,n_M3 );
        const double  u1  = M1*rho_inv;
        const double  u2  = M2*rho_inv;
        const double  u3  = M3*rho_inv;
    
        const double& N11 = Q.get(j,n_N11);
        const double& N12 = Q.get(j,n_N12);
        const double& N13 = Q.get(j,n_N13);
        const double& N22 = Q.get(j,n_N22);
        const double& N23 = Q.get(j,n_N23);
        const double& N33 = Q.get(j,n_N33);
        const double  P11 = N11 - M1*u1;
        const double  P12 = N12 - M1*u2;
        const double  P13 = N13 - M1*u3;
        const double  P22 = N22 - M2*u2;
        const double  P23 = N23 - M2*u3;
        //const double  P33 = N33 - M3*u3;
        
        // -------------------------------------------
        // 1-component of flux function
        // -------------------------------------------
    
        flux.set(j,n_rho,1, M1 );
        flux.set(j,n_M1 ,1, N11 );
        flux.set(j,n_M2 ,1, N12 );
        flux.set(j,n_M3 ,1, N13 );
        flux.set(j,n_N11,1, u1*N11 + 2*u1*P11 ); 
        flux.set(j,n_N12,1, u2*N11 + 2*u1*P12 );
        flux.set(j,n_N13,1, u3*N11 + 2*u1*P13 );
        flux.set(j,n_N22,1, u1*N22 + 2*u2*P12 );
        flux.set(j,n_N23,1, u1*N23 + u3*P12 + u2*P13 );
        flux.set(j,n_N33,1, u1*N33 + 2*u3*P13 );
        
        // -------------------------------------------
        // 2-component of flux function
        // -------------------------------------------
    
        flux.set(j,n_rho,2,  M2 );
        flux.set(j,n_M1 ,2,  N12 );
        flux.set(j,n_M2 ,2,  N22 );
        flux.set(j,n_M3 ,2,  N23 );
        flux.set(j,n_N11,2,  u2*N11 + 2*u1*P12 ); 
        flux.set(j,n_N12,2,  u1*N22 + 2*u2*P12 );
        flux.set(j,n_N13,2,  u2*N13 + u1*P23 + u3*P12 );
        flux.set(j,n_N22,2,  u2*N22 + 2*u2*P22 );
        flux.set(j,n_N23,2,  u3*N22 + 2*u2*P23 );
        flux.set(j,n_N33,2,  u2*N33 + 2*u3*P23 );
    }
}

void TenMomentFluxFunc1(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    using namespace TenMomentComponentID;
    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N11 = n_offset + _N11;
    const int n_N12 = n_offset + _N12;
    const int n_N13 = n_offset + _N13;
    const int n_N22 = n_offset + _N22;
    const int n_N23 = n_offset + _N23;
    const int n_N33 = n_offset + _N33;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        const double& rho = Q.get(j,n_rho);
        assert_gt(rho,0.);
        const double  rho_inv = 1./rho;
        const double& M1  = Q.get(j,n_M1 );
        const double& M2  = Q.get(j,n_M2 );
        const double& M3  = Q.get(j,n_M3 );
        const double  u1  = M1*rho_inv;
        const double  u2  = M2*rho_inv;
        const double  u3  = M3*rho_inv;
    
        const double& N11 = Q.get(j,n_N11);
        const double& N12 = Q.get(j,n_N12);
        const double& N13 = Q.get(j,n_N13);
        const double& N22 = Q.get(j,n_N22);
        const double& N23 = Q.get(j,n_N23);
        const double& N33 = Q.get(j,n_N33);
        const double  P11 = N11 - M1*u1;
        const double  P12 = N12 - M1*u2;
        const double  P13 = N13 - M1*u3;
        
        // 1-component of flux function
        //
        flux.set(j,n_rho, M1 );
        flux.set(j,n_M1 , N11 );
        flux.set(j,n_M2 , N12 );
        flux.set(j,n_M3 , N13 );
        flux.set(j,n_N11, u1*N11 + 2*u1*P11 ); 
        flux.set(j,n_N12, u2*N11 + 2*u1*P12 );
        flux.set(j,n_N13, u3*N11 + 2*u1*P13 );
        flux.set(j,n_N22, u1*N22 + 2*u2*P12 );
        flux.set(j,n_N23, u1*N23 + u3*P12 + u2*P13 );
        flux.set(j,n_N33, u1*N33 + 2*u3*P13 );
    }
}

void TenMomentFluxFunc2(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    using namespace TenMomentComponentID;
    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N11 = n_offset + _N11;
    const int n_N12 = n_offset + _N12;
    const int n_N13 = n_offset + _N13;
    const int n_N22 = n_offset + _N22;
    const int n_N23 = n_offset + _N23;
    const int n_N33 = n_offset + _N33;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        const double& rho = Q.get(j,n_rho);
        assert_gt(rho,0.);
        const double  rho_inv = 1./rho;
        const double& M1  = Q.get(j,n_M1 );
        const double& M2  = Q.get(j,n_M2 );
        const double& M3  = Q.get(j,n_M3 );
        const double  u1  = M1*rho_inv;
        const double  u2  = M2*rho_inv;
        const double  u3  = M3*rho_inv;
    
        const double& N11 = Q.get(j,n_N11);
        const double& N12 = Q.get(j,n_N12);
        const double& N13 = Q.get(j,n_N13);
        const double& N22 = Q.get(j,n_N22);
        const double& N23 = Q.get(j,n_N23);
        const double& N33 = Q.get(j,n_N33);
        const double  P12 = N12 - M1*u2;
        const double  P22 = N22 - M2*u2;
        const double  P23 = N23 - M2*u3;
        
        // 2-component of flux function
        //
        flux.set(j,n_rho,  M2 );
        flux.set(j,n_M1 ,  N12 );
        flux.set(j,n_M2 ,  N22 );
        flux.set(j,n_M3 ,  N23 );
        flux.set(j,n_N11,  u2*N11 + 2*u1*P12 ); 
        flux.set(j,n_N12,  u1*N22 + 2*u2*P12 );
        flux.set(j,n_N13,  u2*N13 + u1*P23 + u3*P12 );
        flux.set(j,n_N22,  u2*N22 + 2*u2*P22 );
        flux.set(j,n_N23,  u3*N22 + 2*u2*P23 );
        flux.set(j,n_N33,  u2*N33 + 2*u3*P23 );
    }
}

void ProjectLeftEig_10moment(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
    assert(Q.getsize(2)==W.getsize(2));
    dTensor2 Lm(10,10);

    using namespace TenMomentComponentID;
    int m_rho = _rho;
    int m_M1  = _M1 ;
    int m_M2  = _M2 ;
    int m_M3  = _M3 ;
    int m_N11 = _N11;
    int m_N12 = _N12;
    int m_N13 = _N13;
    int m_N22 = _N22;
    int m_N23 = _N23;
    int m_N33 = _N33;

    if(ixy!=1)
    {
        assert(ixy==2);

        m_M1  = _M2;
        m_M2  = _M1;

        m_N11  = _N22;
        m_N13  = _N23;
        m_N22  = _N11;
        m_N23  = _N13;
    }

    int n_rho = m_rho + n_offset;
    int n_M1  = m_M1  + n_offset;
    int n_M2  = m_M2  + n_offset;
    int n_M3  = m_M3  + n_offset;
    int n_N11 = m_N11 + n_offset;
    int n_N12 = m_N12 + n_offset;
    int n_N13 = m_N13 + n_offset;
    int n_N22 = m_N22 + n_offset;
    int n_N23 = m_N23 + n_offset;
    int n_N33 = m_N33 + n_offset;

    // Average states
    const double rho = Q_ave.get(n_rho);
    const double u1  = Q_ave.get(n_M1)/rho;
    const double u2  = Q_ave.get(n_M2)/rho;
    const double u3  = Q_ave.get(n_M3)/rho;
    const double P11 = Q_ave.get(n_N11) - rho*u1*u1;
    const double P12 = Q_ave.get(n_N12) - rho*u1*u2;
    const double P13 = Q_ave.get(n_N13) - rho*u1*u3;
    const double P22 = Q_ave.get(n_N22) - rho*u2*u2;
    const double P23 = Q_ave.get(n_N23) - rho*u2*u3;
    const double P33 = Q_ave.get(n_N33) - rho*u3*u3;
    const double c   = sqrt(P11/rho);
  
    // Set non-zero elements of matrix of left eigenvectors
    Lm.set(1,m_rho,   u1*(sq3*P11 + rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,m_M1 ,  -(sq3*P11 + 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,m_N11,   1.0/(6.0*P11*P11) );
    
    Lm.set(2,m_rho,  -(P11*u1*P12 - u2*P11*P11 + P12*u1*u1*c*rho - u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,m_M1 ,   (P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,m_M2 ,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(2,m_N11,  -P12/(2.0*P11*P11) );
    Lm.set(2,m_N12,   1.0/(2.0*P11) );

    Lm.set(3,m_rho,  -(P11*P13*u1 - u3*P11*P11 + P13*u1*u1*c*rho - u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,m_M1 ,   (P11*P13 + 2.0*u1*P13*c*rho - u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,m_M3 ,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(3,m_N11,  -P13/(2.0*P11*P11) );
    Lm.set(3,m_N13,   1.0/(2.0*P11) );
    
    Lm.set(4,m_rho,  -(rho*u1*u1 - 3.0*P11)/(3.0*P11) );
    Lm.set(4,m_M1 ,   2.0*rho*u1/(3.0*P11) );
    Lm.set(4,m_N11,  -rho/(3.0*P11) );

    Lm.set(5,m_rho,  -(u1*u1*P11*P22 - 4.0*u1*u1*P12*P12 + 6.0*u1*P12*u2*P11 - 3.0*u2*u2*P11*P11)/(3.0*P11*P11) );
    Lm.set(5,m_M1 ,   2.0*(P11*P22*u1 - 4.0*u1*P12*P12 + 3.0*u2*P11*P12)/(3.0*P11*P11) );
    Lm.set(5,m_M2 ,   (2.0*u1*P12 - 2.0*u2*P11)/P11 );
    Lm.set(5,m_N11,  -(P11*P22 - 4.0*P12*P12)/(3.0*P11*P11) );
    Lm.set(5,m_N12,  -2.0*P12/P11 );
    Lm.set(5,m_N22,   1.0 );
    
    Lm.set(6,m_rho,  -(u1*u1*P11*P23 - 4.0*P13*u1*u1*P12 + 3.0*u1*P13*u2*P11 + 3.0*u1*P12*u3*P11 
		   - 3.0*u2*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(6,m_M1 ,   (2.0*P11*P23*u1 - 8*u1*P13*P12 + 3.0*u2*P11*P13 + 3.0*u3*P11*P12)/(3.0*P11*P11) );
    Lm.set(6,m_M2 ,   (u1*P13 - u3*P11)/P11 );
    Lm.set(6,m_M3 ,   (u1*P12 - u2*P11)/P11 );
    Lm.set(6,m_N11,  -(P11*P23 - 4.0*P13*P12)/(3.0*P11*P11) );
    Lm.set(6,m_N12,  -P13/P11 );
    Lm.set(6,m_N13,  -P12/P11 );
    Lm.set(6,m_N23,   1.0 );
    
    Lm.set(7,m_rho,  -(u1*u1*P11*P33 - 4.0*u1*u1*P13*P13 + 6.0*u1*P13*u3*P11 - 3.0*u3*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(7,m_M1 ,   2.0*(P11*P33*u1 - 4.0*u1*P13*P13 + 3.0*u3*P11*P13)/(3.0*P11*P11) );
    Lm.set(7,m_M3 ,   (2.0*u1*P13 - 2.0*u3*P11)/P11 );
    Lm.set(7,m_N11,  -(P11*P33 - 4.0*P13*P13)/(3.0*P11*P11) );
    Lm.set(7,m_N13,  -2.0*P13/P11 );
    Lm.set(7,m_N33,  1.0 );

    Lm.set(8,m_rho,   (P11*u1*P12 - u2*P11*P11 - P12*u1*u1*c*rho + u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,m_M1 ,   (-P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,m_M2 ,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(8,m_N11,  -P12/(2.0*P11*P11) );
    Lm.set(8,m_N12,   1.0/(2.0*P11) );

    Lm.set(9,m_rho,   (P11*P13*u1 - u3*P11*P11 - P13*u1*u1*c*rho + u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,m_M1 ,  -(P11*P13 - 2.0*u1*P13*c*rho + u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,m_M3 ,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(9,m_N11,  -P13/(2.0*P11*P11) );
    Lm.set(9,m_N13,   1.0/(2.0*P11) );
    
    Lm.set(10,m_rho,   u1*(rho*u1*c - sq3*P11)/(6.0*c*P11*P11*rho) );
    Lm.set(10,m_M1 ,   (sq3*P11 - 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(10,m_N11,   1.0/(6.0*P11*P11) );    
  
    // eigen indices
    //
    int n1  = 1  + n_offset;
    int n2  = 2  + n_offset;
    int n3  = 3  + n_offset;
    int n4  = 4  + n_offset;
    int n5  = 5  + n_offset;
    int n6  = 6  + n_offset;
    int n7  = 7  + n_offset;
    int n8  = 8  + n_offset;
    int n9  = 9  + n_offset;
    int n10 = 10 + n_offset;

    // Project onto left eigenvectors
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(n1,k,
            Lm.get(1,m_rho)*Q.get(n_rho,k)
          + Lm.get(1,m_M1 )*Q.get(n_M1,k)
          + Lm.get(1,m_N11)*Q.get(n_N11,k) );

        W.set(n2,k,
            Lm.get(2,m_rho)*Q.get(n_rho,k)
          + Lm.get(2,m_M1 )*Q.get(n_M1,k)
          + Lm.get(2,m_M2 )*Q.get(n_M2,k) 
          + Lm.get(2,m_N11)*Q.get(n_N11,k)
          + Lm.get(2,m_N12)*Q.get(n_N12,k) );
        
        W.set(n3,k,
            Lm.get(3,m_rho)*Q.get(n_rho,k)
          + Lm.get(3,m_M1 )*Q.get(n_M1,k)
          + Lm.get(3,m_M3 )*Q.get(n_M3,k) 
          + Lm.get(3,m_N11)*Q.get(n_N11,k)
          + Lm.get(3,m_N13)*Q.get(n_N13,k) );

        W.set(n4,k,
            Lm.get(4,m_rho)*Q.get(n_rho,k)
          + Lm.get(4,m_M1 )*Q.get(n_M1,k)
          + Lm.get(4,m_N11)*Q.get(n_N11,k) );

        W.set(n5,k,
            Lm.get(5,m_rho)*Q.get(n_rho,k)
          + Lm.get(5,m_M1 )*Q.get(n_M1,k)
          + Lm.get(5,m_M2 )*Q.get(n_M2,k) 
          + Lm.get(5,m_N11)*Q.get(n_N11,k)
          + Lm.get(5,m_N12)*Q.get(n_N12,k)
          +              Q.get(n_N22,k) );

        W.set(n6,k,
            Lm.get(6,m_rho)*Q.get(n_rho,k)
          + Lm.get(6,m_M1 )*Q.get(n_M1,k)
          + Lm.get(6,m_M2 )*Q.get(n_M2,k) 
          + Lm.get(6,m_M3 )*Q.get(n_M3,k)
          + Lm.get(6,m_N11)*Q.get(n_N11,k)
          + Lm.get(6,m_N12)*Q.get(n_N12,k)
          + Lm.get(6,m_N13)*Q.get(n_N13,k)
          +              Q.get(n_N23,k) );

        W.set(n7,k,
            Lm.get(7,m_rho)*Q.get(n_rho,k)
          + Lm.get(7,m_M1 )*Q.get(n_M1,k)
          + Lm.get(7,m_M3 )*Q.get(n_M3,k) 
          + Lm.get(7,m_N11)*Q.get(n_N11,k)
          + Lm.get(7,m_N13)*Q.get(n_N13,k)
          +              Q.get(n_N33,k) );

        W.set(n8,k,
            Lm.get(8,m_rho)*Q.get(n_rho,k)
          + Lm.get(8,m_M1 )*Q.get(n_M1,k)
          + Lm.get(8,m_M2 )*Q.get(n_M2,k) 
          + Lm.get(8,m_N11)*Q.get(n_N11,k)
          + Lm.get(8,m_N12)*Q.get(n_N12,k) );

        W.set(n9,k,
            Lm.get(9,m_rho)*Q.get(n_rho,k)
          + Lm.get(9,m_M1 )*Q.get(n_M1,k)
          + Lm.get(9,m_M3 )*Q.get(n_M3,k) 
          + Lm.get(9,m_N11)*Q.get(n_N11,k)
          + Lm.get(9,m_N13)*Q.get(n_N13,k) );

        W.set(n10,k,
            Lm.get(10,m_rho)*Q.get(n_rho,k)
          + Lm.get(10,m_M1 )*Q.get(n_M1,k)
          + Lm.get(10,m_N11)*Q.get(n_N11,k) );
    }
}

void ProjectRightEig_10moment(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
  dTensor2 Rm(10,10);

  using namespace TenMomentComponentID;
  int m_rho = _rho; // m1
  int m_M1  = _M1 ; // m2
  int m_M2  = _M2 ; // m3
  int m_M3  = _M3 ; // m4
  int m_N11 = _N11; // m5
  int m_N12 = _N12; // m6
  int m_N13 = _N13; // m7
  int m_N22 = _N22; // m8
  int m_N23 = _N23; // m9
  int m_N33 = _N33; // m10

  if(ixy!=1)
  {
      assert(ixy==2);

      m_M1  = _M2;
      m_M2  = _M1;

      m_N11  = _N22;
      m_N13  = _N23;
      m_N22  = _N11;
      m_N23  = _N13;
  }

  int n_rho = m_rho + n_offset;
  int n_M1  = m_M1  + n_offset;
  int n_M2  = m_M2  + n_offset;
  int n_M3  = m_M3  + n_offset;
  int n_N11 = m_N11 + n_offset;
  int n_N12 = m_N12 + n_offset;
  int n_N13 = m_N13 + n_offset;
  int n_N22 = m_N22 + n_offset;
  int n_N23 = m_N23 + n_offset;
  int n_N33 = m_N33 + n_offset;

  // Average states
  const double rho = Q_ave.get(n_rho);
  const double u1  = Q_ave.get(n_M1)/rho;
  const double u2  = Q_ave.get(n_M2)/rho;
  const double u3  = Q_ave.get(n_M3)/rho;
  const double P11 = Q_ave.get(n_N11) - rho*u1*u1;
  const double P12 = Q_ave.get(n_N12) - rho*u1*u2;
  const double P13 = Q_ave.get(n_N13) - rho*u1*u3;
  const double P22 = Q_ave.get(n_N22) - rho*u2*u2;
  const double P23 = Q_ave.get(n_N23) - rho*u2*u3;
  const double P33 = Q_ave.get(n_N33) - rho*u3*u3;
  const double c   = sqrt(P11/rho);

  // Set non-zero elements of matrix of right eigenvectors
  Rm.set(m_rho,1,        P11*rho );
  Rm.set(m_M1 ,1,     u1*P11*rho - rho*sq3*c*P11 );
  Rm.set(m_M2 ,1,     u2*P11*rho - rho*sq3*c*P12 );
  Rm.set(m_M3 ,1,     u3*P11*rho - rho*sq3*c*P13 );
  Rm.set(m_N11,1,  u1*u1*P11*rho - 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
  Rm.set(m_N12,1,  u1*u2*P11*rho - rho*u2*sq3*c*P11 - rho*u1*sq3*c*P12 + 3.0*P11*P12 );
  Rm.set(m_N13,1,  u1*u3*P11*rho - rho*u3*sq3*c*P11 - rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
  Rm.set(m_N22,1,  u2*u2*P11*rho - 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
  Rm.set(m_N23,1,  u2*u3*P11*rho - rho*u3*sq3*c*P12 - rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
  Rm.set(m_N33,1,  u3*u3*P11*rho - 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
  
  Rm.set(m_M2 ,2,  -rho*c );
  Rm.set(m_N12,2,  -rho*u1*c + P11 );
  Rm.set(m_N22,2,  -2.0*rho*u2*c + 2.0*P12 );
  Rm.set(m_N23,2,  -rho*u3*c + P13 );

  Rm.set(m_M3 ,3,  -rho*c ); 
  Rm.set(m_N13,3,  -rho*u1*c + P11 );
  Rm.set(m_N23,3,  -rho*u2*c + P12 );
  Rm.set(m_N33,3,  -2.0*rho*u3*c + 2.0*P13 );
  
  Rm.set(m_rho,4,  1.0 );
  Rm.set(m_M1 ,4,  u1 );
  Rm.set(m_M2 ,4,  u2 );
  Rm.set(m_M3 ,4,  u3 );
  Rm.set(m_N11,4,  u1*u1 );
  Rm.set(m_N12,4,  u1*u2 ); 
  Rm.set(m_N13,4,  u1*u3 ); 
  Rm.set(m_N22,4,  u2*u2 );
  Rm.set(m_N23,4,  u2*u3 ); 
  Rm.set(m_N33,4,  u3*u3 );
  
  Rm.set(m_N22,5,  1.0 );
  
  Rm.set(m_N23,6,  1.0 );     

  Rm.set(m_N33,7,  1.0 );
  
  Rm.set(m_M2 ,8,  rho*c );
  Rm.set(m_N12,8,  rho*u1*c + P11 );
  Rm.set(m_N22,8,  2.0*rho*u2*c + 2.0*P12 );
  Rm.set(m_N23,8,  rho*u3*c + P13 );
  
  Rm.set(m_M3 ,9,  rho*c ); 
  Rm.set(m_N13,9,  rho*u1*c + P11 );
  Rm.set(m_N23,9,  rho*u2*c + P12 );
  Rm.set(m_N33,9,  2.0*rho*u3*c + 2.0*P13 );

  Rm.set(m_rho,10,        P11*rho );
  Rm.set(m_M1 ,10,     u1*P11*rho + rho*sq3*c*P11 );
  Rm.set(m_M2 ,10,     u2*P11*rho + rho*sq3*c*P12 );
  Rm.set(m_M3 ,10,     u3*P11*rho + rho*sq3*c*P13 );
  Rm.set(m_N11,10,  u1*u1*P11*rho + 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
  Rm.set(m_N12,10,  u1*u2*P11*rho + rho*u2*sq3*c*P11 + rho*u1*sq3*c*P12 + 3.0*P11*P12 );
  Rm.set(m_N13,10,  u1*u3*P11*rho + rho*u3*sq3*c*P11 + rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
  Rm.set(m_N22,10,  u2*u2*P11*rho + 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
  Rm.set(m_N23,10,  u2*u3*P11*rho + rho*u3*sq3*c*P12 + rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
  Rm.set(m_N33,10,  u3*u3*P11*rho + 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
  
  // eigen indices
  //
  int n1  = 1  + n_offset;
  int n2  = 2  + n_offset;
  int n3  = 3  + n_offset;
  int n4  = 4  + n_offset;
  int n5  = 5  + n_offset;
  int n6  = 6  + n_offset;
  int n7  = 7  + n_offset;
  int n8  = 8  + n_offset;
  int n9  = 9  + n_offset;
  int n10 = 10 + n_offset;

  // Project onto right eigenvectors
  for (int k=1; k<=Q.getsize(2); k++)
  {
      Q.set(n_rho,k,
          Rm.get(m_rho,1) *W.get(n1,k)
        +                  W.get(n4,k)
        + Rm.get(m_rho,10)*W.get(n10,k) );
      Q.set(n_M1,k,
          Rm.get(m_M1,1) *W.get(n1,k)
        + Rm.get(m_M1,4) *W.get(n4,k)
        + Rm.get(m_M1,10)*W.get(n10,k) );
      Q.set(n_M2,k,
          Rm.get(m_M2,1) *W.get(n1,k)
        + Rm.get(m_M2,2) *W.get(n2,k)
        + Rm.get(m_M2,4) *W.get(n4,k) 
        + Rm.get(m_M2,8) *W.get(n8,k)
        + Rm.get(m_M2,10)*W.get(n10,k) );
      Q.set(n_M3,k,
          Rm.get(m_M3,1) *W.get(n1,k)
        + Rm.get(m_M3,3) *W.get(n3,k)
        + Rm.get(m_M3,4) *W.get(n4,k) 
        + Rm.get(m_M3,9) *W.get(n9,k)
        + Rm.get(m_M3,10)*W.get(n10,k) );
      Q.set(n_N11,k,
          Rm.get(m_N11,1) *W.get(n1,k)
        + Rm.get(m_N11,4) *W.get(n4,k)
        + Rm.get(m_N11,10)*W.get(n10,k) );
      Q.set(n_N12,k,
          Rm.get(m_N12,1) *W.get(n1,k)
        + Rm.get(m_N12,2) *W.get(n2,k)
        + Rm.get(m_N12,4) *W.get(n4,k) 
        + Rm.get(m_N12,8) *W.get(n8,k)
        + Rm.get(m_N12,10)*W.get(n10,k) );
      Q.set(n_N13,k,
          Rm.get(m_N13,1) *W.get(n1,k)
        + Rm.get(m_N13,3) *W.get(n3,k)
        + Rm.get(m_N13,4) *W.get(n4,k) 
        + Rm.get(m_N13,9) *W.get(n9,k)
        + Rm.get(m_N13,10)*W.get(n10,k) );
      Q.set(n_N22,k,
          Rm.get(m_N22,1) *W.get(n1,k)
        + Rm.get(m_N22,2) *W.get(n2,k)
        + Rm.get(m_N22,4) *W.get(n4,k) 
        +                  W.get(n5,k)
        + Rm.get(m_N22,8) *W.get(n8,k)
        + Rm.get(m_N22,10)*W.get(n10,k) );
      Q.set(n_N23,k,
          Rm.get(m_N23,1) *W.get(n1,k)
        + Rm.get(m_N23,2) *W.get(n2,k)
        + Rm.get(m_N23,3) *W.get(n3,k) 
        + Rm.get(m_N23,4) *W.get(n4,k)
        +                  W.get(n6,k)
        + Rm.get(m_N23,8) *W.get(n8,k) 
        + Rm.get(m_N23,9) *W.get(n9,k)
        + Rm.get(m_N23,10)*W.get(n10,k) );
      Q.set(n_N33,k,
          Rm.get(m_N33,1) *W.get(n1,k)
        + Rm.get(m_N33,3) *W.get(n3,k)
        + Rm.get(m_N33,4) *W.get(n4,k) 
        +                  W.get(n7,k)
        + Rm.get(m_N33,9) *W.get(n9,k)
        + Rm.get(m_N33,10)*W.get(n10,k) );
  }
}

