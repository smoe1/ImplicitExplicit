#include <cmath>
#include "dogdefs.h"
#include "gas20.h"

void TwentyMomentFluxFunc(
    int n_offset,
    const dTensor2& Q,
    dTensor3& flux)
{
    using namespace TwentyMomentComponentID;
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
    const int n_Q111 = n_offset + _Q111;
    const int n_Q112 = n_offset + _Q112;
    const int n_Q113 = n_offset + _Q113;
    const int n_Q122 = n_offset + _Q122;
    const int n_Q123 = n_offset + _Q123;
    const int n_Q133 = n_offset + _Q133;
    const int n_Q222 = n_offset + _Q222;
    const int n_Q223 = n_offset + _Q223;
    const int n_Q233 = n_offset + _Q233;
    const int n_Q333 = n_offset + _Q333;

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
        const double uu11 = u1*u1;
        const double uu12 = u1*u2;
        const double uu13 = u1*u3;
        const double uu22 = u2*u2;
        const double uu23 = u2*u3;
        const double uu33 = u3*u3;
    
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
        const double  P33 = N33 - M3*u3;
        const double  T11 = P11*rho_inv;
        const double  T12 = P12*rho_inv;
        const double  T13 = P13*rho_inv;
        const double  T22 = P22*rho_inv;
        const double  T23 = P23*rho_inv;
        //const double  T33 = P33*rho_inv;
        const double& QN111 = Q.get(j,n_Q111);
        const double& QN112 = Q.get(j,n_Q112);
        const double& QN113 = Q.get(j,n_Q113);
        const double& QN122 = Q.get(j,n_Q122);
        const double& QN123 = Q.get(j,n_Q123);
        const double& QN133 = Q.get(j,n_Q133);
        const double& QN222 = Q.get(j,n_Q222);
        const double& QN223 = Q.get(j,n_Q223);
        const double& QN233 = Q.get(j,n_Q233);
        const double& QN333 = Q.get(j,n_Q333);

        // *************************************************************
        // Use Gaussian-based closure to compute primitive fourth moment
        //
        //   R = Sym3(rho*Theta*Theta)
        //   Rijkl = Tij*Pkl+Tik*Pjl+Til*Pjk;
        const double R1111 = 3*T11*P11;
        const double R1112 = 3*T11*P12;
        const double R1113 = 3*T11*P13;
        const double R1122 = T11*P22+2*T12*P12;
        const double R1123 = T11*P23+2*T12*P13;
        const double R1133 = T11*P33+2*T13*P13;
        const double R1222 = 3*T12*P22;
        const double R1223 = 2*T12*P23+T13*P22;
        const double R1233 = T12*P33+2*T13*P23;
        const double R1333 = 3*T13*P33;
        // *************************************************************

        // Use primitive moments and lower conserved moments
        // to compute conserved fourth moment:
        //   RN = R + Sym4(u*QN)-Sym6(uu*N)+3*rho*uuuu
        //
        // Sym4uQNijkl = ui*QNjkl+uj*QNikl+uk*QNijl+ul*QNijk;
        //
        const double Sym4uQN1111 = 4*u1*QN111;
        const double Sym4uQN1112 = 3*u1*QN112+u2*QN111;
        const double Sym4uQN1113 = 3*u1*QN113+u3*QN111;
        const double Sym4uQN1122 = 2*u1*QN122+2*u2*QN112;
        const double Sym4uQN1123 = 2*u1*QN123 + u2*QN113 + u3*QN112;
        const double Sym4uQN1133 = 2*u1*QN133+2*u3*QN113;
        const double Sym4uQN1222 = u1*QN222+3*u2*QN122;
        const double Sym4uQN1223 = u1*QN223+2*u2*QN123+u3*QN122;
        const double Sym4uQN1233 = u1*QN233+u2*QN133+2*u3*QN123;
        const double Sym4uQN1333 = u1*QN333+3*u3*QN133;
        //
        // Sym6uuNijkl = uuij*Nkl + uuik*Njl + uuil*Njk
        //             + uujk*Nil + uujl*Nik + uukl*Nij;
        const double Sym6uuN1111 = 6*uu11*N11;
        const double Sym6uuN1112 = 3*uu11*N12+3*uu12*N11;
        const double Sym6uuN1113 = 3*uu11*N13+3*uu13*N11;
        const double Sym6uuN1122 = uu11*N22 + uu22*N11 + 4*uu12*N12;
        const double Sym6uuN1123 = uu11*N23 + 2*uu12*N13 + 2*uu13*N12 + uu23*N11;
        const double Sym6uuN1133 = uu11*N33 + uu33*N11 + 4*uu13*N13;
        const double Sym6uuN1222 = 3*uu12*N22 + 3*uu22*N12;
        const double Sym6uuN1223 = 2*uu12*N23 + uu13*N22 + uu22*N13 + 2*uu23*N12;
        const double Sym6uuN1233 = uu12*N33 + 2*uu13*N23 + 2*uu23*N13 + uu33*N12;
        const double Sym6uuN1333 = 3*uu13*N33 + 3*uu33*N13;
        //
        const double RN1111 = R1111 + Sym4uQN1111 - Sym6uuN1111 + 3*rho*uu11*uu11;
        const double RN1112 = R1112 + Sym4uQN1112 - Sym6uuN1112 + 3*rho*uu11*uu12;
        const double RN1113 = R1113 + Sym4uQN1113 - Sym6uuN1113 + 3*rho*uu11*uu13;
        const double RN1122 = R1122 + Sym4uQN1122 - Sym6uuN1122 + 3*rho*uu11*uu22;
        const double RN1123 = R1123 + Sym4uQN1123 - Sym6uuN1123 + 3*rho*uu11*uu23;
        const double RN1133 = R1133 + Sym4uQN1133 - Sym6uuN1133 + 3*rho*uu11*uu33;
        const double RN1222 = R1222 + Sym4uQN1222 - Sym6uuN1222 + 3*rho*uu12*uu22;
        const double RN1223 = R1223 + Sym4uQN1223 - Sym6uuN1223 + 3*rho*uu12*uu23;
        const double RN1233 = R1233 + Sym4uQN1233 - Sym6uuN1233 + 3*rho*uu12*uu33;
        const double RN1333 = R1333 + Sym4uQN1333 - Sym6uuN1333 + 3*rho*uu13*uu33;
        
        // -------------------------------------------
        // 1-component of flux function
        // -------------------------------------------
    
        flux.set(j,n_rho,1, M1 );
        flux.set(j,n_M1 ,1, N11 );
        flux.set(j,n_M2 ,1, N12 );
        flux.set(j,n_M3 ,1, N13 );
        flux.set(j,n_N11,1, QN111);
        flux.set(j,n_N12,1, QN112);
        flux.set(j,n_N13,1, QN113);
        flux.set(j,n_N22,1, QN122);
        flux.set(j,n_N23,1, QN123);
        flux.set(j,n_N33,1, QN133);
        flux.set(j,n_Q111,1, RN1111);
        flux.set(j,n_Q112,1, RN1112);
        flux.set(j,n_Q113,1, RN1113);
        flux.set(j,n_Q122,1, RN1122);
        flux.set(j,n_Q123,1, RN1123);
        flux.set(j,n_Q133,1, RN1133);
        flux.set(j,n_Q222,1, RN1222);
        flux.set(j,n_Q223,1, RN1223);
        flux.set(j,n_Q233,1, RN1233);
        flux.set(j,n_Q333,1, RN1333);
        
        // -------------------------------------------
        // 2-component of flux function
        // -------------------------------------------

        // Use Gaussian-based closure to compute primitive fourth moment
        //
        //   R = Sym3(rho*Theta*Theta)
        //   Rijkl = Tij*Pkl+Tik*Pjl+Til*Pjk;
        const double R2222 = 3*T22*P22;
        const double R2223 = 3*T22*P23;
        const double R2233 = T22*P33+2*T23*P23;
        const double R2333 = 3*T23*P33;

        // Use primitive moments and lower conserved moments
        // to compute conserved fourth moment:
        //   RN = R + Sym4(u*QN)-Sym6(uu*N)+3*rho*uuuu
        //
        // Sym4uQNijkl = ui*QNjkl+uj*QNikl+uk*QNijl+ul*QNijk;
        //
        const double Sym4uQN2222 = 4*u2*QN222;
        const double Sym4uQN2223 = 3*u2*QN223+u3*QN222;
        const double Sym4uQN2233 = 2*u2*QN233+2*u3*QN223;
        const double Sym4uQN2333 = u2*QN333+3*u3*QN233;
        //
        // Sym6uuNijkl = uuij*Nkl + uuik*Njl + uuil*Njk
        //             + uujk*Nil + uujl*Nik + uukl*Nij;
        const double Sym6uuN2222 = 6*uu22*N22;
        const double Sym6uuN2223 = 3*uu22*N23 + 3*uu23*N22;
        const double Sym6uuN2233 = uu22*N33 + 4*uu23*N23 + uu33*N22;
        const double Sym6uuN2333 = 3*uu23*N33 + 3*uu33*N23;
        //
        const double RN2222 = R2222 + Sym4uQN2222 - Sym6uuN2222 + 3*rho*uu22*uu22;
        const double RN2223 = R2223 + Sym4uQN2223 - Sym6uuN2223 + 3*rho*uu22*uu23;
        const double RN2233 = R2233 + Sym4uQN2233 - Sym6uuN2233 + 3*rho*uu22*uu33;
        const double RN2333 = R2333 + Sym4uQN2333 - Sym6uuN2333 + 3*rho*uu23*uu33;
    
        flux.set(j,n_rho,2,  M2 );
        flux.set(j,n_M1 ,2,  N12 );
        flux.set(j,n_M2 ,2,  N22 );
        flux.set(j,n_M3 ,2,  N23 );
        flux.set(j,n_N11,2,  QN112);
        flux.set(j,n_N12,2,  QN122);
        flux.set(j,n_N13,2,  QN123);
        flux.set(j,n_N22,2,  QN222);
        flux.set(j,n_N23,2,  QN223);
        flux.set(j,n_N33,2,  QN233);
        flux.set(j,n_Q111,2, RN1112);
        flux.set(j,n_Q112,2, RN1122);
        flux.set(j,n_Q113,2, RN1123);
        flux.set(j,n_Q122,2, RN1222);
        flux.set(j,n_Q123,2, RN1223);
        flux.set(j,n_Q133,2, RN1233);
        flux.set(j,n_Q222,2, RN2222);
        flux.set(j,n_Q223,2, RN2223);
        flux.set(j,n_Q233,2, RN2233);
        flux.set(j,n_Q333,2, RN2333);
    }
}

void TwentyMomentFluxFunc1(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    using namespace TwentyMomentComponentID;
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
    const int n_Q111 = n_offset + _Q111;
    const int n_Q112 = n_offset + _Q112;
    const int n_Q113 = n_offset + _Q113;
    const int n_Q122 = n_offset + _Q122;
    const int n_Q123 = n_offset + _Q123;
    const int n_Q133 = n_offset + _Q133;
    const int n_Q222 = n_offset + _Q222;
    const int n_Q223 = n_offset + _Q223;
    const int n_Q233 = n_offset + _Q233;
    const int n_Q333 = n_offset + _Q333;

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
        const double uu11 = u1*u1;
        const double uu12 = u1*u2;
        const double uu13 = u1*u3;
        const double uu22 = u2*u2;
        const double uu23 = u2*u3;
        const double uu33 = u3*u3;
    
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
        const double  P33 = N33 - M3*u3;
        const double  T11 = P11*rho_inv;
        const double  T12 = P12*rho_inv;
        const double  T13 = P13*rho_inv;
        const double  T22 = P22*rho_inv;
        const double  T23 = P23*rho_inv;
        //const double  T33 = P33*rho_inv;
        const double& QN111 = Q.get(j,n_Q111);
        const double& QN112 = Q.get(j,n_Q112);
        const double& QN113 = Q.get(j,n_Q113);
        const double& QN122 = Q.get(j,n_Q122);
        const double& QN123 = Q.get(j,n_Q123);
        const double& QN133 = Q.get(j,n_Q133);
        const double& QN222 = Q.get(j,n_Q222);
        const double& QN223 = Q.get(j,n_Q223);
        const double& QN233 = Q.get(j,n_Q233);
        const double& QN333 = Q.get(j,n_Q333);

        // Use Gaussian-based closure to compute primitive fourth moment
        //
        //   R = Sym3(rho*Theta*Theta)
        //   Rijkl = Tij*Pkl+Tik*Pjl+Til*Pjk;
        const double R1111 = 3*T11*P11;
        const double R1112 = 3*T11*P12;
        const double R1113 = 3*T11*P13;
        const double R1122 = T11*P22+2*T12*P12;
        const double R1123 = T11*P23+2*T12*P13;
        const double R1133 = T11*P33+2*T13*P13;
        const double R1222 = 3*T12*P22;
        const double R1223 = 2*T12*P23+T13*P22;
        const double R1233 = T12*P33+2*T13*P23;
        const double R1333 = 3*T13*P33;

        // Use primitive moments and lower conserved moments
        // to compute conserved fourth moment:
        //   RN = R + Sym4(u*QN)-Sym6(uu*N)+3*rho*uuuu
        //
        // Sym4uQNijkl = ui*QNjkl+uj*QNikl+uk*QNijl+ul*QNijk;
        //
        const double Sym4uQN1111 = 4*u1*QN111;
        const double Sym4uQN1112 = 3*u1*QN112+u2*QN111;
        const double Sym4uQN1113 = 3*u1*QN113+u3*QN111;
        const double Sym4uQN1122 = 2*u1*QN122+2*u2*QN112;
        const double Sym4uQN1123 = 2*u1*QN123 + u2*QN113 + u3*QN112;
        const double Sym4uQN1133 = 2*u1*QN133+2*u3*QN113;
        const double Sym4uQN1222 = u1*QN222+3*u2*QN122;
        const double Sym4uQN1223 = u1*QN223+2*u2*QN123+u3*QN122;
        const double Sym4uQN1233 = u1*QN233+u2*QN133+2*u3*QN123;
        const double Sym4uQN1333 = u1*QN333+3*u3*QN133;
        //
        // Sym6uuNijkl = uuij*Nkl + uuik*Njl + uuil*Njk
        //             + uujk*Nil + uujl*Nik + uukl*Nij;
        const double Sym6uuN1111 = 6*uu11*N11;
        const double Sym6uuN1112 = 3*uu11*N12+3*uu12*N11;
        const double Sym6uuN1113 = 3*uu11*N13+3*uu13*N11;
        const double Sym6uuN1122 = uu11*N22 + uu22*N11 + 4*uu12*N12;
        const double Sym6uuN1123 = uu11*N23 + 2*uu12*N13 + 2*uu13*N12 + uu23*N11;
        const double Sym6uuN1133 = uu11*N33 + uu33*N11 + 4*uu13*N13;
        const double Sym6uuN1222 = 3*uu12*N22 + 3*uu22*N12;
        const double Sym6uuN1223 = 2*uu12*N23 + uu13*N22 + uu22*N13 + 2*uu23*N12;
        const double Sym6uuN1233 = uu12*N33 + 2*uu13*N23 + 2*uu23*N13 + uu33*N12;
        const double Sym6uuN1333 = 3*uu13*N33 + 3*uu33*N13;
        //
        const double RN1111 = R1111 + Sym4uQN1111 - Sym6uuN1111 + 3*rho*uu11*uu11;
        const double RN1112 = R1112 + Sym4uQN1112 - Sym6uuN1112 + 3*rho*uu11*uu12;
        const double RN1113 = R1113 + Sym4uQN1113 - Sym6uuN1113 + 3*rho*uu11*uu13;
        const double RN1122 = R1122 + Sym4uQN1122 - Sym6uuN1122 + 3*rho*uu11*uu22;
        const double RN1123 = R1123 + Sym4uQN1123 - Sym6uuN1123 + 3*rho*uu11*uu23;
        const double RN1133 = R1133 + Sym4uQN1133 - Sym6uuN1133 + 3*rho*uu11*uu33;
        const double RN1222 = R1222 + Sym4uQN1222 - Sym6uuN1222 + 3*rho*uu12*uu22;
        const double RN1223 = R1223 + Sym4uQN1223 - Sym6uuN1223 + 3*rho*uu12*uu23;
        const double RN1233 = R1233 + Sym4uQN1233 - Sym6uuN1233 + 3*rho*uu12*uu33;
        const double RN1333 = R1333 + Sym4uQN1333 - Sym6uuN1333 + 3*rho*uu13*uu33;
        
        // 1-component of flux function
        //
        flux.set(j,n_rho, M1 );
        flux.set(j,n_M1 , N11 );
        flux.set(j,n_M2 , N12 );
        flux.set(j,n_M3 , N13 );
        flux.set(j,n_N11, QN111);
        flux.set(j,n_N12, QN112);
        flux.set(j,n_N13, QN113);
        flux.set(j,n_N22, QN122);
        flux.set(j,n_N23, QN123);
        flux.set(j,n_N33, QN133);
        flux.set(j,n_Q111, RN1111);
        flux.set(j,n_Q112, RN1112);
        flux.set(j,n_Q113, RN1113);
        flux.set(j,n_Q122, RN1122);
        flux.set(j,n_Q123, RN1123);
        flux.set(j,n_Q133, RN1133);
        flux.set(j,n_Q222, RN1222);
        flux.set(j,n_Q223, RN1223);
        flux.set(j,n_Q233, RN1233);
        flux.set(j,n_Q333, RN1333);
    }
}

void TwentyMomentFluxFunc2(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    using namespace TwentyMomentComponentID;
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
    const int n_Q111 = n_offset + _Q111;
    const int n_Q112 = n_offset + _Q112;
    const int n_Q113 = n_offset + _Q113;
    const int n_Q122 = n_offset + _Q122;
    const int n_Q123 = n_offset + _Q123;
    const int n_Q133 = n_offset + _Q133;
    const int n_Q222 = n_offset + _Q222;
    const int n_Q223 = n_offset + _Q223;
    const int n_Q233 = n_offset + _Q233;
    const int n_Q333 = n_offset + _Q333;

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
        const double uu11 = u1*u1;
        const double uu12 = u1*u2;
        const double uu13 = u1*u3;
        const double uu22 = u2*u2;
        const double uu23 = u2*u3;
        const double uu33 = u3*u3;
    
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
        const double  P33 = N33 - M3*u3;
        const double  T11 = P11*rho_inv;
        const double  T12 = P12*rho_inv;
        const double  T13 = P13*rho_inv;
        const double  T22 = P22*rho_inv;
        const double  T23 = P23*rho_inv;
        const double& QN111 = Q.get(j,n_Q111);
        const double& QN112 = Q.get(j,n_Q112);
        const double& QN113 = Q.get(j,n_Q113);
        const double& QN122 = Q.get(j,n_Q122);
        const double& QN123 = Q.get(j,n_Q123);
        const double& QN133 = Q.get(j,n_Q133);
        const double& QN222 = Q.get(j,n_Q222);
        const double& QN223 = Q.get(j,n_Q223);
        const double& QN233 = Q.get(j,n_Q233);
        const double& QN333 = Q.get(j,n_Q333);

        // Use Gaussian-based closure to compute primitive fourth moment
        //
        //   R = Sym3(rho*Theta*Theta)
        //   Rijkl = Tij*Pkl+Tik*Pjl+Til*Pjk;
        const double R1112 = 3*T11*P12;
        const double R1122 = T11*P22+2*T12*P12;
        const double R1123 = T11*P23+2*T12*P13;
        const double R1222 = 3*T12*P22;
        const double R1223 = 2*T12*P23+T13*P22;
        const double R1233 = T12*P33+2*T13*P23;
        const double R2222 = 3*T22*P22;
        const double R2223 = 3*T22*P23;
        const double R2233 = T22*P33+2*T23*P23;
        const double R2333 = 3*T23*P33;

        // Use primitive moments and lower conserved moments
        // to compute conserved fourth moment:
        //   RN = R + Sym4(u*QN)-Sym6(uu*N)+3*rho*uuuu
        //
        // Sym4uQNijkl = ui*QNjkl+uj*QNikl+uk*QNijl+ul*QNijk;
        //
        const double Sym4uQN1112 = 3*u1*QN112+u2*QN111;
        const double Sym4uQN1122 = 2*u1*QN122+2*u2*QN112;
        const double Sym4uQN1123 = 2*u1*QN123 + u2*QN113 + u3*QN112;
        const double Sym4uQN1222 = u1*QN222+3*u2*QN122;
        const double Sym4uQN1223 = u1*QN223+2*u2*QN123+u3*QN122;
        const double Sym4uQN1233 = u1*QN233+u2*QN133+2*u3*QN123;
        const double Sym4uQN2222 = 4*u2*QN222;
        const double Sym4uQN2223 = 3*u2*QN223+u3*QN222;
        const double Sym4uQN2233 = 2*u2*QN233+2*u3*QN223;
        const double Sym4uQN2333 = u2*QN333+3*u3*QN233;
        //
        // Sym6uuNijkl = uuij*Nkl + uuik*Njl + uuil*Njk
        //             + uujk*Nil + uujl*Nik + uukl*Nij;
        const double Sym6uuN1112 = 3*uu11*N12+3*uu12*N11;
        const double Sym6uuN1122 = uu11*N22 + uu22*N11 + 4*uu12*N12;
        const double Sym6uuN1123 = uu11*N23 + 2*uu12*N13 + 2*uu13*N12 + uu23*N11;
        const double Sym6uuN1222 = 3*uu12*N22 + 3*uu22*N12;
        const double Sym6uuN1223 = 2*uu12*N23 + uu13*N22 + uu22*N13 + 2*uu23*N12;
        const double Sym6uuN1233 = uu12*N33 + 2*uu13*N23 + 2*uu23*N13 + uu33*N12;
        const double Sym6uuN2222 = 6*uu22*N22;
        const double Sym6uuN2223 = 3*uu22*N23 + 3*uu23*N22;
        const double Sym6uuN2233 = uu22*N33 + 4*uu23*N23 + uu33*N22;
        const double Sym6uuN2333 = 3*uu23*N33 + 3*uu33*N23;
        //
        const double RN1112 = R1112 + Sym4uQN1112 - Sym6uuN1112 + 3*rho*uu11*uu12;
        const double RN1122 = R1122 + Sym4uQN1122 - Sym6uuN1122 + 3*rho*uu11*uu22;
        const double RN1123 = R1123 + Sym4uQN1123 - Sym6uuN1123 + 3*rho*uu11*uu23;
        const double RN1222 = R1222 + Sym4uQN1222 - Sym6uuN1222 + 3*rho*uu12*uu22;
        const double RN1223 = R1223 + Sym4uQN1223 - Sym6uuN1223 + 3*rho*uu12*uu23;
        const double RN1233 = R1233 + Sym4uQN1233 - Sym6uuN1233 + 3*rho*uu12*uu33;
        const double RN2222 = R2222 + Sym4uQN2222 - Sym6uuN2222 + 3*rho*uu22*uu22;
        const double RN2223 = R2223 + Sym4uQN2223 - Sym6uuN2223 + 3*rho*uu22*uu23;
        const double RN2233 = R2233 + Sym4uQN2233 - Sym6uuN2233 + 3*rho*uu22*uu33;
        const double RN2333 = R2333 + Sym4uQN2333 - Sym6uuN2333 + 3*rho*uu23*uu33;
    
        
        // 2-component of flux function
        //
        flux.set(j,n_rho,  M2 );
        flux.set(j,n_M1 ,  N12 );
        flux.set(j,n_M2 ,  N22 );
        flux.set(j,n_M3 ,  N23 );
        flux.set(j,n_N11,  QN112);
        flux.set(j,n_N12,  QN122);
        flux.set(j,n_N13,  QN123);
        flux.set(j,n_N22,  QN222);
        flux.set(j,n_N23,  QN223);
        flux.set(j,n_N33,  QN233);
        flux.set(j,n_Q111, RN1112);
        flux.set(j,n_Q112, RN1122);
        flux.set(j,n_Q113, RN1123);
        flux.set(j,n_Q122, RN1222);
        flux.set(j,n_Q123, RN1223);
        flux.set(j,n_Q133, RN1233);
        flux.set(j,n_Q222, RN2222);
        flux.set(j,n_Q223, RN2223);
        flux.set(j,n_Q233, RN2233);
        flux.set(j,n_Q333, RN2333);
    }
}

