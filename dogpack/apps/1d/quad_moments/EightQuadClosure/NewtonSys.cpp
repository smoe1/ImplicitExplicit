#include "dogdefs.h"
#include "dog_math.h"

void NewtonSys(const double rho, double* Qt, double* mut0, double* mutval)

  {
    void Feqns(const double rho, double* Qt, double* v, double* F);

    double TOL = 1.0e-13;
    int MaxIts = 1000;
    int NumIts = 0;
    int mstop = 0;
    double temp[3];
    double det,dis,fnorm;

    double Qt111 = Qt[0];
    double Qt112 = Qt[1];
    double Qt113 = Qt[2];
    double Qt122 = Qt[3];
    double Qt123 = Qt[4];
    double Qt133 = Qt[5];
    double Qt222 = Qt[6];
    double Qt223 = Qt[7];
    double Qt233 = Qt[8];
    double Qt333 = Qt[9];

    mutval[0] = mut0[0];
    mutval[1] = mut0[1];
    mutval[2] = mut0[2];

    double J11=2.0*mutval[0]-Qt111;
    double J12=-Qt112;
    double J13=-Qt113;
    double J21=-Qt122;
    double J22=2.0*mutval[1]-Qt222;
    double J23=-Qt223;
    double J31=-Qt133;
    double J32=-Qt233;
    double J33=2.0*mutval[2]-Qt333;

    double f[3];
    Feqns(rho,Qt,mutval,f);


    fnorm = sqrt(pow(f[0],2)+pow(f[1],2)+pow(f[2],2));

    if (fnorm<TOL)
       {
          mstop = 1;
       }
      
   
    while (mstop == 0)
      {
         NumIts = NumIts +1;

         det = J13*J22*J31-J12*J23*J31-J13*J21*J32+J11*J23*J32+J12*J21*J33-J11*J22*J33;

         temp[0] = mutval[0]-1.0/det*((J23*J32-J22*J33)*f[0]
                   -(J13*J32-J12*J33)*f[1]+(J13*J22-J12*J23)*f[2]);
         temp[1] = mutval[1]-1.0/det*(-(J23*J31-J21*J33)*f[0]
                   +(J13*J31-J11*J33)*f[1]-(J13*J21-J11*J23)*f[2]);
         temp[2] = mutval[2]-1.0/det*((J22*J31-J21*J32)*f[0]
                   -(J12*J31-J11*J32)*f[1]+(J12*J21-J11*J22)*f[2]);

         dis = sqrt(pow((temp[0]-mutval[0]),2)+pow((temp[1]-mutval[1]),2)+pow((temp[2]-mutval[2]),2));

         Feqns(rho,Qt,temp,f);
 
         fnorm = sqrt(pow(f[0],2)+pow(f[1],2)+pow(f[2],2));

         
         mutval[0] = temp[0];
         mutval[1] = temp[1];
         mutval[2] = temp[2];

         if(dis<=TOL && fnorm<=TOL)
            {
              mstop = 1;
            }
         else
            {
              if (NumIts>=MaxIts)
                {
                   printf("Does not converge\n");
                   printf("\n");
                   exit(1);
                }
              J11 = 2.0*mutval[0]-Qt111;
              J22 = 2.0*mutval[1]-Qt222;
              J33 = 2.0*mutval[2]-Qt333;
            }
       }
}



void Feqns(const double rho, double* Qt, double* v, double* F)
{
    const double Qt111 = Qt[0];
    const double Qt112 = Qt[1];
    const double Qt113 = Qt[2];
    const double Qt122 = Qt[3];
    const double Qt123 = Qt[4];
    const double Qt133 = Qt[5];
    const double Qt222 = Qt[6];
    const double Qt223 = Qt[7];
    const double Qt233 = Qt[8];
    const double Qt333 = Qt[9];
    
    double v1 = v[0];
    double v2 = v[1];
    double v3 = v[2];
    
    F[0] = pow(v1,2)-1.0/rho-Qt111*v1-Qt112*v2-Qt113*v3;
    F[1] = pow(v2,2)-1.0/rho-Qt122*v1-Qt222*v2-Qt233*v3;
    F[2] = pow(v3,2)-1.0/rho-Qt133*v1-Qt233*v2-Qt333*v3;

}