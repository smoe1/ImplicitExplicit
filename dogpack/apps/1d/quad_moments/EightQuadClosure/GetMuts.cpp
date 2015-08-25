#include "dogdefs.h"
#include "dog_math.h"

void GetMus(const double rho, const double* u,
            const double* L, double* Qt,
            double mut[][3], double mu[][3])
{
   void NewtonSys(const double rho, double* Qt, double* mut0, double* mutval);

   const double u1    = u[0];
   const double u2    = u[1];
   const double u3    = u[2];
   const double L11   = L[0];
   const double L22   = L[1];
   const double L33   = L[2];
   const double L21   = L[3];
   const double L31   = L[4];
   const double L32   = L[5];


   double mut0[3];
   double mutval[3]; 
   double n = pow(rho,-0.5);


   double mutguess[8][3] = {{n,n,n},{n,n,-n},{n,-n,n},{n,-n,-n},
                           {-n,n,n},{-n,-n,n},{-n,n,-n},{-n,-n,-n}}; 


   for (int i=0;i<8;i++)
     {
        for (int j=0;j<3;j++)
          {
             mut0[j] = mutguess[i][j];
          }
   
        NewtonSys(rho,Qt,mut0,mutval);
       
        for (int j=0;j<3;j++)
          {
             mut[i][j] = mutval[j];
            // printf(" mut[%d][%d] = %f  \n",i,j,mut[i][j]);
          }

        mu[i][0]  = L11*mutval[0]+u1;
        mu[i][1]  = L21*mutval[0]+L22*mutval[1]+u2;
        mu[i][2]  = L31*mutval[0]+L32*mutval[1]+L33*mutval[2]+u3;
     }



}
