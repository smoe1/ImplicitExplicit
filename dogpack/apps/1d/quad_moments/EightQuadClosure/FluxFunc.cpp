#include "dogdefs.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Qin,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();
  
  for (int i=1; i<=numpts; i++)
    {
      const double M000 = Qin.get(i,1);
      const double M100 = Qin.get(i,2);
      const double M010 = Qin.get(i,3);
      const double M001 = Qin.get(i,4);
      const double M200 = Qin.get(i,5);
      const double M110 = Qin.get(i,6);
      const double M101 = Qin.get(i,7);
      const double M020 = Qin.get(i,8);
      const double M011 = Qin.get(i,9);
      const double M002 = Qin.get(i,10);
      const double M300 = Qin.get(i,11);
      const double M210 = Qin.get(i,12);
      const double M201 = Qin.get(i,13);
      const double M120 = Qin.get(i,14);
      const double M111 = Qin.get(i,15);
      const double M102 = Qin.get(i,16);
      const double M030 = Qin.get(i,17);
      const double M021 = Qin.get(i,18);
      const double M012 = Qin.get(i,19);
      const double M003 = Qin.get(i,20);

      const double rho  = M000;
      const double u1   = M100/rho;
      const double u2   = M010/rho;
      const double u3   = M001/rho;
      const double P11  = M200 - rho*u1*u1;
      const double P22  = M020 - rho*u2*u2;
      const double P33  = M002 - rho*u3*u3;
      const double P12  = M110 - rho*u1*u2;
      const double P13  = M101 - rho*u1*u3;
      const double P23  = M011 - rho*u2*u3;
      const double Q111 = M300 - 3.0*u1*P11 - rho*u1*u1*u1;
      const double Q112 = M210 - 2.0*u1*P12 - u2*P11 - rho*u1*u1*u2;
      const double Q113 = M201 - 2.0*u1*P13 - u3*P11 - rho*u1*u1*u3;
      const double Q122 = M120 - 2.0*u2*P12 - u1*P22 - rho*u1*u2*u2;
      const double Q123 = M111 - u1*P23 - u2*P13 - u3*P12 - rho*u1*u2*u3;
      const double Q133 = M102 - 2.0*u3*P13 - u1*P33 - rho*u1*u3*u3;
      const double Q222 = M030 - 3.0*u2*P22 - rho*u2*u2*u2;
      const double Q223 = M021 - 2.0*u2*P23 - u3*P22 - rho*u2*u2*u3;
      const double Q233 = M012 - 2.0*u3*P23 - u2*P33 - rho*u2*u3*u3;
      const double Q333 = M003 - 3.0*u3*P33 - rho*u3*u3*u3;


      const double detP= P11*P22*P33-P11*pow(P23,2)-P22*pow(P13,2)-P33*pow(P12,2)+2.0*P12*P13*P23;


      if ( rho<=0.0 || P11<=0.0 || P22<=0.0 || P33<=0.0 || P11*P22-pow(P12,2)<=0.0 || detP<=0.0)
        {
          printf("              rho = %22.16e \n",rho);
          printf("              P11 = %22.16e \n",P11);
          printf("              P22 = %22.16e \n",P22);
          printf("              P33 = %22.16e \n",P33);
          printf("              P12 = %22.16e \n",P12);
          printf("              P13 = %22.16e \n",P13);
          printf("              P23 = %22.16e \n",P23);
          printf("  P11*P22-P12*P12 = %22.16e \n",P11*P22-pow(P12,2));
          printf("             detP = %22.16e \n",detP);
          printf("\n");
          exit(1);
        }

      //Decomposition of P tensor
      const double L11 = sqrt(P11);
      const double L22 = sqrt(P22-P12*P12/P11);
      const double L33 = sqrt(P33-P13*P13/P11-pow((P11*P23-P13*P13),2)/(P11*P11*P22-P11*P12*P12));
      const double L21 = P12/sqrt(P11);
      const double L31 = P13/sqrt(P11);
      const double L32 = (P11*P23-P12*P13)/sqrt(P11*P11*P22-P11*P12*P12);


      const double u[3]   = {u1,u2,u3};
      const double L[6]   = {L11,L22,L33,L21,L31,L32};
      const double Q[10]  = {Q111,Q112,Q113,Q122,Q123,Q133,Q222,Q223,Q233,Q333};
      double Qt[10];
      

      void LinTrans(const double* Q, const double* L, double* Qt);
      LinTrans(Q,L,Qt);

/*
    
    printf(" u1 = %22.16e \n",u1);
    printf(" u2 = %22.16e \n",u2);
    printf(" u3 = %22.16e \n",u3);

    printf(" P11 = %22.16e \n",P11);
    printf(" P22 = %22.16e \n",P22);
    printf(" P33 = %22.16e \n",P33);
    printf(" P12 = %22.16e \n",P12);
    printf(" P13 = %22.16e \n",P13); 
    printf(" P23 = %22.16e \n",P23);
    

    printf(" L11 = %22.16e \n",L11);
    printf(" L22 = %22.16e \n",L22);
    printf(" L33 = %22.16e \n",L33);
    printf(" L21 = %22.16e \n",L21);
    printf(" L31 = %22.16e \n",L31); 
    printf(" L32 = %22.16e \n",L32);    
  
    printf(" Qt111 = %22.16e \n",Qt[0]);
    printf(" Qt112 = %22.16e \n",Qt[1]);
    printf(" Qt113 = %22.16e \n",Qt[2]);
    printf(" Qt122 = %22.16e \n",Qt[3]);
    printf(" Qt123 = %22.16e \n",Qt[4]); 
    printf(" Qt133 = %22.16e \n",Qt[5]);
    printf(" Qt222 = %22.16e \n",Qt[6]);
    printf(" Qt223 = %22.16e \n",Qt[7]);  
    printf(" Qt233 = %22.16e \n",Qt[8]); 
    printf(" Qt333 = %22.16e \n",Qt[9]);
*/
      
      double mut[8][3];
      double mu[8][3];
      double w[8];

      void GetMus(const double rho, const double* u,
                  const double* L, double* Qt,
                  double mut[][3], double mu[][3]);
      GetMus(rho,u,L,Qt,mut,mu);
/*
      for(int i=0;i<8;i++)
       {
        for(int j=0;j<3;j++)
          {
             printf("  mut[%d][%d] = %f \t",i+1,j+1,mut[i][j]);
          }
        printf("\n");
       }
*/

      void GetWeights(const double rho, double Qt123, double mut[][3], double* w);
      GetWeights(rho,Qt[4],mut,w);
/*
      for(int i=0;i<8;i++)
       {
            printf(" w[%d] = %f \n",i+1,w[i]);
       }
      exit(1);
*/
      double M400=0.0;
      double M310=0.0;
      double M301=0.0;
      double M220=0.0;
      double M211=0.0;
      double M202=0.0;
      double M130=0.0;
      double M121=0.0;
      double M112=0.0;
      double M103=0.0;

      for(int j=0;j<8;j++)
        {
          M400 += w[j]*pow(mu[j][0],4);
          M310 += w[j]*pow(mu[j][0],3)*mu[j][1];
          M301 += w[j]*pow(mu[j][0],3)*mu[j][2];
          M220 += w[j]*pow(mu[j][0],2)*pow(mu[j][1],2);
          M211 += w[j]*pow(mu[j][0],2)*mu[j][1]*mu[j][2];
          M202 += w[j]*pow(mu[j][0],2)*pow(mu[j][2],2);
          M130 += w[j]*mu[j][0]*pow(mu[j][1],3);
          M121 += w[j]*mu[j][0]*pow(mu[j][1],2)*mu[j][2];
          M112 += w[j]*mu[j][0]*mu[j][1]*pow(mu[j][2],2);
          M103 += w[j]*mu[j][0]*pow(mu[j][2],3);
        }

      

      flux.set(i,1, M100 );
      flux.set(i,2, M200 );
      flux.set(i,3, M110 ); 
      flux.set(i,4, M101 );
      flux.set(i,5, M300 );
      flux.set(i,6, M210 );
      flux.set(i,7, M201 );
      flux.set(i,8, M120 );
      flux.set(i,9, M111 );
      flux.set(i,10,M102 );
      flux.set(i,11,M400 );
      flux.set(i,12,M310 );
      flux.set(i,13,M301 );
      flux.set(i,14,M220 );
      flux.set(i,15,M211 );
      flux.set(i,16,M202 );
      flux.set(i,17,M130 );
      flux.set(i,18,M121 );
      flux.set(i,19,M112 );
      flux.set(i,20,M103 );
    }

}
